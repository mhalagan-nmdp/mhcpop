
from itertools import repeat
from scipy.stats import chi2
import pandas as pd
import argparse
import logging
import glob
import sys
import os

logging.basicConfig(format='%(asctime)s - %(name)s - %(levelname)s - %(message)s',
                    datefmt='%m/%d/%Y %I:%M:%S %p',
                    level=logging.DEBUG)

loci_order = {"A": 0, "C": 1, "B": 2, "DRB1": 3, "DRBX": 4,
              "DQA1": 5, "DQB1": 6, "DPA1": 7, "DPB1": 8}


def get_pval(row, hap, ab, ab_dict, ba_dict, total_count):
    ab_hap = row[hap]
    a_allele, b_allele = ab_hap.split("~")
    total_n = ab[ab_hap]['Count']
    total_wo = total_count - total_n
    b_total = sum(map(lambda j: ab["~".join([a_allele, j])]['Count'],
                  filter(lambda i: i != b_allele, ab_dict[a_allele])))
    a_total = sum(map(lambda j: ab["~".join([j, b_allele])]['Count'],
                  filter(lambda i: i != a_allele, ba_dict[b_allele])))
    chisq_stat = ((((total_n*total_wo) - (b_total*a_total))**2) *
                 (total_n + b_total+a_total+total_wo)) / \
                 ((total_n+a_total) * (b_total+total_wo) *
                  (a_total + total_wo) * (total_n + b_total))
    return 1 - chi2.cdf(x=chisq_stat, df=1)


def main():
    """This is run if file is directly executed, but not if imported as
    module. Having this in a separate function  allows importing the file
    into interactive python, and still able to execute the
    function for testing"""
    parser = argparse.ArgumentParser()
    parser.add_argument("-f", "--freq",
                        required=True,
                        help="frequency directory",
                        type=str)

    parser.add_argument("-o", "--outdir",
                        required=True,
                        help="Ouput directory",
                        type=str)

    parser.add_argument("-t", "--title",
                        required=True,
                        help="Title of excel file",
                        type=str)

    parser.add_argument("-v", "--verbose",
                        help="Option for running in verbose",
                        action='store_true')

    args = parser.parse_args()

    verbose = False
    if args.verbose:
        verbose = True

    freqdir = args.freq
    outdir = args.outdir
    xlsx_name = outdir + "/" + args.title + ".xlsx"
    freq_glob = freqdir + "/*.freqs"

    # Check that frequency directory exists
    if not os.path.isdir(freqdir):
        logging.warning(freqdir + " is not a valid directory!")
        parser.print_help()
        sys.exit(2)

    # Check that output directory exists
    if not os.path.isdir(outdir):
        logging.warning(outdir + " is not a valid directory!")
        parser.print_help()
        sys.exit(2)

    # Check that there are frequency files
    if len(glob.glob(freq_glob)) < 1:
        logging.warning("No frequency file found in " + freqdir)
        parser.print_help()
        sys.exit(2)

    # Create excel file
    writer = pd.ExcelWriter(xlsx_name, engine='xlsxwriter')

    # Loop through frequency files
    for f in glob.glob(freq_glob):

        if verbose:
            logging.info("Loading " + f + "...")

        # Read freq file into pandas df
        df = pd.read_csv(f)

        if verbose:
            logging.info("Finished loading " + f)
            hapcnt = df['Haplo'].count()
            logging.info(str(hapcnt) + " unique haplotypes")

        # get population name
        freq_name = f.split("/")[len(f.split("/"))-1].split(".")[0]

        # get loci in frequency file
        loci_df = df['Haplo'].apply(lambda x: pd.Series(x.split('~')))
        loci = list(map(lambda x: x.split("*")[0], df['Haplo'][0].split("~")))

        loci_dic = {i: x for i, x in enumerate(loci)}
        freq_df = loci_df.filter(loci_df.columns.values, axis=1)

        freq_df['Freq'], freq_df['Count'] = df['Freq'], df['Count']
        freq_df = freq_df.rename(index=str,
                                 columns=loci_dic)

        sum_dic = dict(zip(loci, map(lambda loc:
                                     freq_df.groupby([loc], as_index=True).sum()
                                     .to_dict(orient="index"), loci)))

        sorted_loci = sorted(loci, key=lambda x: loci_order[x])
        loci_combos = list(map(lambda x: [sorted_loci[x], sorted_loci[x+1]],
                               range(0, len(sorted_loci)-1)))

        for loc_comb in loci_combos:
            a, b = loc_comb[0], loc_comb[1]
            ab_hap = "~".join(loc_comb)

            a_sum = sum_dic[a]
            b_sum = sum_dic[b]

            # get all 2-loc haplotypes
            ab_sum = freq_df.groupby(loc_comb, as_index=False).sum()
            ab_df = ab_sum.apply(lambda row: pd.Series(["~".join([row[a],
                                                        row[b]]),
                                                        row['Count']]),
                                 axis=1)\
                          .rename(index=str,
                                  columns={0: ab_hap, 1: "Count"})

            # sum the counts across all 2-loc haplotypes
            ab = ab_df.groupby([ab_hap], as_index=True)\
                      .sum().to_dict(orient="index")

            # Calculate LD
            ab_sum['D'] = ab_sum.apply(lambda row:
                                       row['Freq'] -
                                       (a_sum[row[a]]['Freq']*b_sum[row[b]]['Freq']),
                                       axis=1)
            dmin = min(ab_sum['D'])
            dmax = max(ab_sum['D'])
            nhaps = ab_sum['D'].count()

            # Calculate LD prime (D')
            ab_sum['D\''] = ab_sum.apply(lambda row:
                                            row['D'] / dmax if row['D'] > 0 else row['D'] / dmin,
                                            axis=1)
            ab_dict = ab_sum.groupby([a])[b].apply(list).to_dict()
            ba_dict = ab_sum.groupby([b])[a].apply(list).to_dict()
            total_count = freq_df['Count'].sum()

            # Calculate p_val of LD
            ab_df['p_val'] = ab_df.apply(lambda row: get_pval(row, ab_hap, ab,
                                                              ab_dict,
                                                              ba_dict,
                                                              total_count),
                                         axis=1)
            ab_sum = ab_sum.apply(lambda row: pd.Series(["~".join([row[a],
                                                                   row[b]]),
                                                         row['Freq'],
                                                         row['Count'],
                                                         row['D'],
                                                         row['D\'']]),
                                  axis=1)\
                           .rename(index=str,
                                   columns={0: ab_hap, 1: "Freq", 2: "Count",
                                            3: "D", 4: "D\'"})
            ab_df['significant'] = ab_df.apply(lambda row:
                                               '*' if row['p_val'] < .05 else '',
                                               axis=1)
            sign_cnt = ab_df[ab_df['significant'] == "*"]['significant']\
                .count()

            if verbose:
                log ="{:10}: # Haplos = {:6d}, Min D = {:.4f}, Max D = {:.4f}, # Significant = {:6d}".format(ab_hap, nhaps, dmin, dmax, sign_cnt)
                logging.info(log)

            tmp = pd.merge(ab_df, ab_sum, on=ab_hap)
            final_df = tmp.filter([ab_hap, 'Freq', 'Count_x', 'D',
                                   'D\'', 'p_val',
                                  'significant'],
                                  axis=1)
            final_df = final_df.rename(index=str,
                                       columns={"Count_x": "Count"})
            sheet = freq_name + "_LD-" + "".join(loc_comb)
            final_df.to_excel(writer, sheet_name=sheet, index=False)

        if verbose:
            logging.info("Finished " + freq_name)
            logging.info("".join(repeat("-", 110)))
    if verbose:
        logging.info("FINISHED ANALAYSIS")
        logging.info("Writing to " + xlsx_name)
    writer.save()

if __name__ == '__main__':
    """The following will be run if file is executed directly,
    but not if imported as a module"""
    main()


