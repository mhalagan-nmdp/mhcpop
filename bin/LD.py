
from scipy.stats import chi2
import pandas as pd
import argparse

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
                        help="frequency file (ex. AAFA.freqs)",
                        type=str)

    parser.add_argument("-o", "--outdir",
                        required=True,
                        help="Ouput directory",
                        type=str)

    parser.add_argument("-v", "--verbose",
                        help="Option for running in verbose",
                        default=False,
                        type=bool)

    args = parser.parse_args()
    data_file = args.freq
    outdir = args.outdir
    df = pd.read_csv(data_file)
    freq_name = data_file.split("/")[len(data_file.split("/"))-1].split(".")[0]
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

    if 'A' in loci and 'B' in loci and ['A', 'B'] not in loci_combos:
        loci_combos.append(['A', 'B'])

    for loc_comb in loci_combos:
        a, b = loc_comb[0], loc_comb[1]
        ab_hap = "~".join(loc_comb)
        a_sum = sum_dic[a]
        b_sum = sum_dic[b]
        ab_sum = freq_df.groupby(loc_comb, as_index=False).sum()
        ab_df = ab_sum.apply(lambda row: pd.Series(["~".join([row[a], row[b]]),
                                                    row['Count']]), axis=1)\
                      .rename(index=str,
                              columns={0: ab_hap, 1: "Count"})
        ab = ab_df.groupby([ab_hap], as_index=True)\
                  .sum().to_dict(orient="index")
        ab_sum['D'] = ab_sum.apply(lambda row:
                                   row['Freq'] -
                                   (a_sum[row[a]]['Freq']*b_sum[row[b]]['Freq']),
                                   axis=1)
        ab_dict = ab_sum.groupby([a])[b].apply(list).to_dict()
        ba_dict = ab_sum.groupby([b])[a].apply(list).to_dict()
        total_count = freq_df['Count'].sum()
        ab_df['p_val'] = ab_df.apply(lambda row: get_pval(row, ab_hap, ab,
                                                          ab_dict,
                                                          ba_dict,
                                                          total_count),
                                     axis=1)
        ab_sum = ab_sum.apply(lambda row: pd.Series(["~".join([row[a],
                                                               row[b]]),
                                                    row['Freq'],
                                                    row['Count'],
                                                    row['D']]), axis=1)\
                       .rename(index=str,
                               columns={0: ab_hap, 1: "Freq", 2: "Count",
                                        3: "D"})
        ab_df['significant'] = ab_df.apply(lambda row:
                                           '*' if row['p_val'] < .05 else '',
                                           axis=1)
        tmp = pd.merge(ab_df, ab_sum, on=ab_hap)
        final_df = tmp.filter([ab_hap, 'Freq', 'Count_x', 'D', 'p_val',
                               'significant'],
                              axis=1)
        final_df = final_df.rename(index=str,
                                   columns={"Count_x": "Count"})
        outfile = outdir + "/" + freq_name + "_" + "".join(loc_comb) + ".csv"
        final_df.to_csv(outfile, index=False)


if __name__ == '__main__':
    """The following will be run if file is executed directly,
    but not if imported as a module"""
    main()


