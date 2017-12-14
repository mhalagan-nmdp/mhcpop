# mhcpop
Population genetics scripts for the MHC

## Setup

```
virtualenv -p python3 venv
source venv/bin/activate
pip install -r requirements.txt
```

## Example

```
python bin/LD.py -f /path/to/frequencies -o /path/to/output/ -t ExampleTitle -v

12/14/2017 02:57:49 PM - root - INFO - Loading /path/to/frequencies/ExamplePop1.freqs...
12/14/2017 02:57:49 PM - root - INFO - Finished loading /path/to/frequencies/ExamplePop1.freqs
12/14/2017 02:57:49 PM - root - INFO - 25629 unique haplotypes
12/14/2017 02:57:54 PM - root - INFO - A~C       : # Haplos =   1212, Min D = -0.0102, Max D = 0.0231, # Significant =    759
12/14/2017 02:57:54 PM - root - INFO - C~B       : # Haplos =    351, Min D = -0.0112, Max D = 0.0731, # Significant =    322
12/14/2017 02:57:57 PM - root - INFO - B~DRB1    : # Haplos =   3456, Min D = -0.0046, Max D = 0.0259, # Significant =   2154
12/14/2017 02:57:58 PM - root - INFO - DRB1~DQB1 : # Haplos =    288, Min D = -0.0210, Max D = 0.0895, # Significant =    256
12/14/2017 02:57:58 PM - root - INFO - Finished ExamplePop1
12/14/2017 02:57:58 PM - root - INFO - ---------------------------------------------------------------------------------------
12/14/2017 02:57:58 PM - root - INFO - Loading /path/to/frequencies/ExamplePop2.freqs...
12/14/2017 02:57:58 PM - root - INFO - Finished loading /path/to/frequencies/ExamplePop2.freqs
12/14/2017 02:57:58 PM - root - INFO - 41639 unique haplotypes
12/14/2017 02:58:07 PM - root - INFO - A~C       : # Haplos =   1465, Min D = -0.0103, Max D = 0.0185, # Significant =    865
12/14/2017 02:58:07 PM - root - INFO - C~B       : # Haplos =    486, Min D = -0.0124, Max D = 0.0567, # Significant =    419
12/14/2017 02:58:10 PM - root - INFO - B~DRB1    : # Haplos =   3993, Min D = -0.0048, Max D = 0.0212, # Significant =   2681
12/14/2017 02:58:11 PM - root - INFO - DRB1~DQB1 : # Haplos =    417, Min D = -0.0238, Max D = 0.0884, # Significant =    356
12/14/2017 02:58:11 PM - root - INFO - Finished ExamplePop2
12/14/2017 02:58:11 PM - root - INFO - ---------------------------------------------------------------------------------------
12/14/2017 02:59:17 PM - root - INFO - FINISHED ANALAYSIS
12/14/2017 02:59:17 PM - root - INFO - Writing to /path/to/output/ExampleTitle.xlsx...

```

