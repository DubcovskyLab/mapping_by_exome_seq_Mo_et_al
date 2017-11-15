# mapping_by_exome_seq_Mo_et_al

## BSA Simulation
This script takes two files, one with names for individuals in the high bulk and one with names for the low bulk. Both should be sorted is descending order (phenotypic ranking) and the script will remove the bottom individual from the high bulk and the top individual from the low bulk for each iteration.

usage: step_simulation.py [-h] -m MPILEUP -1 HIGHBULK -2 LOWBULK [-s SNPS_OF_INTEREST] [-r REPORT_THRESHOLD] [-d DISCARD_THRESHOLD] [-o OUT_PREFIX]

optional arguments: -h, --help show this help message and exit -s SNPS_OF_INTEREST, --SNPs_of_interest SNPS_OF_INTEREST SNP or SNPs to report rankings for (Comma separated for multiple SNPs). Format CONTIG:POSITION -r REPORT_THRESHOLD, --report_threshold REPORT_THRESHOLD Threshold proportion for |SNP-index| to output. Default reports all |SNP-index| values [0.0] -d DISCARD_THRESHOLD, --discard_threshold DISCARD_THRESHOLD Threshold to discard SNPS with extreme values in both bulks. Value is distance from 1 or 0 [0.2] -o OUT_PREFIX, --out_prefix OUT_PREFIX Prefix for the output files [mpileup_file_name]

required named arguments: -m MPILEUP, --mpileup MPILEUP Parsed mpileup (run through generateCountsOfBasesForFilteredMpilup.py) -1 HIGHBULK, --highBulk HIGHBULK List of individuals in bulk group 1, "High bulk". Sorted in decending order (highest at top) -2 LOWBULK, --lowBulk LOWBULK List of individuals in bulk group 2, "Low bulk". Sorted in decending order (lowest at bottom)

## generate Counts from MAPS output
This script takes in a parental mutation file (--m) which can be downloaded from the Dubcovsky Lab's TILLING project and a filtered mpileup file (--p) for all mapped samples for the called mutations in the parental line. The script will then print a TSV file with counted NTs at each mutation position and called SNPs.  

usage: generateCountsOfBasesForFilteredMpilup.py [-h] -m MUTATIONFILE -p PILEUPFILE [optional arguments]

Optional Arguments
  -d LOWDEPTH, --depth=LOWDEPTH
                        Minimum depth to be considered in the analysis;
                        Default 5 reads
  -o HOMPRCNT, --hom=HOMPRCNT
                        Minimum percentage to be considered homozygous;
                        Default 0.9
  -u UPPERHET, --upper_het=UPPERHET
                        Upper limit to be considered heterozygous; Default 0.7
  -l LOWERHET, --lower_het=LOWERHET
                        Lower limit to be considered heterozygous; Default 0.3
