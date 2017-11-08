# mapping_by_exome_seq_Mo_et_al

## BSA Simulation
This script takes two files, one with names for individuals in the high bulk and one with names for the low bulk. Both should be sorted is descending order (phenotypic ranking) and the script will remove the bottom individual from the high bulk and the top individual from the low bulk for each iteration.

usage: step_simulation.py [-h] -m MPILEUP -1 HIGHBULK -2 LOWBULK [-s SNPS_OF_INTEREST] [-r REPORT_THRESHOLD] [-d DISCARD_THRESHOLD] [-o OUT_PREFIX]

optional arguments: -h, --help show this help message and exit -s SNPS_OF_INTEREST, --SNPs_of_interest SNPS_OF_INTEREST SNP or SNPs to report rankings for (Comma separated for multiple SNPs). Format CONTIG:POSITION -r REPORT_THRESHOLD, --report_threshold REPORT_THRESHOLD Threshold proportion for |SNP-index| to output. Default reports all |SNP-index| values [0.0] -d DISCARD_THRESHOLD, --discard_threshold DISCARD_THRESHOLD Threshold to discard SNPS with extreme values in both bulks. Value is distance from 1 or 0 [0.2] -o OUT_PREFIX, --out_prefix OUT_PREFIX Prefix for the output files [mpileup_file_name]

required named arguments: -m MPILEUP, --mpileup MPILEUP Parsed mpileup (run through generateCountsOfBasesForFilteredMpilup.py) -1 HIGHBULK, --highBulk HIGHBULK List of individuals in bulk group 1, "High bulk". Sorted in decending order (highest at top) -2 LOWBULK, --lowBulk LOWBULK List of individuals in bulk group 2, "Low bulk". Sorted in decending order (lowest at bottom)
