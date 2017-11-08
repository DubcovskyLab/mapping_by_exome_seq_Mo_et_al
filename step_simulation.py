#!/usr/bin/env python
import argparse
import sys
from subprocess import call
import os
import operator
from scipy.stats import rankdata

# Written by Tyson Howell at the University of California, Davis July 2017.
# All information obtained/inferred with this script is without any
# implied warranty of fitness for any purpose or use whatsoever.
# Use at your own risk.

parser = argparse.ArgumentParser()
requiredNamed = parser.add_argument_group('required named arguments')
requiredNamed.add_argument("-m", "--mpileup", required=True, help="Parsed mpileup (run through generateCountsOfBasesForFilteredMpilup.py)")
requiredNamed.add_argument("-1", "--highBulk", required=True, help="List of individuals in bulk group 1, \"High bulk\". Sorted in decending order (highest at top)")
requiredNamed.add_argument("-2", "--lowBulk", required=True, help="List of individuals in bulk group 2, \"Low bulk\". Sorted in decending order (lowest at bottom)")
parser.add_argument("-s", "--SNPs_of_interest", default=None, help="SNP or SNPs to report rankings for (Comma separated for multiple SNPs). Format CONTIG:POSITION")
parser.add_argument("-r", "--report_threshold", default=0.0, help="Threshold proportion for |SNP-index| to output. Default reports all |SNP-index| values [0.0]")
parser.add_argument("-d", "--discard_threshold", default=0.2, help="Threshold to discard SNPS with extreme values in both bulks. Value is distance from 1 or 0 [0.2]")
parser.add_argument("-o", "--out_prefix", help="Prefix for the output files [mpileup_file_name]")
parser.add_argument("-q", "--quiet", default=False, action="store_true", help="Turn off extra messages to STDOUT when program is running (current iteration number and header)")
parser.add_argument("-R", "--ranks", default=False, action="store_true", help="Report competitive rankings in output")
args = parser.parse_args()

#Set the output file prefix. If no argument given, set it to the input without the extension
out = args.out_prefix if args.out_prefix else os.path.plitext(args.mpileup)[0]

with open(args.highBulk) as f:
    high_length = sum(1 for _ in f)

with open(args.lowBulk) as f:
    low_length = sum(1 for _ in f)

#print "High bulk length = " + str(high_length)

if (high_length != low_length):
    print "Bulks 1 and 2 must be of equal size!"
    print "Bulk 1 length = " + str(high_length)
    print "Bulk 2 length = " + str(low_length)
    sys.exit
else:
    #print "Both bulks of size " + str(high_length)
    None

bulk1 = []
bulk2 = []
for line in open (args.highBulk):
    line = line.rstrip()
    bulk1.append(line)

for line in open (args.lowBulk):
    line = line.rstrip()
    bulk2.append(line)

#print bulk1
#print bulk2

mpileup = open (args.mpileup)

header = mpileup.readline().split()
#print header

#First 6 fields (0-5) are always the same, individuals come after
individuals = header[6:]
#print individuals

#Remove callsnp columns from data
#baseind = [x for x in individuals if "callsnp" not in x]
#print baseind

counts = {"dictindex":[]}
for el in individuals:
    counts[el] = []
    counts["dictindex"].append(el)

positions = []
#print counts
for line in mpileup:
    line = line.rstrip ()
    position = line.split("\t")
    #some files have extra spaces; remove them
    position = [x.strip(' ') for x in position]
    positions.append(position[0:5])
    #print positions
    index = 6
    #Add count data for each individual to the correct dictionary entry
    for el in counts["dictindex"[:]]:
        counts[el].append(position[index])
        index += 1

#print positions
#print counts

#Calculate coverage and SNP-index for each SNP of selected individuals
current_num = high_length
current_start = 0
current_pos = 0
high_index = {"dictindex":[]}
low_index = {"dictindex":[]}
deltasnp = {}
for num in range(current_num):
    if not args.quiet: print "Currently calculating for %s individuals in each bulk: " % (current_num)
    for position in positions:
        mutbase = position[4]
        wtbase = position[3]
        if mutbase == "*": 
            #print "Skipping " + str(position)
            continue
        snppos = position[0] + ":" + position[1]
        high_index["dictindex"].append(snppos)
        low_index["dictindex"].append(snppos)
        mutcounthigh = 0
        wtcounthigh = 0
        totcounthigh = 0
        mutcountlow = 0
        wtcountlow = 0
        totcountlow = 0
        #print "Current position in file: " + str(current_pos)
    
        #Calculate SNP-index for high bulk
        for plant in bulk1[current_start:]:
            #print "\"" + plant + "\""
            #Calculate mutant base coverage
            plantmutbase = plant + "-" + mutbase
            plantwtbase = plant + "-" + wtbase
            mutcounthigh = mutcounthigh + int(counts[plantmutbase][current_pos])
            wtcounthigh = wtcounthigh + int(counts[plantwtbase][current_pos])
            ##Calculate total position coverage
            ##Removed, Jun uses just the mutant and WT base, not all bases!
            #planta = plant + "-A"
            #plantc = plant + "-C"
            #plantg = plant + "-G"
            #plantt = plant + "-T"
            #totcounthigh = totcounthigh + int(counts[planta][current_pos]) + int(counts[plantc][current_pos]) + int(counts[plantg][current_pos]) + int(counts[plantt][current_pos])
            totcounthigh = mutcounthigh + wtcounthigh
            #print "%s: mut:%s WT:%s totcounthigh = %s + %s" % (snppos, plantmutbase, plantwtbase, mutcounthigh, wtcounthigh)
            #print totcounthigh
        
        #Calculate SNP-index for low bulk
        for plant in bulk2[:current_num]:
            #print plant
            #Calculate mutant base coverage
            plantmutbase = plant + "-" + mutbase
            plantwtbase = plant + "-" + wtbase
            mutcountlow = mutcountlow + int(counts[plantmutbase][current_pos])
            wtcountlow = wtcountlow + int(counts[plantwtbase][current_pos])

            ##Calculate total position coverage
            #planta = plant + "-A"
            #plantc = plant + "-C"
            #plantg = plant + "-G"
            #plantt = plant + "-T"
            #totcountlow = totcountlow + int(counts[planta][current_pos]) + int(counts[plantc][current_pos]) + int(counts[plantg][current_pos]) + int(counts[plantt][current_pos])
            totcountlow = mutcountlow + wtcountlow
            #print "%s: mut:%s WT:%s totcountlow = %s + %s" % (snppos, plantmutbase, plantwtbase, mutcountlow, wtcountlow)
            #print totcountlow
        #print totcountlow 

        #print snppos + ": totcounthigh: " + str(totcounthigh) + ", totcountlow: " + str(totcountlow)
        if (totcounthigh > 0 and totcountlow > 0):
            high_index[snppos] = float(mutcounthigh) / float(totcounthigh)
            low_index[snppos] = float(mutcountlow) / float(totcountlow)
            #print "%s: Low index: %s, High index: %s" % (snppos, low_index[snppos], high_index[snppos])
            if ((high_index[snppos] > (1 - args.discard_threshold) and low_index[snppos] > (1 - args.discard_threshold)) or (high_index[snppos] < args.discard_threshold and low_index[snppos] < args.discard_threshold)):
                #print "deleted " + snppos + ", high_index = " + str(high_index[snppos]) + " low_index = " + str(low_index[snppos])
                del high_index[snppos]
                del low_index[snppos]
                high_index["dictindex"].remove(snppos)
                low_index["dictindex"].remove(snppos)
    
            #elif mutbase == "*":
            #    print "deleted " + snppos + ", mutbase = " + mutbase
            #    del high_index[snppos]
            #    del low_index[snppos]
            #    del high_index["dictindex"][-1]
            #    del low_index["dictindex"][-1]
    
            else:
                #print "Kept " + snppos
                None
        else: #Remove positions that have zero coverage in both bulks from the dictindex
            high_index["dictindex"].remove(snppos)
            low_index["dictindex"].remove(snppos) 
    
        #print high_index
        #print low_index
        if snppos in high_index: 
            #print "calculation for " + snppos + ": abs(" + str(high_index[snppos]) + " - " + str(low_index[snppos]) + ")"
            #print "snppos " + snppos + " in high_index: " + str(high_index[snppos])
            deltasnp[snppos] = abs(high_index[snppos] - low_index[snppos])
        #print snppos + " " + str(deltasnp)
        current_pos += 1

    #Sort the DeltaSNP-index values by value to get the most significant
    sdeltasnp = sorted(deltasnp.items(), key=operator.itemgetter(1), reverse=True)
            
    if args.SNPs_of_interest:
        snps = args.SNPs_of_interest.split(",")
        #print [i[1] for i in sdeltasnp]
        rankings = len(sdeltasnp) - rankdata([i[1] for i in sdeltasnp], method = 'min') + 1
        #print "Here are the ranks:"
        #for i in range(len(sdeltasnp)):
        #    print sdeltasnp[i][0], rankings[i]
        if not args.quiet: print "SNP:position\t|SNP-index|\tRank"
        for snp in snps:
            if snp in deltasnp:
                #Calculate rank of SNP in sorted list of tupes, add 1 to make it one-based
                index = [item[0] for item in sdeltasnp].index(snp)
                print "%s\t%s\t%s" % (snp, deltasnp[snp], int(rankings[index]))
                #print "%s\t%s\t%s\t%s\t%s" % (snp, deltasnp[snp], rank, low_index[snp], high_index[snp])
            else:
                print "%s is not a valid SNP, please check that it is present in both High and Low bulks at and meets the specified threshold values" % (snp)

    current_out = "%s.%d.tsv" % (out, current_num)
    with open (current_out, 'w') as f:
        index = 0
        for el in sdeltasnp:
            if float(el[1]) >= float(args.report_threshold):
                #print "%s is >= %s" % (el[1], args.report_threshold)
                if args.ranks is True:
                    f.write("%s\t%s\t%s\n" % (el[0], el[1], int(rankings[index])))
                    index += 1
                else: f.write("%s\t%s\n" % (el[0], el[1]))
                #f.write("%s\t%s\t%s\t%s\n" % (el[0], el[1], low_index[el[0]], high_index[el[0]]))
            else: continue


    current_num -=1
    current_start += 1
    current_pos = 0
#print len(high_index)
#print len(high_index["dictindex"])
#Output calculated deltaSNP-index values for all valid positions
#for pos in high_index["dictindex"]:
#    print pos + "\t" + str(deltasnp[pos])
