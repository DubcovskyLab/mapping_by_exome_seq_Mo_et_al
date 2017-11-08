#!/usr/bin/env python
## By: Hans Vasquez-Gross   
## Last Updated: 11/07/2017
from __future__ import division
from optparse import OptionParser
import pprint
import re

parser = OptionParser(usage="usage: %prog [options]\n\nPrints to STDOUT a TSV file of counted NTs and called SNPs", version="%prog 1.6")
parser.add_option("-p", "--pileup", help="filtered mpileup file based off only positions from the M2")
parser.add_option("-m", "--mutation", help="MAPS2 file filtered for the particular line")
parser.add_option("-d", "--depth", type="int", dest="lowdepth", help="Minimum depth to be considered in the analysis; Default 5 reads", default=5)
parser.add_option("-o", "--hom", type="float", dest="homprcnt", help="Minimum percentage to be considered homozygous; Default 0.9", default=0.9)
parser.add_option("-u", "--upper_het", type="float", dest="upperhet", help="Upper limit to be considered heterozygous; Default 0.7", default=0.7)
parser.add_option("-l", "--lower_het", type="float", dest="lowerhet", help="Lower limit to be considered heterozygous; Default 0.3", default=0.3)
(options, args) = parser.parse_args()


##Function Definitions
def initDictFromHeader(header):

    dictionary = dict()
    sampleOrderList = list()
    headerList = header.split("\t")
    contig = headerList.pop(0)
    pos = headerList.pop(0)
    ref = headerList.pop(0)
    #print "Num entries:", len(headerList), "Num samples:", len(headerList) / 3
    num_samples = int(len(headerList) / 3)
    for i in range(1, num_samples + 1):
        cov = headerList.pop(0)
        call = headerList.pop(0)
        qual = headerList.pop(0)
        sample_name = cov.replace("Cov-", "")
        #print i, sample_name, cov, call, qual
        dictionary.setdefault(sample_name, dict())
        sampleOrderList.append(sample_name)

    return (dictionary, sampleOrderList)


def uniqify(seq, idfun=None): 
   # order preserving
   if idfun is None:
       def idfun(x): return x
   seen = {}
   result = []
   for item in seq:
       marker = idfun(item)
       if marker in seen: continue
       seen[marker] = 1
       result.append(item)
   return result

##Main Program
if __name__ == "__main__":
    if(options.pileup):
        mpileup_fh = open(options.pileup, "r")
    else:
        raise Exception('Mpileup file open failed!')

    if(options.mutation):
        mutation_fh = open(options.mutation, "r")
    else:
        raise Exception('Mutation file open failed!')


    ##BUILD the data of mutation
    #header = mutation_fh.readline()
    mapsMutationDict = dict()
    for line in mutation_fh:
        lineList = line.split("\t")
	##FOR downloaded TSV
        #contig = lineList[0]
        #pos = lineList[2]
        #ref = lineList[3]
        #wt = lineList[4]
        #mut = lineList[5]
        #hethom = lineList[6]

	##FOR MAPS grepped TSV
        contig = lineList[0]
        pos = lineList[1]
        ref = lineList[2]
        wt = lineList[4]
        mut = lineList[5]
        hethom = lineList[7]
        mapsMutationDict.setdefault(contig, dict())
        mapsMutationDict[contig][pos] = (ref, wt, mut, hethom)


    
    ##Build the data structure to organize the data
    header = mpileup_fh.readline()
    (dataDict, sampleOrderList) = initDictFromHeader(header)

    #print "SampleOrderList"
    #pprint.pprint(sampleOrderList)

    contigOrderList = list()
    contigDataDict = dict()
    for line in mpileup_fh:
        line = line.strip()
        lineList = line.split("\t")

        contig = lineList.pop(0)
        pos = lineList.pop(0)
        ref = lineList.pop(0)
        num_samples = int(len(lineList) / 3)
        #print "Num entries:", len(lineList), "Num samples:", num_samples
        contigOrderList.append(contig)
        contigDataDict.setdefault(contig, set())
        contigDataDict[contig].add((pos, ref))

        for i in range(0, num_samples ):
            sample_name = sampleOrderList[i]
            cov = lineList.pop(0)
            call = lineList.pop(0).upper()
            qual = lineList.pop(0)

            dataDict[sample_name].setdefault(contig, dict())
            dataDict[sample_name][contig][pos] = (cov, call, qual)


    ##Uniq the list but maintain order for looping
    contigOrderList = uniqify(contigOrderList)

    ##print header
    out_line = "contig\tpos\tref\twt\tmut\thethom"
    for sample_name in sampleOrderList:
        out_line += "\t%s-A\t%s-C\t%s-G\t%s-T\t%s-callsnp" % (sample_name, sample_name, sample_name, sample_name, sample_name)
    print out_line


    for contig in contigOrderList:
        for values in contigDataDict[contig]:
            pos = values[0]
            ref = values[1]
            wt = mapsMutationDict[contig][pos][1]
            mut = mapsMutationDict[contig][pos][2]
            hethom = mapsMutationDict[contig][pos][3]
            out_line = "%s\t%s\t%s\t%s\t%s\t%s" % (contig, pos, ref, wt,mut,hethom)
            for sample_name in sampleOrderList:
                ##init values
                count_A = 0
                count_T = 0
                count_G = 0
                count_C = 0
                if(contig in dataDict[sample_name]):

                    ##Get cov and call, and remove deletion data
                    if(pos in dataDict[sample_name][contig]):
                        cov = dataDict[sample_name][contig][pos][0]
                        call = dataDict[sample_name][contig][pos][1]
                    else:
                        cov = 0
                        call = "NA"
                    call = re.sub(r'[-\+]\d+[ACGT]+', "", call) ##Remove the indels
                    call = re.sub(r'\^\w+', "", call) ## Remove the anchored read starts with quals

                    ##COUNT data
                    count_A = call.count('A')
                    count_T = call.count('T')
                    count_G = call.count('G')
                    count_C = call.count('C')
                    count_N = call.count('N')
                    count_comma = call.count(',')
                    count_period = call.count('.')
                    count_dollar = call.count('$')
                else:
                    count_A = 0
                    count_T = 0
                    count_G = 0
                    count_C = 0
                    count_N = 0
                    count_comma = 0
                    count_period = 0
                    count_dollar = 0

                count_ref = count_comma + count_period 
                if(ref == 'A'):
                    count_A += count_ref
                elif(ref == 'C'):
                    count_C += count_ref
                elif(ref == 'G'):
                    count_G += count_ref
                elif(ref == 'T'):
                    count_T += count_ref
                total_count = count_A + count_T + count_G + count_C + count_N
       
                #Compute percentage
                if total_count >= options.lowdepth:
                    percnt_A = count_A / total_count
                    percnt_C = count_C / total_count
                    percnt_G = count_G / total_count
                    percnt_T = count_T / total_count
        
                    if(percnt_A >= options.homprcnt):
                        call_snp = "AA"
                    elif(percnt_C >= options.homprcnt):
                        call_snp = "CC"
                    elif(percnt_G >= options.homprcnt):
                        call_snp = "GG"
                    elif(percnt_T >= options.homprcnt):
                        call_snp = "TT"
                    elif(percnt_A >= options.lowerhet and percnt_A <= options.upperhet):
                    #if between 30% to 70%, call het
                        call_snp = "HET"
                    elif(percnt_C >= options.lowerhet and percnt_C <= options.upperhet):
                        call_snp = "HET"
                    elif(percnt_G >= options.lowerhet and percnt_G <= options.upperhet):
                        call_snp = "HET"
                    elif(percnt_T >= options.lowerhet and percnt_T <= options.upperhet):
                        call_snp = "HET"
                    else:
                        call_snp = "AMBIG"
                else:
                    call_snp = "LOW_DEPTH"

                out_line += "\t%s\t%s\t%s\t%s\t%s" % (count_A, count_C, count_G, count_T, call_snp)

            print out_line

