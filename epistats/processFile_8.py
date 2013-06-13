#!/usr/bin/env python

import sys
import getopt
import string
from optparse import OptionParser
from collections import Counter
import random

def processFile(inFile, sampleFile, maxHetRate, minMAF, maxMAF, AFD):
    minAF = float(minMAF)
    maxAF = float(maxMAF)
    AFD = float(AFD)
    maxHetRate = float(maxHetRate)
    
# load the sample data
    file = open(sampleFile,'r')
    sample = []
    for line in file:
        sample.append(line.strip().split('\t'))

# get vectors which are the indices of the samples for each of ASN, AMR, AFR & EUR superpopulations;
    amr = set()
    afr = set() 
    asn = set()
    eur = set()

    for i in range(0,len(sample)):
        if sample[i][2] == 'AMR':
            amr.add(i);
        elif sample[i][2] == 'AFR':
            afr.add(i);
        elif sample[i][2] == 'ASN':
            asn.add(i);
        elif sample[i][2] == 'EUR':
            eur.add(i)

# load the snp data
    if inFile == "stdin":
        data = sys.stdin
    else:
        data = open(inFile, 'r')
    snp = []
    for line in data:
        line = line.strip().split('\t')
        numHet = 0

        # the number of samples that are not ./. at the locus
        numInformative = 0;
        for i in range(8,len(line)):
            if line[i] == "./.":
                line[i] = -1
            else:
                line[i] = line[i].replace(line[3],'0')
                line[i] = line[i].replace(line[4],'1')
                line[i] = str(int(line[i][0]) + int(line[i][2]))
                numInformative += 1

            # number of heterozygotes at locus.
            if line[i] == '1':
                numHet += 1

        gen = line[8:]
        maxHetCount = maxHetRate * numInformative

# now calculate the allele frequencies for the total set, and for each subpopulation
# check with original VCF file to make sure that we get the same numbers
# first get the total AF of snp1
        tCount = float(0)
        for i in range(0,len(gen)):
            if gen[i] != -1:
                tCount = tCount + float(gen[i])
        tAF = tCount / (2*numInformative)

# now calculate the total variance:
# how to calculate variance: 
# 1) first calculate the mean M
# 2) for each value Vi, Di = (Vi-M)^2
# 3) variance is the mean of the D values
        # first the mean:
        # tMean = tCount / numSamples
        # now the difference:
        #t = 0
        #for i in range(0,len(gen)):
        #    t = t + ((float(gen[i]) - tMean)**2)
        #Vt = t / numInformative;

# now get af for subpopulations
        eurCount = float(0)
        # number of informative eur at that locus
        eurInf = 0
        for val in eur:
            if gen[val] != -1:
                eurCount = eurCount + float(gen[val])
                eurInf += 1
        if eurInf != 0:
            eurAF = eurCount / (2*eurInf)
        else: eurAF = -1

        afrCount = float(0)
        afrInf = 0
        for val in afr:
            if gen[val] != -1:
                afrCount = afrCount + float(gen[val])
                afrInf += 1
        if afrInf != 0:
            afrAF = afrCount / (2*afrInf)
        else: afrAF = -1

        asnCount = float(0)
        asnInf = 0
        for val in asn:
            if gen[val] != -1:
                asnCount = asnCount + float(gen[val])
                asnInf += 1
        if asnInf != 0:
            asnAF = asnCount / (2*asnInf)
        else: asnAF = -1

        amrCount = float(0)
        amrInf = 0
        for val in amr:
            if gen[val] != -1:
                amrCount = amrCount + float(gen[val])
                amrInf += 1
        if amrInf != 0:
            amrAF = amrCount / (2*amrInf)
        else: amrAF = -1

# now calculate variance among subpopulations
        #if tAF >= minAF and tAF <= maxAF and (abs(tAF-eurAF)/tAF) <= AFD and  (abs(tAF-afrAF)/tAF) <= AFD and (abs(tAF-asnAF)/tAF) <= AFD and (abs(tAF-amrAF)/tAF) <= AFD:
        #    print "tAF=" + str(tAF) + "\t" + "eurAF=" + str(eurAF) + "\t" + "afrAF=" + str(afrAF) + "\t" + "asnAF=" + str(asnAF) + "\t" + "amrAF=" + str(amrAF)
        # + "\t" + "totalVariance=" + str(Vt)

        if tAF >= minAF and tAF <= maxAF and (abs(tAF-eurAF)/tAF) <= AFD and  (abs(tAF-afrAF)/tAF) <= AFD and (abs(tAF-asnAF)/tAF) <= AFD and (abs(tAF-amrAF)/tAF) <= AFD and numHet <= maxHetCount:
            print line[0] + "\t" + line[2] + "\t" + line[5] + '\t' + \
            str(tAF)[0:4] + "\t" + str(amrAF)[0:4] + "\t" + str(afrAF)[0:4] + "\t" + str(asnAF)[0:4] + "\t" + str(eurAF)[0:4] + "\t",
            for item in line[8:]:
                print str(item) + "\t",
            print

          
#---------------------------------------------------------------------------
# argument parsing

class Usage(Exception):
    def __init__(self, msg):
        self.msg = msg		

def main():	
    usage = """%prog -i <file> 

getFrequency
Author: Ira Hall	
Description: 
    """
    parser = OptionParser(usage)

    parser.add_option("-i", "--file", dest="inFile", 
        help="A tab delimited VCF-like file with genotypes starting at column 8 or standard input (-i stdin).",
        metavar="FILE")
        
    parser.add_option("-s", "--sampleFile", dest="sampleFile", 
        help="A tab delimited file with sample, population, and superpopulation",
        metavar="FILE")        

    parser.add_option('-m', '--maxHetRate', dest='maxHetRate', default=0.8, help='maximum heterozygosity rate for a variant (default: 0.8)')

    parser.add_option('-p', '--minMAF', dest='minMAF', default=0.4, help='minimum minor allele frequency (default: 0.4')

    parser.add_option('-q', '--maxMAF', dest='maxMAF', default=0.6, help='maximum minor allele frequency (default: 0.6)')

    parser.add_option('-a', '--alleleFreqDiff', dest='alleleFreqDiff', default=0.1, help='maximum difference in allele frequency between subpopulations abs(totalAlleleFreq - subpopFreq) / totalAllelFreq < alleleFreqDiff  (default 0.1)')

    (opts, args) = parser.parse_args()

    if opts.inFile is None:
        parser.print_help()
        print
    else:
        processFile(opts.inFile, opts.sampleFile, opts.maxHetRate, opts.minMAF, opts.maxMAF, opts.alleleFreqDiff)

if __name__ == "__main__":
    sys.exit(main()) 
    
