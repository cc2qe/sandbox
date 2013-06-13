#!/usr/bin/env python

import sys
import getopt
import string
from optparse import OptionParser
from collections import Counter
import random

def processFile(inFile, sampleFile, maxHetRate, minMAF, maxMAF):
    minAF = float(minMAF)
    maxAF = float(maxMAF)
    AFD = 0.1
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
        for i in range(7,len(line)):
            line[i] = line[i].replace(line[3],'0')
            line[i] = line[i].replace(line[4],'1')
            line[i] = str(int(line[i][0]) + int(line[i][2]))

            # number of heterozygotes at locus.
            if line[i] == '1':
                numHet += 1

        gen = line[7:]
        numSamples = len(gen)
        maxHetCount = maxHetRate * numSamples

# now calculate the allele frequencies for the total set, and for each subpopulation
# check with original VCF file to make sure that we get the same numbers
# first get the total AF of snp1
        tCount = float(0)
        for i in range(0,len(gen)):
            tCount = tCount + float(gen[i])
        tAF = tCount / (2*numSamples)

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
        #Vt = t / numSamples;

# now get af for subpopulations
        eurCount = float(0)
        for val in eur:
            eurCount = eurCount + float(gen[val])
        eurAF = eurCount / (2*len(eur))

        afrCount = float(0)
        for val in afr:
            afrCount = afrCount + float(gen[val])
        afrAF = afrCount / (2*len(afr))

        asnCount = float(0)
        for val in asn:
            asnCount = asnCount + float(gen[val])
        asnAF = asnCount / (2*len(asn))

        amrCount = float(0)
        for val in amr:
            amrCount = amrCount + float(gen[val])
        amrAF = amrCount / (2*len(amr))

# now calculate variance among subpopulations
        #if tAF >= minAF and tAF <= maxAF and (abs(tAF-eurAF)/tAF) <= AFD and  (abs(tAF-afrAF)/tAF) <= AFD and (abs(tAF-asnAF)/tAF) <= AFD and (abs(tAF-amrAF)/tAF) <= AFD:
        #    print "tAF=" + str(tAF) + "\t" + "eurAF=" + str(eurAF) + "\t" + "afrAF=" + str(afrAF) + "\t" + "asnAF=" + str(asnAF) + "\t" + "amrAF=" + str(amrAF)
        # + "\t" + "totalVariance=" + str(Vt)

        if tAF >= minAF and tAF <= maxAF and (abs(tAF-eurAF)/tAF) <= AFD and  (abs(tAF-afrAF)/tAF) <= AFD and (abs(tAF-asnAF)/tAF) <= AFD and (abs(tAF-amrAF)/tAF) <= AFD and numHet <= maxHetCount:
            print line[0] + "\t" + line[1] + "\t" + \
            str(tAF)[0:4] + "\t" + str(eurAF)[0:4] + "\t" + str(afrAF)[0:4] + "\t" + str(asnAF)[0:4] + "\t" + str(amrAF)[0:4] + "\t",
            for item in line[7:]:
                print item + "\t",
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

    parser.add_option('-m', '--maxHetRate', dest='maxHetRate', default=0.8, help='maximum heterozygosity rate for a variant (default0.8)')

    parser.add_option('-p', '--minMAF', dest='minMAF', default=0.4, help='minimum minor allele frequency (default: 0.4')

    parser.add_option('-q', '--maxMAF', dest='maxMAF', default=0.6, help='maximum minor allele frequency (default: 0.6)')

    (opts, args) = parser.parse_args()

    if opts.inFile is None:
        parser.print_help()
        print
    else:
        processFile(opts.inFile, opts.sampleFile, opts.maxHetRate, opts.minMAF, opts.maxMAF)

if __name__ == "__main__":
    sys.exit(main()) 
    
