#!/usr/bin/env python

import sys
import getopt
import string
from optparse import OptionParser
from collections import Counter
import random

def getFrequency(inFile, numIter, minDist, tag, yates, sampleFile):
    # load the sample data
    sampleFile = open(sampleFile, 'r')
    sample = []
    for line in sampleFile:
        sample.append(line.rstrip().split('\t'))

    # make list of the different types of subpopulations
    subpop_list = ['amr', 'afr', 'asn', 'eur']

    # get vectors which are the indices of the samples for each of ASN, AMR, AFR & EUR superpopulations;
    subpop_set = dict()
    for subpop in subpop_list:
        subpop_set[subpop] = set()

    for i in range(0,len(sample)):
        if sample[i][2] == 'AMR':
            subpop_set['amr'].add(i);
        elif sample[i][2] == 'AFR':
            subpop_set['afr'].add(i);
        elif sample[i][2] == 'ASN':
            subpop_set['asn'].add(i);
        elif sample[i][2] == 'EUR':
            subpop_set['eur'].add(i)

    # make the 27 genotype combinations for each population.
    temp = ["0","1","2"]
    combos= []
    for a in temp:
        for b in temp:
            for c in temp:
                combos.append([a,b,c])

    # load the genotype data 
    if inFile == "stdin":
        data = sys.stdin
    else:
        data = open(inFile, 'r')
    snp = []
    for line in data:
        snp.append(line.strip().split('\t'))
    numSamples = float(len(snp[1][8:]))
    numSnps = float(len(snp))
    count = 0
    while count < numIter:
        index = random.sample(range(0,int(numSnps)),3)
        snp1 = snp[index[0]]
        snp2 = snp[index[1]]
        snp3 = snp[index[2]]
        if (snp1[0]!=snp2[0] or abs(float(snp1[1])-float(snp2[1]))>=minDist) and (snp1[0]!=snp3[0] or abs(float(snp1[1])-float(snp3[1]))>=minDist) and (snp2[0]!=snp3[0] or abs(float(snp2[1])-float(snp3[1]))>=minDist):
            count = count + 1
            id = str(tag) + "_" + str(count)
            gen1 = snp1[8:]
            gen2 = snp2[8:]
            gen3 = snp3[8:]
            # make the counter for the entire population of samples
            cnt1=Counter()
            cnt2=Counter()
            cnt3=Counter()

            # make the counter for each of the subpopulations
            subpop_cnt = dict()
            for subpop in subpop_list:
                subpop_cnt[subpop] = [Counter(), Counter(), Counter()]

            # do the counting for the entire population, as well as each of the subpopulations.
            for t in range(len(gen1)):
                word = gen1[t]
                # iterate the total pop count
                cnt1[word] += 1
                # iterate the subpop counts
                if t in subpop_set['amr']:
                    subpop_cnt['amr'][0][word] += 1
                elif t in subpop_set['afr']:
                    subpop_cnt['afr'][0][word] += 1
                elif t in subpop_set['asn']:
                    subpop_cnt['asn'][0][word] += 1
                elif t in subpop_set['eur']:
                    subpop_cnt['eur'][0][word] += 1

            for t in range(len(gen2)):
                word = gen2[t]
                # iterate the total pop count
                cnt2[word] += 1
                # iterate the subpop counts
                if t in subpop_set['amr']:
                    subpop_cnt['amr'][1][word] += 1
                elif t in subpop_set['afr']:
                    subpop_cnt['afr'][1][word] += 1
                elif t in subpop_set['asn']:
                    subpop_cnt['asn'][1][word] += 1
                elif t in subpop_set['eur']:
                    subpop_cnt['eur'][1][word] += 1

            for t in range(len(gen3)):
                word = gen3[t]
                # iterate the total pop count
                cnt3[word] += 1
                # iterate the subpop counts
                if t in subpop_set['amr']:
                    subpop_cnt['amr'][2][word] += 1
                elif t in subpop_set['afr']:
                    subpop_cnt['afr'][2][word] += 1
                elif t in subpop_set['asn']:
                    subpop_cnt['asn'][2][word] += 1
                elif t in subpop_set['eur']:
                    subpop_cnt['eur'][2][word] += 1

			
            # now calculate the frequency of each genotype. af is allele frequency, ct is count
            af1 = []
            af2 = []
            af3 = []
            ct1 = []
            ct2 = []
            ct3 = []
            for i in range(0,3):
                 af1.append(float(cnt1[str(i)]) / numSamples)
                 af2.append(float(cnt2[str(i)]) / numSamples)
                 af3.append(float(cnt3[str(i)]) / numSamples)
                 ct1.append(float(cnt1[str(i)]))
                 ct2.append(float(cnt2[str(i)]))
                 ct3.append(float(cnt3[str(i)]))
            
            # now calculate the frequency of each genotype in each of the subpopulations.
            subpop_af = dict()
            subpop_ct = dict()
            for subpop in subpop_list:
                subpop_af[subpop] = [[], [], []]
                subpop_ct[subpop] = [[], [], []]

                for i in range(3):
                    # clean this up with a for loop if it works!
                    myCnt1 = subpop_cnt[subpop][0]
                    myCnt2 = subpop_cnt[subpop][1]
                    myCnt3 = subpop_cnt[subpop][2]

                    subpop_af[subpop][0].append(float(myCnt1[str(i)]) / len(subpop_set[subpop]))
                    subpop_af[subpop][1].append(float(myCnt2[str(i)]) / len(subpop_set[subpop]))
                    subpop_af[subpop][2].append(float(myCnt3[str(i)]) / len(subpop_set[subpop]))

                    subpop_ct[subpop][0].append(float(myCnt1[str(i)]))
                    subpop_ct[subpop][1].append(float(myCnt2[str(i)]))
                    subpop_ct[subpop][2].append(float(myCnt3[str(i)]))
                    
            # now calculate the expected frequency of each combination
            # make counter for the total population
            comboCount=Counter()
            # Make a counter for each of the subpopulations
            subpop_comboCount = dict()
            for subpop in subpop_list:
                subpop_comboCount[subpop] = Counter()
            for i in range(0,int(numSamples)):
                # gen1, gen2, and gen3 are strings. They are concatenated to something
                # like '102' for 'het,hom_ref, hom_alt'
                comboCount[gen1[i] + gen2[i] + gen3[i]] += 1
                if i in subpop_set['amr']:
                    subpop_comboCount['amr'][gen1[i] + gen2[i] + gen3[i]] += 1
                elif i in subpop_set['afr']:
                    subpop_comboCount['afr'][gen1[i] + gen2[i] + gen3[i]] += 1
                elif i in subpop_set['asn']:
                    subpop_comboCount['asn'][gen1[i] + gen2[i] + gen3[i]] += 1
                elif i in subpop_set['eur']:
                    subpop_comboCount['eur'][gen1[i] + gen2[i] + gen3[i]] += 1
               	

            # dSum is the sum of the probability differences between obs and expected
            dSum = 0
            # pSum is the sum of the chi squared values at each genotype
            pSum = 0

            # dSum and pSum for each of the subpops
            subpop_dSum = dict()
            subpop_pSum = dict()
            for subpop in subpop_list:
                subpop_dSum[subpop] = 0
                subpop_pSum[subpop] = 0

            for c in combos:
                # get the observed and expected for the total population.
                # c2 is the full genotype combo (e.g. 201 for hom_alt, hom_ref, het
                c2 = c[0] + c[1] + c[2]
                obs = float(comboCount[c2])/numSamples
                exp = float((af1[int(c[0])] * af2[int(c[1])] * af3[int(c[2])]))
                obsMinusExp = obs - exp
                dSum = dSum + abs(obsMinusExp)

                # calculate pearson's chi-square test stat
                # before i was doing frequency, now put the numbers
                # chi = ((obs-exp)**2) / exp

                # if yates flag is set, calculate yate's chi-square
                # correction for continuity
                # chi yates = ((abs(obs-exp)-0.5)**2)/exp
                expNum = exp * numSamples
                obsNum = float(comboCount[c2])
                if expNum == 0:
                    # print "DIVIDE BY ZERO ERROR HERE" 
                    chi = 0
                else:
                    # print ct1
                    # print ct2
                    # print ct3
                    # print str(obs) + "\t" + str(exp) + "\t" + str(expNum)
                    if yates:
                        chi = (((abs(obsNum-expNum))-0.5)**2) / expNum
                    else:
                        chi = ((obsNum-expNum)**2) / expNum
                pSum = pSum + chi

                # calculate chi-squared for subpopulations
                subpop_obs = dict()
                subpop_exp = dict()
                subpop_obsNum = dict()
                subpop_expNum = dict()
                subpop_chi = dict()
                for subpop in subpop_list:
                    subpop_obs[subpop] = float(subpop_comboCount[subpop][c2]) / len(subpop_set[subpop])
                    subpop_exp[subpop] = float(subpop_af[subpop][0][int(c[0])] * subpop_af[subpop][1][int(c[1])] * subpop_af[subpop][2][int(c[2])])
                    subpop_obsNum[subpop] = float(subpop_comboCount[subpop][c2])
                    subpop_expNum[subpop] = subpop_exp[subpop] * len(subpop_set[subpop])

                    # calc chi-squared over the subpops
                    if subpop_expNum[subpop] == 0:
                        # avoid divide by zero error
                        subpop_chi[subpop] = 0
                    elif yates:
                        subpop_chi[subpop] = (((abs(subpop_obsNum[subpop]-subpop_expNum[subpop]))-0.5)**2) / subpop_expNum[subpop]
                    else:
                        subpop_chi[subpop] = ((subpop_obsNum[subpop] - subpop_expNum[subpop])**2) / subpop_expNum[subpop]
                    
                    subpop_pSum[subpop] += subpop_chi[subpop]
                
                print '\t'.join(map(str, ['COMBO', id, snp1[0], snp1[1], snp1[2], round(af1[0],2), round(af1[1],2), round(af1[2],2), snp2[0], snp2[1], snp2[2], round(af2[0],2), round(af2[1],2), round(af2[2],2), snp3[0], snp3[1], snp3[2], round(af3[0],2), round(af3[1],2), round(af3[2],2), c2, round(obs,2), round(exp,2), round(obsMinusExp,2), round(subpop_chi['amr'],2), round(subpop_chi['afr'],2), round(subpop_chi['asn'],2), round(subpop_chi['eur'],2), round(chi,2)]))


            print '\t'.join(map(str, ['SUM', id, snp1[0], snp1[1], snp1[2], round(af1[0],2), round(af1[1],2), round(af1[2],2), snp2[0], snp2[1], snp2[2], round(af2[0],2), round(af2[1],2), round(af2[2],2), snp3[0], snp3[1], snp3[2], round(af3[0],2), round(af3[1],2), round(af3[2],2), 'NA', 'NA', 'NA', round(dSum,2), round(subpop_pSum['amr'],2), round(subpop_pSum['afr'],2), round(subpop_pSum['asn'],2), round(subpop_pSum['eur'],2), round(pSum,2)]))
            
        else: continue

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
        help="A tab delimited file with chrom, pos as cols 1-2 and genotypes starting at column 8 (or standard input, -i stdin).",
        metavar="FILE")

    parser.add_option("-n", "--numIter", dest="numIter", default=1000, type = "int",
        help="The number of iterations for sampling; Default = 1000.",
        metavar="INT")

    parser.add_option("-d", "--minDist", dest="minDist", default=1000000, type = "float",
        help="The minimum distance between loci; Default = 1000000.",
        metavar="FLOAT")

    parser.add_option("-t", "--tag", dest="tag", default="1", type = "str",
        help="The number of iterations for sampling; Default = 1.",
        metavar="STR")

    parser.add_option("-y", "--yates", dest="yates", action="store_true", default=False, help='Yates correction on chi-squared probability test')

    parser.add_option('-s', '--sampleFile', dest='sampleFile', help='a tab delimited file with sample, population, superpopulation', metavar='FILE')
        
    (opts, args) = parser.parse_args()

    if opts.inFile is None:
        parser.print_help()
        print
    else:
        getFrequency(opts.inFile, opts.numIter, opts.minDist, opts.tag, opts.yates, opts.sampleFile)

if __name__ == "__main__":
    sys.exit(main()) 
