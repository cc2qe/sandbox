#!/usr/bin/env python

import sys, subprocess, argparse
from string import *

#runs a command in the bash shell
def bash(cmd,cwd=None):
	# print(cmd)	# diagnostic
	retVal = subprocess.Popen(cmd, shell=True, \
		stdout=subprocess.PIPE, cwd=cwd).stdout.read().strip('\n').split('\n')
	if retVal==['']:
		return(0)
	else:
		return(retVal)

# gets the reverse complement of a DNA strand
def revcomp(dna):
    """ reverse complement of a DNA sequence """
    
    comp = dna.translate(maketrans("AGCTagct", "TCGAtcga"))
    lcomp = list(comp)
    lcomp.reverse()

    return ''.join(lcomp)

def GCcontent(dna):
	return ( dna.count('G') + dna.count('C') ) / float(len(dna))

def TM(dna):
	return float(bash('oligotm ' + dna)[0])
	#return 2 * (dna.count('A') + dna.count('T')) + 4 * (dna.count('G') + dna.count('C'))

#self anneal temperature (lower is better)
def selfAnneal(dna):
	primer3Dir = '/scr/talkowski/tools/bin/primer3-2.2.3/src/primer3_config/'
	return float(bash('ntthal -r -path ' + primer3Dir + ' -s1 ' + dna + ' -s2 ' + dna)[0])

#1	6	HWI-ST177_0158:7:7:2454:166639#0	OAR11	HD_transgene	47789106	11307	37	0	25	25	113	177	AATCTTATATTTGTAAAGAGTAAAG	ACTACGGTCTGCTCTCCTGCTTCCG	HHFHHHHHHHHHHHHHHHHHHHHHH	HHHHHHHFHHGHFHHHHHHHHHHHH

class readPair():
	def parse(self, line):
		(self.clusterID, self.clusterCount, self.readName, self.chrA, self.chrB, self.posA, self.posB, self.mapQA, self.mapQB, self.lenA, self.lenB, self.flagA, self.flagB, self.seqA, self.seqB, self.qualA, self.qualB) = line.rstrip().split('\t')
		
		self.posA = int(self.posA)
		self.posB = int(self.posB)
		self.flagA = int(self.flagA)
		self.flagB = int(self.flagB)
		self.lenA = int(self.lenA)
		self.lenB = int(self.lenB)
	
	def getStrand(self, AorB):
		if AorB == 'A':
			myFlag = self.flagA
		elif AorB == 'B':
			myFlag = self.flagB
		
		if myFlag % 32 < 16:
			strand = '+'
		elif myFlag % 32 >= 16:
			strand = '-'
		
		if orientation == 'inward':
			return strand
		elif orientation == 'outward':
			if strand == '+': strand = '-'
			elif strand == '-': strand = '+'
			
			return strand
	
	def export(self):
		print '\t'.join(map(str, (self.readName, self.chrA, self.chrB, self.posA, self.posB, self.mapQA, self.mapQB, self.lenA, self.lenB, self.flagA, self.flagB, self.seqA, self.seqB, self.qualA, self.qualB)))

parser = argparse.ArgumentParser(description='convert BamStat clustering file to primers. (input is stdin)')
parser.add_argument('-o', '--orientation', required=True, choices=['inward', 'outward'], help='read orientation of alignment')
parser.add_argument('-min', default=20, help='min primer length (base pairs)')
parser.add_argument('-max', default=30, help='max primer length (base pairs)')
parser.add_argument('-t', '--targetTm', default=60, help='minimum TM (will add up to 5 bases to reach this TM)')

args = parser.parse_args()
orientation = args.orientation
primerMin = int(args.min)
primerMax = int(args.max)
targetTm = int(args.targetTm)

f = sys.stdin


print '\t'.join(['#chrA', 'chrB', 'posA', 'posB', 'strA', 'strB', 'mapQA', 'mapQB', 'primerA', 'primerB', 'len A', 'len B', 'Tm A', 'Tm B', 'GC A', 'GC B', 'self Tm A', 'self Tm B'])

for l in f:
	if l != '\n' and list(l)[0] != '\t':
		current = readPair()
		current.parse(l)
		#current.export()

		if current.getStrand('A') == '-':
			primerA = revcomp(current.seqA)
		else:
			primerA = current.seqA
		if current.getStrand('B') == '-':
			primerB = revcomp(current.seqB)
		else:
			primerB = current.seqB
		
		primerLenA = primerMin
		primerLenB = primerMin
		
		if orientation == 'inward':
			while TM(primerA[:primerLenA]) < targetTm and primerLenA < primerMax:
				primerLenA += 1
			while TM(primerB[:primerLenB]) < targetTm and primerLenB < primerMax:
				primerLenB += 1
			
			if abs(TM(primerA[:primerLenA-1]) - targetTm) < abs(TM(primerA[:primerLenA]) - targetTm) and primerLenA > primerMin:
				primerA = primerA[:primerLenA-1]
			else:
				primerA = primerA[:primerLenA]

			if abs(TM(primerB[:primerLenB-1]) - targetTm) < abs(TM(primerB[:primerLenB]) - targetTm) and primerLenB > primerMin:
				primerB = primerB[:primerLenB-1]
			else:
				primerB = primerB[:primerLenB]
		
		elif orientation == 'outward':
			while TM(primerA[-primerLenA:]) < targetTm and primerLenA < primerMax:
				primerLenA += 1
			
			if abs(TM(primerA[-(primerLenA-1):]) - targetTm) < abs(TM(primerA[-(primerLenA):]) - targetTm and primerLenA > primerMin):
				primerA = primerA[-(primerLenA-1):]
			else:
				primerA = primerA[-(primerLenA):]
			
			while TM(primerB[-primerLenB:]) < targetTm and primerLenB < primerMax:
				primerLenB += 1
			
			if abs(TM(primerB[-(primerLenB-1):]) - targetTm) < abs(TM(primerB[-primerLenB:]) - targetTm and primerLenB > primerMin):
				primerB = primerB[-(primerLenB-1):]
			else:
				primerB = primerB[-primerLenB:]
		
#		print 'primerA', primerA
#		print 'primerB', primerB
		
		print '\t'.join(map(str, [current.chrA, current.chrB, current.posA, current.posB, current.getStrand('A'), current.getStrand('B'), current.mapQA, current.mapQB, primerA, primerB, len(primerA), len(primerB), TM(primerA), TM(primerB), GCcontent(primerA), GCcontent(primerB), selfAnneal(primerA), selfAnneal(primerB) ]))
	
	else: print
		
