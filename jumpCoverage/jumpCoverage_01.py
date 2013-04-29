#!/usr/bin/env python

import sys, argparse

f = sys.stdin


class bamRec():
	def parse(self, line):
		(self.readName, self.flag, self.chrom, self.pos, self.mapQ, self.CIGAR, self.mateChrom, self.matePos, self.pairDist, self.seq, self.baseQual) = line.rstrip().split('\t')[:11]
		
		self.pos = int(self.pos)
		self.flag = int(self.flag)
		self.length = len(self.seq)
		self.pairDist = int(self.pairDist)

class region():
	def __init__(self, chrom, start, end, value):
		self.chrom = chrom
		self.start = start
		self.end = end
		self.value = value
		self.length = end - start

class bedChrom():
	def __init__(self, chrom, myRegion):
		self.chrom = chrom
		self.graph = [myRegion]
	
	def addRegion(self, newRegion):
#		print len(self.graph)
		i = 0
		while i < len(self.graph):
			mySlice = self.graph[i]
			if newRegion.start > mySlice.start and newRegion.start <= mySlice.end:
				leftInsert = self.graph[i]
				break
			i += 1
		
		j = i
		while j < len(self.graph):
			mySlice = self.graph[j]
			if newRegion.end > mySlice.start and newRegion.end <= mySlice.end:
				rightInsert = self.graph[j]
				break
			j += 1
		
#		print self.chrom, newRegion.chrom, newRegion.start, newRegion.end, i, j
		
		if leftInsert == rightInsert:
			addSlices = []
			addSlices.append( region(leftInsert.chrom, leftInsert.start, newRegion.start, leftInsert.value) )
			addSlices.append( region(leftInsert.chrom, newRegion.start, newRegion.end, leftInsert.value + newRegion.value) )
			addSlices.append( region(leftInsert.chrom, newRegion.end, leftInsert.end, leftInsert.value) )
			
			for m in addSlices:
				if m.length == 0: addSlices.remove(m)
			
			self.graph = self.graph[:i] + addSlices + self.graph[i+1:]
		
		else:
			addSlicesLeft = []
			addSlicesRight = []
			
			addSlicesLeft.append( region(leftInsert.chrom, leftInsert.start, newRegion.start, leftInsert.value) )
			addSlicesLeft.append( region(leftInsert.chrom, newRegion.start, leftInsert.end, leftInsert.value + newRegion.value) )
			for k in range(i+1,j):
				self.graph[k].value += newRegion.value
			addSlicesRight.append( region(rightInsert.chrom, rightInsert.start, newRegion.end, rightInsert.value + newRegion.value) )
			addSlicesRight.append( region(rightInsert.chrom, newRegion.end, rightInsert.end, rightInsert.value) )
			
			for m in addSlicesLeft:
				if m.length == 0: addSlicesLeft.remove(m)
			
			for m in addSlicesRight:
				if m.length == 0: addSlicesRight.remove(m)
			
			self.graph = self.graph[:i] + addSlicesLeft + self.graph[i+1:j] + addSlicesRight + self.graph[j+1:]
	
	# dump all regions to output before a certain position
	def dumpComplete(self, maxChrom, maxPos):

#		print maxPos, maxChrom
		i = 0
		while i < len(self.graph) and (self.graph[i].end < maxPos or self.graph[i].chrom != maxChrom):
			r = self.graph[i]
			print '\t'.join(map(str, [r.chrom, r.start, r.end, r.value]))
			i += 1
		self.graph = self.graph[i:]
			

parser = argparse.ArgumentParser(description='Generate a BED graph of jumping library coverage')
parser.add_argument('-g', '--genome', type=file, help='Tab delimited file of chromosome sizes in genome\n<chromName><TAB><chromSize>')

args = parser.parse_args()

genome = args.genome

#nitRegion = region()
#myChrom = bedChrom()

currentChrom = None
myChrom = None

for line in f:
	myRec = bamRec()
	myRec.parse(line)
	
	if myRec.chrom != currentChrom:
#		if myChrom != None:
#			for r in myChrom.graph:
#				print '\t'.join(map(str, [currentChrom, r.start, r.end, r.value]))

		if myChrom != None:
			myChrom.dumpComplete("dumpAll", 999999999)

		for c in genome:
			c = c.rstrip().split('\t')
			if c[0] == myRec.chrom:
				currentChrom = myRec.chrom
				initRegion = region(currentChrom, 0, int(c[1]), 0)
				myChrom = bedChrom(currentChrom, initRegion)
				
		genome.seek(0)

	if myRec.pairDist > 0:
		newR = region(myRec.chrom, myRec.pos, myRec.pos + myRec.pairDist, 1)
		myChrom.addRegion(newR)
	
		myChrom.dumpComplete(newR.chrom, newR.start)


#for r in myChrom.graph:
#	print '\t'.join(map(str, [currentChrom, r.start, r.end, r.value]))


genome.close()
f.close()





#
#
#
#
#
#	if myRec.chrom != indexChrom:
#		while len(valueArray) > 0:
#			i = 0
#			while valueArray[i] == valueArray[i + 1]:
#				i += 1
#			print myRec.chrom, indexStart, indexStart + i, valueArray[i-1]
#			indexStart += i
#			valueArray = valueArray[i+1:]
#		print "HIHI", len(valueArray), valueArray
#		indexStart = 0
#		indexChrom = myRec.chrom
#
#	offset = myRec.pos - indexStart
#
#	for j in range(offset, offset + myRec.pairDist):
#		while len(valueArray) < (j + 1):
#			valueArray.append(0)
#		
#		valueArray[j] += 1
##		print (valueArray)
#	
#	i = 0
#	while  len(valueArray) > i and valueArray[i] == valueArray[i + 1]:
#		i += 1
#		if i + 1 == len(valueArray): break
#	if len(valueArray) > i and indexStart < myRec.pos:
#		print myRec.chrom, indexStart, indexStart + i, valueArray[i-1]
#		indexStart += i
#		valueArray = valueArray[i+1:]
#
#
#
#
#
#
#
#
#

