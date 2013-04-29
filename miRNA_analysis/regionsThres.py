#!/shr/home/cc472/bin/python/bin/python

import sys

def getMedian(numericValues):
	theValues = sorted(numericValues)

	if len(theValues) % 2 == 1:
		return theValues[(len(theValues)+1)/2-1]
	else:
		lower = theValues[len(theValues)/2-1]
		upper = theValues[len(theValues)/2]
		return (float(lower + upper)) / 2

if len(sys.argv) != 3:
	print "\n\tusage: %s [threshold] [bedFile]\n" % sys.argv[0]
	exit(1)


threshold = int(sys.argv[1])
filename = sys.argv[2]

f = open(filename, 'r')

extend = False

for l in f:
	lineData = l.rstrip().split('\t')
	#print lineData
	#print extend
	lineData[3] = int(lineData[3])


	if lineData[3] >= threshold:
		if extend and lineData[1] == prevLineData[2] and lineData[0] == prevLineData[0]:
			regEnd = lineData[2]
			regVals.append(lineData[3])
		
		elif extend:
			print regChr, regStart, regEnd, sum(regVals)/len(regVals), getMedian(regVals)
			
                        regChr = lineData[0]
                        regStart = lineData[1]
                        regEnd = lineData[2]
                        regVals = [lineData[3]]

                else:
                        regChr = lineData[0]
                        regStart = lineData[1]
                        regEnd = lineData[2]
                        regVals = [lineData[3]]
			
			extend = True
	
	elif extend:
		print regChr, regStart, regEnd, sum(regVals)/len(regVals), getMedian(regVals)

		extend = False
		
	prevLineData = lineData


