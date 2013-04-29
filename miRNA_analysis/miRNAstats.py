#!/scr/talkowski/tools/bin/python/bin/python

import sys

def getMedian(numericValues):
  theValues = sorted(numericValues)
  if len(theValues) % 2 == 1:
    return theValues[(len(theValues)+1)/2-1]
  else:
    lower = theValues[len(theValues)/2-1]
    upper = theValues[len(theValues)/2]
    return (float(lower + upper)) / 2  

def getStdDev(values):
	n = len(values)
	m = sum(values)/n
	
	mySum = 0
	for x in values:
		mySum += (x-m)**2
	s = (mySum/n)**(0.5)
	return s
		

if len(sys.argv) != 2:
	print '\tusage: %s <coverageBed>' % sys.argv[0]
	exit()

f = open(sys.argv[1], 'r')

miRNA = None
covArray = []

print '#name\tminCov\tmaxCov\tmedianCov\tmeanCov\tstdDev'

for l in f:
	data = l.rstrip().split('\t')

	data[1] = int(data[1])
	data[2] = int(data[2])
	data[3] = int(data[3])
	
	if data[0] != miRNA:
		if miRNA:
			print '\t'.join(map(str, [miRNA, min(covArray), max(covArray), getMedian(covArray), sum(covArray)/len(covArray), getStdDev(covArray)]) )
		miRNA = data[0]
		covArray = []

	for i in range(data[1], data[2]):
		covArray.append(data[3])
f.close()
