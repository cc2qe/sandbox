#!/scr/talkowski/tools/bin/python/bin/python

import sys

if len(sys.argv) != 3:
	print '\tusage: %s <genomeChromSizes> <bed file>' % sys.argv[0]
	exit()

f1 = open(sys.argv[1], 'r')
f2 = open(sys.argv[2], 'r')

print f2.readline().rstrip() + '\tlength'	# burn header of bed file

l1 = f1.readline().rstrip()
l2 = f2.readline().rstrip()

while l1 != '':
	data1 = l1.split('\t')
	data2 = l2.split('\t')
	
	if data1[0] == data2[0]:
		print '%s\t%s' % (l2, data1[1])
		l1 = f1.readline().rstrip()
		l2 = f2.readline().rstrip()
	else:
		print '%s\t0\t0\t0\t0\t0\t%s' % (data1[0], data1[1])
		l1 = f1.readline().rstrip()

f1.close()
f2.close()
