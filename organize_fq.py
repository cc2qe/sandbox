#!/usr/bin/env python

import argparse, sys
import os
import re

# --------------------------------------
# define functions

def organizeFastqs(fastqList):
    rootDir = '/mnt/thor_pool1/user_data/cc2qe/projects/human_diversity/1kg/batch1'   
    for line in fastqList:
        v = line.rstrip().split('\t')

        sampleId = v[10]
        runId = v[3]
        libraryId = v[15]
        pairInsert = v[18]
        center = v[6]
        platform = v[13]
        studyId = v[4]
        file = v[0]

        strand = int( re.sub('.*_(.)\.filt\.fastq\.gz', '\\1', file) )

        destDir = rootDir + '/' + sampleId + '/' + libraryId

        print sampleId, runId, file, libraryId, strand
        RGstring = '@RG\tID:%s\tLB:%s\tSM:%s\tPI:%s\tCN:%s\tPL:%s\tDS:%s' % (runId, libraryId, sampleId, pairInsert, center, platform, studyId)

        # check if directory exists
#        if not os.path.isdir(rootDir + '/' + sampleId):
            #print 'mkdir ' + rootDir + '/' + sampleId
 #           os.mkdir(rootDir + '/' + sampleId)
  #      if not os.path.isdir(destDir):
   #         print 'mkdir ' + destDir
    #        os.mkdir(destDir)
        
        # make readgroup sam flag
        

        #print 'moving %s to %s' % (file, destDir)
        #os.rename(file, destDir + '/' + file)
    return

# --------------------------------------
# argument parsing

def main():
    parser = argparse.ArgumentParser(description="Organize 1kg fastq files by sample")
    parser.add_argument('-i', '--input', required=True, type=argparse.FileType('r'), help='fastq filelist')
    
    # parse the arguments
    args = parser.parse_args()

    # store into global values
    fastqList = args.input
    
    organizeFastqs(fastqList)

    fastqList.close()

# initialize the script
if __name__ == '__main__':
    sys.exit(main())
