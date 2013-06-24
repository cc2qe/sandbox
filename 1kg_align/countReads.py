#!/usr/bin/env python

import argparse, sys
import os
import re

# --------------------------------------
# define functions

def countReads(fastqList, destRoot):
    # destRoot = '/mnt/thor_pool1/user_data/cc2qe/projects/human_diversity/1kg/batch1'   
    for line in fastqList:
        v = line.rstrip().split('\t')

        sampleId = v[10]
        runId = v[3]
        libraryId = v[15]
        pairInsert = v[18]
        center = v[6]
        platform = v[13]
        studyId = v[4]
        readCount = v[24]
        file = v[0]

        strand = int( re.sub('.*_(.)\.filt\.fastq\.gz', '\\1', file) )

        destDir = destRoot + '/' + sampleId

        print sampleId, runId, file, libraryId, strand
        RGstring = '@RG\\tID:%s\\tLB:%s\\tSM:%s\\tPI:%s\\tCN:%s\\tPL:%s\\tDS:%s' % (runId, libraryId, sampleId, pairInsert, center, platform, studyId)

        # check if directory exists
        if not os.path.isdir(destDir):
            print 'mkdir ' + destDir
            os.mkdir(destRoot + '/' + sampleId)

        # print read count
        readCountFile = open(destDir + '/' + runId + '_count.txt', 'w')
        readCountFile.write(readCount)
        readCountFile.close()

        # make readgroup sam flag
        # readGroupFile = open(destDir + '/' + runId + '_readgroup.txt', 'w')
        # readGroupFile.write(RGstring)
        # readGroupFile.close()
        
        # verbose
        # print 'moving %s to %s' % (file, destDir)
        
        # the actual move command
        # os.rename(file, destDir + '/' + file)
    return

# --------------------------------------
# argument parsing

def main():
    parser = argparse.ArgumentParser(description="count reads in the fastq files")
    parser.add_argument('-i', '--input', required=True, type=argparse.FileType('r'), help='fastq filelist (fqbind)')
    parser.add_argument('-d', '--dest', required=True, type=str, help='destination directory')
    
    # parse the arguments
    args = parser.parse_args()

    # store into global values
    fastqList = args.input
    destRoot = args.dest
    
    countReads(fastqList, destRoot)

    fastqList.close()

# initialize the script
if __name__ == '__main__':
    sys.exit(main())
