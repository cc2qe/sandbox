#!/bin/bash

#testarray=( $(ls .) )

# Source global definitions
if [ -f /etc/bashrc ]; then
	. /etc/bashrc
fi

nJobs=`bjobs -q normal -r -w | sed -e 1d | awk '{ print $7 }'`

#echo $nJobs

for i in $nJobs
do
	bswitch -J $i -q normal long;
done


