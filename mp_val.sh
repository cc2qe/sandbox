#!/bin/bash

############################################################
#  Program:
#  Author :
############################################################


## BEGIN SCRIPT
usage()
{
    cat << EOF

usage: $0 OPTIONS

OPTIONS can be:
    -h      Show this message
    -i      Ira file name
    -s      slop (default: 0)
    -m      Moleculo file (Defualt /mnt/hall13_local/cc2qe/na12878_pacbio...)
    -p      Pacbio file (Default /mnt/hall13_local/cc2qe/na12878_moleculo...)

EOF
}

# Show usage when there are no arguments.
if test -z "$1"
then
    usage
    exit
fi

VERBOSE=
FILENAME=

PB="/mnt/hall13_local/cc2qe/na12878_pacbio/NA12878.pacbio.splitreads.excldups.breakpoint.bedpe.gz"
MO="/mnt/hall13_local/cc2qe/na12878_moleculo/NA12878.moleculo.splitreads.excldups.breakpoint.bedpe.gz"
L=
SLOP=0

# Check options passed in.
while getopts "h m:p:i:s:" OPTION
do
    case $OPTION in
        h)
            usage
            exit 1
            ;;
        m)
            MO=$OPTARG
            ;;
        p)
            PB=$OPTARG
            ;;
        i)
            L=$OPTARG
            ;;
	s)
	    SLOP=$OPTARG
	    ;;
        ?)
            usage
            exit
            ;;
    esac
done

pairToPair -type both -is -slop $SLOP -a $PB -b $L \
    | awk '$18==$29' \
    | zapdups -u \
    > $L.p.slop$SLOP.tmp

pairToPair -type both -is -slop $SLOP -a $MO -b $L \
    | awk '$18==$29' \
    | zapdups -u \
    > $L.m.slop$SLOP.tmp

cat $L.p.slop$SLOP.tmp \
    | cut -f 19- \
    | sort -k7,7n \
    | groupBy -g 1,2,3,4,5,6 -c 1 -o count -full \
    | zjoin -r -a $L -b stdin -1 7 -2 7 \
    | cut -f -22,45 \
    | awk '{ if ($23=="NA") { $23=0 } print $0 }' OFS="\t" \
    | zjoin -r -1 7 -2 7 -a stdin -b <(cat $L.m.slop$SLOP.tmp | cut -f 19- | sort -k7,7n | groupBy -g 1,2,3,4,5,6 -c 1 -o count -full) \
    | cut -f -23,46 \
    | awk '{ if ($24=="NA") { $24=0 } print $0 }' OFS="\t" \
    > $L.slop$SLOP.val

rm $L.p.slop$SLOP.tmp $L.m.slop$SLOP.tmp
