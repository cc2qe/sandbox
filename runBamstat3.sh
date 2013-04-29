if [ $# -eq 4 ]; then

BAMSTAT=/scr/talkowski/tools/code/bamstat/bamstat
CLUSTER=/scr/talkowski/tools/code/readPairCluster/readPairCluster

BAMFILE=$1
STATSDIR=$2
DISTANCE=$3
MAPQ=$4

echo Saving stats on `pwd`/$BAMFILE to `pwd`/$STATSDIR

mkdir $STATSDIR
cd $STATSDIR

$BAMSTAT -b -i ../$BAMFILE

sort -k3,3 -k6,6 -k4n,4n translocPairs.txt > translocPairs.sorted.txt
$CLUSTER -u -d $DISTANCE -s 3 -q $MAPQ -r translocPairs.sorted.txt > translocClusters_d${DISTANCE}_q${MAPQ}.txt

sort -k3,3 -k6,6 -k4n,4n inversionPairs.txt > inversionPairs.sorted.txt
$CLUSTER -u -d $DISTANCE -s 3 -q $MAPQ -r inversionPairs.sorted.txt > inversionClusters_d${DISTANCE}_q${MAPQ}.txt

sort -k3,3 -k6,6 -k4n,4n deletionPairs.txt > deletionPairs.sorted.txt
$CLUSTER -u -d $DISTANCE -s 3 -q $MAPQ -r deletionPairs.sorted.txt > deletionClusters_d${DISTANCE}_q${MAPQ}.txt

MAPQ=20

$CLUSTER -u -d $DISTANCE -s 3 -q $MAPQ -r translocPairs.sorted.txt > translocClusters_d${DISTANCE}_q${MAPQ}.txt
$CLUSTER -u -d $DISTANCE -s 3 -q $MAPQ -r inversionPairs.sorted.txt > inversionClusters_d${DISTANCE}_q${MAPQ}.txt
$CLUSTER -u -d $DISTANCE -s 3 -q $MAPQ -r deletionPairs.sorted.txt > deletionClusters_d${DISTANCE}_q${MAPQ}.txt



else
 echo "usage:"
 echo "  runBamstat.sh [bamfilename] [statsOuputDir] [distance] [mapQ]"
 echo ""
fi

