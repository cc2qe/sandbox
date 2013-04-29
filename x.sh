RS=""
for i in {1..10}
do
	R=`echo ls $i | qsub -N x$i`
	RS="$RS,$R"
done
echo $RS
