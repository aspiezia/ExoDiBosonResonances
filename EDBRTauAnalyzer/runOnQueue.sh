#!/bin/bash   

dataset=$1
numJob=$2
queue=$3


mkdir $dataset
cd $dataset
cp ../createConfigFile.sh .
cp ../createRunFile.sh .
cp ../data/$dataset\_fileList.txt .

x=$(cat $dataset\_fileList.txt | wc -l)
x=$(($x+1))
maxJob=$(($x / $numJob)) 
maxJob=$(($maxJob+1))

for ((i = 0; i < $maxJob ; i++)) ;
do
    min=$(($i*$numJob));
    max=$(($(($i+1))*$numJob));
    if [ "$min" -lt "$x" ]; then
	./createConfigFile.sh $min $max $i $dataset &> $dataset\_$i\_cfg.py
	./createRunFile.sh $dataset\_$i\_cfg.py $dataset &> run$i.sh
	chmod 777 run$i.sh
	echo "bsub -q 8nh run$i.sh"
	bsub -q $queue run$i.sh
    fi
done
