#!/bin/bash   

dataset=$1
numJob=$2
queue=$3
analysis=$4
folder=$5


mkdir $folder
cd $folder
cp ../script/createConfigFileFullyLep.sh .
cp ../script/createConfigFileEleID.sh .
cp ../script/createConfigFileJetID.sh .
cp ../script/createConfigFilePt.sh .
cp ../script/createConfigFileMuTau.sh .
cp ../script/createConfigFileNoCleanTau.sh .
cp ../script/createConfigFileCleanTau.sh .
cp ../script/createConfigFileNoCleanETau.sh .
cp ../script/createConfigFileCleanETau.sh .
cp ../script/createRunFile.sh .
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
	if [ "$analysis" == "fullyLeptonic" ]; then
	    ./createConfigFileFullyLep.sh $min $max $i $dataset &> $dataset\_$i\_cfg.py
	fi
	if [ "$analysis" == "eleID" ]; then
	    ./createConfigFileEleID.sh $min $max $i $dataset &> $dataset\_$i\_cfg.py
	fi
	if [ "$analysis" == "jetID" ]; then
	    ./createConfigFileJetID.sh $min $max $i $dataset &> $dataset\_$i\_cfg.py
	fi
	if [ "$analysis" == "ptStudy" ]; then
	    ./createConfigFilePt.sh $min $max $i $dataset &> $dataset\_$i\_cfg.py
	fi
	if [ "$analysis" == "MuonTau" ]; then
	    ./createConfigFileMuTau.sh $min $max $i $dataset &> $dataset\_$i\_cfg.py
	fi
	if [ "$analysis" == "NoCleanTau" ]; then
	    ./createConfigFileNoCleanTau.sh $min $max $i $dataset &> $dataset\_$i\_cfg.py
	fi
	if [ "$analysis" == "CleanTau" ]; then
	    ./createConfigFileCleanTau.sh $min $max $i $dataset &> $dataset\_$i\_cfg.py
	fi
	if [ "$analysis" == "NoCleanETau" ]; then
	    ./createConfigFileNoCleanETau.sh $min $max $i $dataset &> $dataset\_$i\_cfg.py
	fi
	if [ "$analysis" == "CleanETau" ]; then
	    ./createConfigFileCleanETau.sh $min $max $i $dataset &> $dataset\_$i\_cfg.py
	fi
	./createRunFile.sh $dataset\_$i\_cfg.py $folder &> run$i.sh
	chmod 777 run$i.sh
	echo "bsub -q 8nh run$i.sh"
	bsub -q $queue run$i.sh
    fi
done
