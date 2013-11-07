####### TAKE THE CODE FORM GIT #######
cmsrel CMSSW_5_3_12
cd CMSSW_5_3_12
cmsenv
git cms-merge-topic -u cms-tau-pog:CMSSW_5_3_X_HighPt
cd src
cmsenv
git clone git@github.com:YOURNAME/ExoDiBosonResonances.git
cvs co -r V00-02-03s TauAnalysis/CandidateTools
scram b


###### RUN ON THE QUEUES #####
#1) Inside data/ create a txt file with name dataset_fileList.txt and put the list of root file that you want to process
#2) In createConfigFile.sh correct the python config file that you want to use
#3) Run on the queues with the following command:
./runOnQueue.sh dataset N queue

#where:
# - is the dataset name from dataset_fileList.txt
# - N gives the number of jobs to run through the following formula N(jobs) = N(root file)/N
# - queue is the queue on which you want to run 