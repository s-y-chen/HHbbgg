for dataset in `cat testfile/bbgg_nanoAOD_file_lists_2018.txt`; do
    echo "get a list of dataset:", $dataset

    #./analyzeHHbbgg_w testfile/lists_2018/$dataset.txt ${dataset}_sumOfgenW.root $dataset F T 2018
    ./analyzeHHbbgg testfile/skim_lists_2018/$dataset.txt $dataset.root $dataset F T 2018

    #test command to run 2018 SM HH signal 
    #./analyzeHHbbgg testfile/skim_lists_2018/GluGluToHHTo2B2G_node_cHHH1_TuneCP5_PSWeights_13TeV-powheg-pythia8.txt out.root GluGluToHHTo2B2G_node_cHHH1_TuneCP5_PSWeights_13TeV-powheg-pythia8 F T 2018

    #test command to run 2017 SM HH signal 
    #./analyzeHHbbgg testfile/lists_2017/GluGluToHHTo2B2G_node_cHHH1_TuneCP5_PSWeights_13TeV-powheg-pythia8.txt out.root GluGluToHHTo2B2G_node_cHHH1_TuneCP5_PSWeights_13TeV-powheg-pythia8 F T 2017

    #test command to run 2016 SM HH signal
    #./analyzeHHbbgg testfile/lists_2016/GluGluToHHTo2B2G_node_cHHH1_TuneCUETP8M1_PSWeights_13TeV-powheg-pythia8.txt out.root GluGluToHHTo2B2G_node_cHHH1_TuneCUETP8M1_PSWeights_13TeV-powheg-pythia8 F T 2016
done
