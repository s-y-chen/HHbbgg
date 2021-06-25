for dataset in `cat testfile/bbgg_nanoAOD_file_lists.txt`; do
    echo "get a list of dataset:", $dataset

    ./analyzeHHbbgg_w testfile/lists/$dataset.txt ${dataset}_sumOfgenW.root $dataset F T 2018
    #./analyzeHHbbgg testfile/skim_lists/$dataset.txt $dataset.root $dataset F T 2018 
    #./analyzeHHbbgg testfile/skim_lists/GluGluToHHTo2B2G_node_cHHH1_TuneCP5_PSWeights_13TeV-powheg-pythia8.txt out.root GluGluToHHTo2B2G_node_cHHH1_TuneCP5_PSWeights_13TeV-powheg-pythia8 F T 2018

done
