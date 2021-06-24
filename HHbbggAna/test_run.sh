for dataset in `cat testfile/bbgg_nanoAOD_file_lists.txt`; do
    echo "get a list of dataset:", $dataset

    ./analyzeHHbbgg testfile/lists/$dataset.txt $dataset.root $dataset F F 2018
    #./analyzeHHbbgg runList.txt out.root mc F F 2018

done
