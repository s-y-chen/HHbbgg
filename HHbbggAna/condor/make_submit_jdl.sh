NJOB=0
date=ntuple0625v1
for dataset in `cat ../testfile/bbgg_nanoAOD_file_lists.txt`; do
    echo "get a list of dataset:", $dataset
    NJOB=$((NJOB + 1))

    echo "job index: "$NJOB

    rm -Rf output/job_${NJOB}_${date}
    mkdir output/job_${NJOB}_${date}
    mkdir output/job_${NJOB}_${date}/logs
    cd output/job_${NJOB}_${date}
    #switch
    cp ../../bin/analyzeHHbbgg .
    cp ../../example_job.sh .
    echo "cat ../../../testfile/skim_lists/${dataset}.txt > Job${NJOB}_list.txt"
    #switch
    cat ../../../testfile/skim_lists/${dataset}.txt > Job${NJOB}_list.txt
    pwd > path.txt
    jobname=job${NJOB}_${date}
    cat ../../example_job_tmp.jdl > submit_$jobname.jdl
    sed -i "s/IJOB/$NJOB/g" submit_$jobname.jdl
    sed -i "s/NAME/$dataset/g" submit_$jobname.jdl
    sed -i "s/DATE/$data/g" submit_$jobname.jdl
    echo "Queue" >> submit_$jobname.jdl
    echo >> submit_$jobname.jdl
    condor_command="condor_submit submit_"$jobname".jdl"
    echo $condor_command
    condor_submit "submit_"$jobname".jdl"
    cd ../../

done
