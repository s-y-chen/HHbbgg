NJOB=0
NCHUNK=5
year=2017
ref=2018
date=ntuple${year}0802v1
for dataset in `cat ../testfile/bbgg_nanoAOD_file_lists_${year}.txt`; do
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
    if [ $year -eq $ref ]
    then
        echo "cat ../../../testfile/skim_lists_${year}/${dataset}.txt > Job${NJOB}_list.txt"
    	cat ../../../testfile/skim_lists_${year}/${dataset}.txt > Job${NJOB}_list.txt
    else
        echo "cat ../../../testfile/lists_${year}/${dataset}.txt > Job${NJOB}_list.txt"
	cat ../../../testfile/lists_${year}/${dataset}.txt > Job${NJOB}_list.txt
    fi
    pwd > path.txt
    
    #Run N different random chunks per dataset
    echo "Preparing job chunks"
    split -l$NCHUNK Job${NJOB}_list.txt jobfiles_split.txt.

    #Split on line, not on space
    IFS=$'\n'
    ifile=0
    for f in `\ls -1 jobfiles_split.txt.*`; do
        echo "f: ", $f
        mv $f Job${NJOB}ifile${ifile}_list.txt
        jobname=job${NJOB}ifile${ifile}_${date}
        cat ../../example_job_tmp.jdl > submit_$jobname.jdl
        sed -i "s/IJOB/${NJOB}ifile${ifile}/g" submit_$jobname.jdl
        sed -i "s/NAME/$dataset/g" submit_$jobname.jdl
        sed -i "s/DATE/$data/g" submit_$jobname.jdl
        sed -i "s/YEAR/$year/g" submit_$jobname.jdl
        echo "Queue" >> submit_$jobname.jdl
        echo >> submit_$jobname.jdl
        condor_command="condor_submit submit_"$jobname".jdl"
        echo $condor_command
        condor_submit "submit_"$jobname".jdl"
    
        ifile=$((ifile + 1))
    done
    cd ../../
done
