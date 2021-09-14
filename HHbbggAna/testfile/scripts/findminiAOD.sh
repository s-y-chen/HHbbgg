NJOB=0
for dataset in `cat nanoAODlist_all.txt`; do
    #echo "nanoAOD dataset:", $dataset
    NJOB=$((NJOB + 1))
    dasgoclient --query="parent dataset=$dataset"
done
echo $NJOB