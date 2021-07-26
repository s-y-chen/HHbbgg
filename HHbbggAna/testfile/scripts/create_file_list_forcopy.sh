for dataset in `cat nanoAODlist.txt`; do
	echo "get a list of dataset:", $dataset
	dasgoclient --query="file dataset=$dataset" >> file.txt
done
