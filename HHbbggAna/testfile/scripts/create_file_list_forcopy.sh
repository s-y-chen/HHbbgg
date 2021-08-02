for dataset in `cat nanoAODlist_data.txt`; do
	echo "get a list of dataset:", $dataset
	dasgoclient --query="file dataset=$dataset" >> file.txt
done
