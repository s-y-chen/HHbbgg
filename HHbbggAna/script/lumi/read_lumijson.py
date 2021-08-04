# Python program to read
# json file
  
  
import json

def main():
  
	# Opening JSON file
	f = open('../../data/GoodRunJSON/Cert_271036-284044_13TeV_ReReco_07Aug2017_Collisions16_JSON.txt',)
  
	# returns JSON object as a dictionary
	data = json.load(f)
  
	# print the keys and values
	for key in data:
    		value = data[key]
    		print(key," ",value)
 
	# Closing file
	f.close()

if __name__ == '__main__':
    main()
