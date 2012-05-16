#Code to randomly generate a dataset to run trueTree on.


from random import randrange

#numWords = 500 #number of worods to use in each chunk
#numChunks = 10 #number of chunks to generate. The code can only generate 702 chunks with current code, due to difficulties generating unique strings as chunk names. This should enough for any testing though
#minWordCount = 0
#maxWordCount = 6

def randomGen(fileName, numChunks, numWords, minWordCount, maxWordCount):
	output = open(fileName, 'w') #open a file to write to

	#create the first line which is a list of the words
	#each word is defined as being a number from 0 to numWords - 1
	line = "" #start with a blank line
	for i in range(0,numWords): 
			line = line + "\t" + str(i) #add a tab to go the the next column then add the word (The first column is supposed to be blank)
			
	line = line + "\n" #add the end line character
	output.write(line)

	#create each chunk
	for i in range(0,numChunks): #for each chunk
			line = ""
			
			#each chunk needs a unique name using only alphabetical characters. As such chunks are named A, B, C .... Z. If this is not enough the code will go on to
			#AA, AB, AC, ... AZ, BA, BB, ......
			
			if(i >= 676):
					line = line + chr((i / 26 / 26) + 64)
			
			if(i >= 26): #if i is greater then 26 the name is two letters long. The number of complete sets of 26 chunks determines the first letter. (If i is 26 first letter is A, if i is 52 first letter is B etc)
					line = line + chr((i / 26) % 26 + 65)
				
			line = line + chr(i % 26 + 65) #convert the chunk number into a letter. %26 is to seperate portion of number determining the first letter with part determing the second letter

			if(line == "NA"):
				line = "SPECIALFIX"
			
			for j in range(0,numWords): #for each word
					n = randrange(minWordCount,maxWordCount + 1,1) #pick a number between minWordCount and maxWordCount inclusive
					line = line + "\t" + str(n) #add tab to go to next column then add count for word to line
					
			line = line + "\n" #add end of line character
			output.write(line) 
		
	output.close()

randomGen("genomicsTest1000.tsv", 1000, 136, 0, 500)
