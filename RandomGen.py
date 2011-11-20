#try shuffling word counts
#change to maxWordCount 5  X
#set numWords to DAZ set maxWordCount to max of overall chunks in DAZ set numChunks to 11 
#check if transpose is working correctly


from random import randrange

numWords = 100
numChunks = 4 #can only generate 702 chunks with current code, due to difficulties generating unique strings as chunk names. This should enougth for any testing though
maxWordCount = 5

output = open("inputTest.tsv", 'w')

#create the first line which is a list of the words
line = ""
for i in range(0,numWords):
        line = line + "\t" + str(i)
		
line = line + "\n"
output.write(line)

#create each chunk
for i in range(0,numChunks):
        line = ""
        if(i > 26):
                line = chr(i / 26 + 65)
			
        line = line + chr(i % 26 + 65)

        for j in range(0,numWords):
                n = randrange(0,maxWordCount,1)
                line = line + "\t" + str(n)
				
        line = line + "\n"
        output.write(line)

output.close()
