#try shuffling word counts
#change to maxWordCount 5  X
#set numWords to DAZ set maxWordCount to max of overall chunks in DAZ set numChunks to 11 
#check if transpose is working correctly


from random import randrange

numWords = 500
numChunks = 10 #can only generate 702 chunks with current code, due to difficulties generating unique strings as chunk names. This should enougth for any testing though
minWordCount = 0
maxWordCount = 6

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
                n = randrange(minWordCount,maxWordCount + 1,1)
                line = line + "\t" + str(n)
				
        line = line + "\n"
        output.write(line)

i = 10
line = ""
if(i > 26):
		line = chr(i / 26 + 65)
	
line = line + chr(i % 26 + 65)

for j in range(0,250):
		n = randrange(1,maxWordCount + 1,1)
		line = line + "\t" + str(n)
		
for j in range(250,500):
		line = line + "\t" + str(0)
		
line = line + "\n"
output.write(line)
		
output.close()

