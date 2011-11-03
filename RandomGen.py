from random import randrange

numWords = 1000
numChunks = 25 #can only generate 702 chunks with current code, due to difficulties generating unique strings as chunk names. This should enougth for any testing though
maxWordCount = 10000

output = open("inputTest.tsv", 'w')

#create the first line which is a list of the words
line = "\t"
for i in range(0,numWords):
        line = line + str(i) + "\t"

line = line + "\n"
output.write(line)

#create each chunk
for i in range(0,numChunks):
        line = ""
        if(i > 26):
                line = chr(i / 26 + 65)
			
        line = line + chr(i % 26 + 65) + "\t"

        for j in range(0,numWords):
                n = randrange(0,maxWordCount,1)
                line = line + str(n) + "\t"

        line = line + "\n"
        output.write(line)

output.close()
