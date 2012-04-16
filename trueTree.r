source ( 'pvclust.R' )
source ( 'pvclust-internal.R' )

#library(pvclust)

# Idea for text labeling input
# textlabs is a list of character strings for the text part of the the label
# chunksize is a list of numeric for the number of chunks in each text
# textlabs and chunksize must have same length and are assumed to correspond elementwise
# so that first textlabs has first chunksize number of chunks

trueTree <- function(input.file, outputFilename = NULL, main = NULL, textlabs = NULL , chunksize = NULL , labelFileName = NULL, 
					distMetric = "euclidean" , clustMethod = "average" , use.cor="pairwise.complete.obs", input.transposed = TRUE, nboot = 100, runParallel = FALSE,
					numCPUs = 2, clusterType = 'SOCK', confidenceInterval = .95, seed = NULL, cladeChunkIn=NULL, storehClust=FALSE, storeChunks=FALSE, rowSample=FALSE, r=seq(.5,1.4,by=.1), 
					height=800, width=800, highlightFileName=NULL, metadata = FALSE, plotOut = TRUE, logFileName = "")
{
	#input.file is a character vector containing the name of the file to use as input.
	
	#outputFilename is a character vector used to name an output file. If outputFilename is NULL the program will display a plot of the dendrogram. If outputFilename is not null that plot will instead be saved as a png file named
	#outputFilename
	
	#main is a character vector that will be used as a title for the the dendrogram drawn. If main is left set to null no title will be added

	## List of possible distance metrics
	## METHODS <- c("euclidean", "maximum", "manhattan", "canberra",
	## "binary", "minkowski", "correlation", "uncentered", "abscor")

	## List of possible cluster-distance methods
	## METHODS <- c("ward", "single", "complete", "average", "mcquitty",
	## "median", "centroid")

	## If the input has text names on top, set input.transposed <- FALSE

	#nboot is the number of bootstrap replications to preform. Since pvclust does multiscale bootstrapping it generates this number of bootstraps 10 different times using different sample sizes
	#so the total number of new hclust objects generated is 10*nboot
	
	#runParallel sets whether the program does the bootstrapping on a single core or runs it in parallel. If runParallel is false the code is executed on a single processer and there are no special requirements for 
	#being able to run the code. The numCPUs, and clusterType parameters will also be ignored. If runParallel is set to true, the snow and snowfall librarys must be installed and the numCPUs,
    #and clusterType parameters will be used. Currently there is one other discrepency that needs to be fixed. When runParallel is false each run will produce a new set of results, and when runParallel is true, each run will
	#produce the same set of results due to parpvclust giving each processor a fixed seed.
	
	#clusterType is the type of cluster to create. If set to 'SOCK' the cluster will use raw sockets, which due to security concerns will not work over a network with some computers. If set to MPI will use the MPI protocal
	#which requries that the computer has a MPI implementation, installed, properly set up and running, as well as having the r package rmpi installed. PVM (Parallel Virtual Machine) and NWS(networkspaces) 
	#are theoretically supported but not tested. See the snow/snowfall documentation for more information about the necessary steps to use those.
	
	#numCPUs is the number of processers to use in the cluster. If the clusterType is set to SOCK all the processers used will be from the machine running this code. If clusterType is set to MPI the processers used 
	#can come from any machine mpi is set up to access. If numCPUs is set to a value greater then the number of processers actually available the code will still run, but some processers will end up running multiple instances
	#of the program which will increase the total run time.
	
	#confidence interval is the confidence interval for the cophenetic correlations given as a percent entered in decimal form. The program will give the output the cophenetic correlations at the lower and upper ends of the interval. 
	#The lower bound can be computed as (100% - confidence interval) / 2 and the upper bound as (100% - lower bound). For example if the confidence interval is 95% the lowerbound will be 2.5% and the upper bound 97.5%

	#the seed parameter is used to set the seed used for random number generation. By default seed is set to null which will cause the program to use a random seed. If a value for seed is given the program will instead use that
	#value as the seed for random number generation. When not running in parallel a single value should be entered for seed. When running in parallel a list the same length as numCPUs should be entered which each item in the list
	#being the seed for one particular processor. If seed is set to null and you want to check what the number generated for use as a seed was, it is stored as the attribute seed of the pvclust object returned
	#For example if you named the pvclust object returned result the new chunks can be accessed as follows: result$seed.
	
	#cladeChunkIn is a parameter to allow resampling by clades. Normally the resampling happens on a per chunk basis, meaning that the words in the new chunks created from the bootstraping are taken from the list of words in that chunk
	#when resampling is done by clades each chunk inside of one of the clades used in the resampling draws from all the words in the clade when resampling is perfomed.  cladeChunkIn allows the user to define the clades being used. 
	#It is a list the same length as the number of chunks in the input file. The first number in the list corrosponds to the first chunk and so on with the first chunk being the first chunk listed in the input file. Each number
	#specifies which clade that chunk belongs in. For instance if the list passed was [1,2,2,3] the first chunk goes into clade 1, chunks 2 and 3 go into clade 2, and chunk 4 goes into clade 3.
	#If this parameter is left to null the program will resample intraclade instead.
	#WARNING. It is up to the user to make sure the clades defined by this parameter make sense. It is quite possible to input clades that don't exist in the actual dendrogram which may produce weird results. T\
	
	#storehClust is a boolean that sets pvclusts store parameter. If storehClust is set to true then the result objected returned by the function will contain a storehClust attribute which is a list holding all the hclust objects created during
	#the bootstrapping. The hclust obects will be stored in a list as the attribute storehClust of the of the pvclust object returned
	#For example if you named the pvclust object returned result the new chunks can be accessed as follows: result$storehClust. The format of storehClust is as a list containing a number of sublists with there being one sublist for each
	#r value pvclust used, and each sublist having all the hclust objects created for that r value
	#WARNING: setting storehClust to true will likely cause R to run out of memory if nboot is set to a high value like 100,000. As such this parameter should be set to false when preforming actual validation
	#and only used on smaller test runs to see what kind of dendrograms are being created. 
	
	#storeChunks is a boolean that sets pvclusts store parameter. If store is set to true then the result objected returned by the function will contain a store attribute which is a list holding all the new chunks created during
	#the bootstrapping. The new chunks	will be stored in a list as the attribute storeChunk of the pvclust object returned
	#For example if you named the pvclust object returned result the new chunks can be accessed as follows: result$storeChunks. The format of storeChunks is as a list containing
	#a number of sublists with there being one sublist for each r value pvclust used, and each sublist having all the new chunks created for that r value
	#WARNING: setting store to true will likely cause R to run out of memory if nboot is set to a high value like 100,000. As such this parameter should be set to false when preforming actual validation
	#and only used on smaller test runs to see what kind of dendrograms are being created. 
	
	#rowSample is a boolean that if true allows the use of pvclusts original resampling method, where instead of basing the resampling around the number of times each individual word appeared in each individual clade, it builds
	#the new data table by picking out a new set of words to use for all the chunks keeping the counts for the number of times those words appeared in the orignal text for those chunks. If rowSample is set to true then the 
	#cladeNumber parameter will be ignored. WARNING. It is currently unclear whether this is actually a good method of resampling. By default this parameter is set to false and for now it is advised that it be left to false
	#unless you have a specific reason for wanting to use it.
	
	#r lets you set the r values that will be used in the multiscale bootstrapping. The format is a vector/list containing all the r values to use
	
	#height and width set the height and width of the dendrogam plottted, when the dendrogram is saved to a file
	
	#highlightFileName is the name of a file containing of a list of chunks to highlight in the dendrogram. The format is one chunk per line, and each line contains the exacat label for the chunk that should be
	#highlighted
	
	#metadata is a boolean that states whether the file contains metadata along with the counts for different words. Currently this only works with one specific metadata format. 
	
	#This function will return a list containing a modified pvclust object with a number of additional attributes. The storeCop attribute is a list of lists with each sublist containing the cophenetic correlation values
	#for all the hclust objects generated for a particular r value. The first sublist contains the values for the first r value used and so on. There is also a seed attribute and a storeChunks attributes which are described
	#in the documentation for the parameters with the same names. If you want to see what the r values where they are stored as a list in the attribute r
	
	library(stats)
	Sys.time()->startTotal; #holds start time of program so time to run entire program can be calculated
	Sys.time()->startSection; #start timing the bootstraping
	
	#set inital values for certain variables
	storeData <- NULL #
	
	#change this for the text you'd like to input
	input.data <- read.table(as.character(input.file), header=T, comment.char="", row.names=1, sep="\t", quote="")

	metaTable <- NULL #assume there is no metadata
		
	#if there is metadata
	if(metadata)
	{
		metaTable <- input.data[,-(1:136)] #put metadata into a seperate table
		input.data <- input.data[,1:136] #keep the non metadata in the original table
	}
	
	if ( input.transposed ) #if the input
		tTable <- t(input.data)
	else
		tTable <- ( input.data )


	if( !is.null(textlabs) && !is.null(chunksize)) 
	{
		if(length(textlabs) != length(chunksize)) 
		{
				stop("number of texts and corresponding chunk numbers must match")
		}
		else 
		{		# check that sum(chunksize) == dim(relFreq)[1] , total number of chunks equals number of rows in relFreq
				L <- length(chunksize)
				temp <- NULL
				for(i in 1:L) 
				{
					for(k in 1:chunksize[i])
					{
						temp <- c(temp,paste(textlabs[i],as.character(k),sep=""))
					}
				}
				#print(temp)
				colnames(tTable) <- temp
		}
	}
	
	if(!is.null(labelFileName)) #if a file of label names to highlight was specified
	{
		chunkLabels <- scan(labelFileName, what = "character", sep = "\n")
		colnames(tTable) <- chunkLabels
	}


	copValues <- numeric(0)

	write("Inital runtime:", file = logFileName, append=TRUE)
    write(format(Sys.time()-startSection, usetz=TRUE), file = logFileName, append=TRUE);
	Sys.time()->startSection; #start timing the bootstraping
	
	if(!runParallel)
	{
		pCluster <- pvclust(tTable, nboot=nboot, method.hclust=clustMethod, method.dist=distMetric, use.cor = use.cor, storeCop=TRUE, seed=seed, cladeChunkIn=cladeChunkIn, store=storehClust, storeChunks=storeChunks, rowSample=rowSample, r=r)
	}

	else
	{
		#ADD MORE ERROR CHECKING
		library(snowfall) #I will write more about next line missing ip argument
		
		sfInit(parallel=TRUE, cpus=numCPUs,type=clusterType) #This creates the cluster. If clusterType is SOCK all the processors will be taken from the current computer. If the type is MPI, the processors can be taken
		                                                           #from any computer running an MPI implementation, and the MPI program will handle choosing what processors to use
																   
		if(!sfIsRunning())
		{
			write("Error. Cluster not created successfully. Will use nonparallel version of pvclust instead")
			pCluster <- pvclust(tTable, nboot=nboot, method.hclust=clustMethod, method.dist=distMetric, use.cor = use.cor, storeCop=TRUE, seed=seed, cladeChunkIn=cladeChunkIn, store=storehClust, storeChunks=storeChunks, rowSample=rowSample, r=r)
		}
		cl <- sfGetCluster()
		sfExportAll()
		## parallel version of pvclust
		pCluster <- parPvclust(cl,tTable, nboot=nboot, method.hclust=clustMethod, method.dist=distMetric, use.cor = use.cor, storeCop=TRUE, normalize=TRUE, seed=seed, cladeChunkIn=cladeChunkIn, store=storehClust, storeChunks=storeChunks, rowSample=rowSample, r=r)
	}
	
	#add the metadata table
	pCluster <- c(pCluster, metaTable = list())
	pCluster$metaTable[[1]] <- metaTable #if we tried to add metaTable directly to pCluster, the c command would break it into a number of different componenets
	class(pCluster) <- "trueTree" #set object to right class

	write("Bootstrap runtime:", file = logFileName, append=TRUE)
    write(format(Sys.time()-startSection, usetz=TRUE), file = logFileName, append=TRUE);
	Sys.time()->startSection; #start timing the analysis of the cophenetic correlations

	#find cophenetic correlation of orignal hclustering	
	orginalCopDist <- cophenetic(pCluster$hclust) #get cophenetic distance matrix for original hclust object
	
	originalCor <- cor(orginalCopDist, pCluster$distance) #get correlation with original distance matrix

	#the cophenentic correlations are stored in a series of lists which are themselves stored in a list, so merge them all into one list
	listH <- pCluster$storeCop

	copValues <- unlist(listH)


	if(runParallel)
	{
		#copValues <- parallelSort(copValues)
	}
	
	#write(copValues)
	copValues <- unlist(copValues)
	copValues <- sort(copValues) #sort the list
	copSize <- length(copValues)
	
	lowerBound = (1-confidenceInterval) / 2 #to get the lower bound subtract the confidence interval from 1 to get the precentage outside the interval and divide by 2 to get the precentage below the interval
	upperBound = 1 - lowerBound #subtract the lowerbound from 100% to get upper end of ranger
	
	write(paste("Original Cophenetic correlation", originalCor, sep=" "), file = logFileName, append=TRUE)
	write(paste("Number of Cophenetic correlation values", length(copValues), sep=" "), file = logFileName, append=TRUE)
	write(paste("Minimum Cophenetic correlation", min(copValues), sep=" "), file = logFileName, append=TRUE)
    write(paste(lowerBound * 100, "% interval Copheneticcorrelation ", copValues[round(copSize * lowerBound)], sep=""), file = logFileName, append=TRUE)	
    write(paste("Median Cophenetic correlation", median(copValues), sep=" "), file = logFileName, append=TRUE)
	write(paste(upperBound * 100, "% interval Cophenetic correlation ", copValues[round(copSize * upperBound)], sep=""), file = logFileName, append=TRUE)	
    write(paste("Maximum Cophenetic correlation", max(copValues), sep=" "), file = logFileName, append=TRUE)
	
	write(paste("Mean Cophenetic correlation", mean(copValues), sep=" "), file = logFileName, append=TRUE)
    write(paste("Standard Deviation Cophenetic correlation", sd(copValues), sep=" "), file = logFileName, append=TRUE)
    
	write("Cophenetic analysis runtime:", file = logFileName, append=TRUE)
    write(format(Sys.time()-startSection, usetz=TRUE), file = logFileName, append=TRUE);
	
	#find range between 2.5 % in and 97.5 % sorted make parameters use order/sort
	
	specialLabels = NULL
	
	if(plotOut)
	{
		if(!is.null(highlightFileName)) #if a file of label names to highlight was specified
		{
			specialLabels <- scan(highlightFileName, what = "character", sep = "\n")
		}
		
		if(!is.null(outputFilename))
		{
			plot(pCluster, filename=outputFilename, main = main, height=height, width=width, specialLabels=specialLabels)
		}
		
		else
		{
			plot(pCluster, main = main, specialLabels=specialLabels)
			par(ask=TRUE)
		}
		
		#plot a histogram showing distribution of cophenetic correlation values)
		
		if(!is.null(outputFilename)) #if the dendrogram is being output to file, output the histogram as well
		{
			  if(.Platform$OS.type == "windows")
			  {
					png(paste(outputFilename, "histogram.png", sep=""), width=800, height=800)
			  }
			  
			  else if(.Platform$OS.type == "unix")
			  {
					png(paste(outputFilename, "histogram.png", sep=""), width=800, height=800, type="Xlib")
			  }
		}
		
		hist(copValues)
		abline(v = originalCor, col = "red")
		
		if(!is.null(outputFilename)) #if outputing to a file close the file
		{
			dev.off()
		}
	}
	
	write("Total runtime:", file = logFileName, append=TRUE)
    write(format(Sys.time()-startTotal, usetz=TRUE), file = logFileName, append=TRUE);

	if(runParallel)
	{
		sfStop() #end the cluster
	}
	
	#pCluster <- c(pCluster, metaTable)
	
	return(pCluster)
}

varianceTest <- function(input.file, distMetric = "euclidean" , clustMethod = "average" , input.transposed = TRUE, nboot = 10, runParallel = FALSE,
					numCPUs = 2, clusterType = 'SOCK', cladeChunkIn=NULL, rowSample=FALSE, r=seq(.5,1.4,by=.1), 
					metadata = FALSE, testRuns = 5)
{
		#do first test run
		result <- myCluster(input.file, distMetric = distMetric, clustMethod = clustMethod, input.transposed = input.transposed, nboot = nboot, runParallel = runParallel,
		numCPUs = numCPUs, clusterType = clusterType, cladeChunkIn = cladeChunkIn, rowSample = rowSample, r = r, metadata = metadata, plotOut = FALSE, logFileName = "dummy.txt")
		
		auvalues <- result$edge[,"au"]
		bpvalues <- result$edge[,"bp"]

		#create matrice for AU and BP values
		allAU <- matrix(byrow=T, ncol=length(auvalues), nrow = testRuns)
		allBP <- matrix(byrow=T, ncol=length(auvalues), nrow = testRuns)
		
		allAU[1,] <- auvalues
		allBP[1,] <- bpvalues

		#do rest of test runs
		for(i in 2:testRuns)
		{
			result <- myCluster(input.file, distMetric = distMetric, clustMethod = clustMethod, input.transposed = input.transposed, nboot = nboot, runParallel = runParallel,
			numCPUs = numCPUs, clusterType = clusterType, cladeChunkIn = cladeChunkIn, rowSample = rowSample, r = r, metadata = metadata, plotOut = FALSE, logFileName = "dummy.txt")
			
			auvalues <- result$edge[,"au"]
			bpvalues <- result$edge[,"bp"]

			#create matrice for AU and BP values
			
			allAU[i,] <- auvalues
			allBP[i,] <- bpvalues
		}
		
		minAUs <- apply(allAU, 2, min) #find the minimum au value in each column
		maxAUs <- apply(allAU, 2, max) #find the maximum au value in each column
		
		auDiffs <- (maxAUs - minAUs) * 100 #get difference between min and max au values and convert to percentages (values were decimals)
		
		minBPs <- apply(allBP, 2, min) #find the minimum bp value in each column
		maxBPs <- apply(allBP, 2, max) #find the maximum bp value in each column
		
		bpDiffs <- (maxBPs - minBPs) * 100 #get difference between min and max au values and convert to percentages (values were decimals)
		
		print(paste("Greatest AU Difference", max(auDiffs)), sep=" ")
		#print(paste("Greatest BP Difference", max(bpDiffs)), sep=" ")
		
}

#print("10")
#varianceTest("danile-azarius.txt", nboot = 10, input.transposed = FALSE, runParallel = TRUE, testRuns = 20)
#print("100")
#varianceTest("danile-azarius.txt", nboot = 100, input.transposed = FALSE, runParallel = TRUE, testRuns = 20)
#print("500")
#varianceTest("danile-azarius.txt", nboot = 500, input.transposed = FALSE, runParallel = TRUE, testRuns = 20)
#print("1000")
#varianceTest("danile-azarius.txt", nboot = 1000, input.transposed = FALSE, runParallel = TRUE, testRuns = 20)
#print("5000")
#varianceTest("danile-azarius.txt", nboot = 5000, input.transposed = FALSE, runParallel = TRUE, testRuns = 20)
#print("10000")
#varianceTest("danile-azarius.txt", nboot = 10000, input.transposed = FALSE, runParallel = TRUE, testRuns = 20)
#print("20000")
#varianceTest("danile-azarius.txt", nboot = 20000, input.transposed = FALSE, runParallel = TRUE, testRuns = 20)

logFile <- ""
#write("3000 \n", file = logFile)
result <- trueTree("fedpapers.tsv", outputFilename = "dummy.png", nboot=5, distMetric = "euclidean", runParallel = TRUE, input.transposed = TRUE, numCPUs = 2, clusterType = "SOCK", plotOut=FALSE, logFileName = logFile, metadata=FALSE)
#write("4000 \n", file = logFile, append=TRUE)
#result <- myCluster("genomicsTest4000.tsv", outputFilename = "Dummy.png", nboot=2, main="Genomics Data", distMetric = "euclidean", runParallel = TRUE, input.transposed = TRUE, numCPUs = 2, clusterType = "SOCK", height = 3000, width = 100000, plotOut=FALSE, logFileName = logFile, metadata=FALSE)
#write("5000 \n", file = logFile, append=TRUE)
#result <- myCluster("genomicsTest5000.tsv", outputFilename = "Dummy.png", nboot=2, main="Genomics Data", distMetric = "euclidean", runParallel = TRUE, input.transposed = TRUE, numCPUs = 2, clusterType = "SOCK", height = 3000, width = 100000, plotOut=FALSE, logFileName = logFile, metadata=FALSE)
#write("6000 \n", file = logFile, append=TRUE)
#result <- myCluster("genomicsTest6000.tsv", outputFilename = "Dummy.png", nboot=2, main="Genomics Data", distMetric = "euclidean", runParallel = TRUE, input.transposed = TRUE, numCPUs = 2, clusterType = "SOCK", height = 3000, width = 100000, plotOut=FALSE, logFileName = logFile, metadata=FALSE)


#result <- myCluster("GT.txt", outputFilename = "ABCDEF", nboot=1, main="Genomics Data", distMetric = "euclidean", runParallel = FALSE, input.transposed = TRUE, numCPUs = 2, clusterType = "SOCK", height = 3000, width = 100000, plotOut=TRUE, logFileName = logFile, metadata=TRUE)
#textlabs = c("Az", "Dan"), chunksize = c(1,10))
