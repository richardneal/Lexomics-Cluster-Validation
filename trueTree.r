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
#The main function in trueTree that is used to perform cluster validation. trueTree will output
#a dendrogram as well as an analysis of the cophenetic correlation coefficients from the bootstrap
#process and return an object that stores all of the important information for later use.

# input.file - Character vector containing the name of the file containing the input data

# outputFilename - Character vector containing the name of the file to save the dendrogram
# showing the results. The dendrogram will be saved as a .png file but the file extension should
# be left out of the file name as trueTree will add it automatically. Additionally a histogram
# showing the distribution of cophenetic correlation coefficients will be saved to outputFilenamehistogram.
# png. If outputFilename is set to NULL, which it is by default, R will instead create a
# window to display the results.

# main - Optional character vector containing a title for the dendrogram. This title will be
# displayed at the top of the dendrogram.

# textLabs - See description for chunksize

# chunkSize - textLabs and chunkSize are a pair of optional vectors that are used to relabel all
# the chunks in the dataset. These parameters are designed to be used when there the dataset
# contains a number of texts which have been split into multiple chunks. textLabs is a character
# vector and should contain the names of all the texts in the order in which they appear in the
# dataset. chunksize is a vector that should contain then number of chunks in each of those texts
# in the same order as the texts were listed. The function will then give the initial chunk the label
# for first chunk in the first text, the second chunk as being the second chunk in the first text and
# so forth until the number of chunks listed in chunkSize has been seen, at which point it will start
# labeling chunks as being part of the second text. For instance if there were five chunks total and
# textLabs = ("Bob", "George") and chunkSize = (2,3) the chunks will be labeled Bob 1, Bob 2,
# George 1, George 2, George 3. textLabs and chunkSize need to have the same length and the
# sum of the chunksizes needs to be equal to the total number of chunks.

# distMetric - Character vector containing the method to use when calculating distances between
# chunks. Supports "euclidean", "maximum", "manhattan", "canberra", "binary", "minkowski",
# "correlation", "uncentered", and "abscor" methods.

# clustMethod - Character vector containing the method to use when calculating distances
# between chunks. Supports "euclidean", "maximum", "manhattan", "canberra", "binary",
# "minkowski", "correlation", "uncentered", and "abscor" methods.

# use.cor - correlation method to use with the "correlation", "uncentered", and "abscor" distance methods
#supports "all.obs", "complete.obs" and "pairwise.complete.obs".

# input.transposed - Boolean that tells trueTree what input format to expect. input.transposed
# should equal TRUE if each row holds a chunk and each column a word, or FALSE if each column
# holds a chunk and each row a word.

# nboot - Integer that says how many resampled datasets should be generated for each r value

# runParallel - Boolean that whether to run the bootstrap in parallel. If TRUE the bootstrap
# will be run over multiple processors with the details of how specified by the numCPUs and
# clusterType parameters. If FALSE then the bootstrap will be run on a single processor.

# numCPUs - Integer denoting how many processors should the code be run over. This parameter
# has no effect if runParallel is set to false.

# clusterType - Character vector containing what type of machine cluster is being used. Valid
# types are "sock" for raw sockets, "pvm" for parallel virtual machine, "nws" for NetWorkSpaces
# and "mpi" for Message Passing Interface. "sock" is the only option that can be used without any
# initial setup but only allows the bootstrap to be split across the cores in the machine the code
# is being run on. "mpi" supports running the bootstrap across multiple machines but requires
# an existing mpi runtime. "pvm" and "nws" are untested. See the snow documentation Tierney
# et al. (2011) for more information on how to use these options.

# confidenceInterval - Numeric value containing is the confidence interval to use for the cophenetic
# correlation coefficient analysis. The function will return the cophenetic correlation coefficients
# that lie at the upper and lower bounds of this interval. The interval should be given as a
# decimal not a percentage.

# seed - Parameter that gives a seed to use for the random number generator. If runParallel is
# FALSE, then seed should be a single integer. If runParallel is true then seed should be a vector
# containing one seed for each processor with each seed being an integer. If seed is left to NULL
# then a seed or seeds will be automatically generated based on the current time.

# cladeChunkIn - Vector that allows the use of intra-clade based resampling. cladeChunkIn
# should contain 1 integer for each chunk, with the integer giving the number of the clade that
# chunk is in. If cladeChunkIn is left to NULL then intra-chunk based resampling will be used
# instead.

# storehClust - Boolean specifying whether or not to store all the generated hclust objects from
# the bootstrap. If TRUE the hclust objects will be saved, if FALSE they won't. WARNING:
# Storing the hclust objects consumes a lot of memory. If this option is used in conjunction with
# a high nboot value, or a large number of r values, the computer running trueTree may run out
# of RAM, causing the code to crash.

# storeChunks - Boolean specifying whether or not to store all the resampled datasets from the
# bootstrap. If TRUE the resampled datasets will be saved, if FALSE they won't. WARNING:
# Storing the resampled datasets consumes a lot of memory. If this option is used in conjunction
# with a high nboot value, or a large number of r values, the computer running trueTree may run
# out of RAM, causing the code to crash.

# rowSample - Boolean specifying whether or not to use pvclusts original resampling method. If
# TRUE interrow based resampling will be used, where rows are picked from the original dataset
# to place in the new resampled dataset. WARNING: This resampling method is poorly suited to
# some domains like text mining. Before this option is used the user should determine whether it
# is appropriate for the domain where it is being applied.

# r - Vector specifying what r values should be used to run the multiscale bootstrapping. r must
# contain at least two values for any AU values to be generated.

# height - Vector specifying the height in pixels of the dendrogram that will be created. This
# parameter only works when outputting to a file.

# width - Vector specifying the width in pixels of the dendrogram that will be created. This
# parameter only works when outputting to a file.

# highlightFileName - Character vector containing the name of a file containing a list of chunks
# that should be highlighted in the dendrogram. The file should list one chunk per line, and the
# line should contain the exact label of the chunk to highlight.

# metadata - Boolean specifying whether the input file contains metadata. Metadata is used to
# color a dendrogram. Currently only one metadata format is recognized. See Appendix A for
# details.

# plotOut - Boolean that specifies whether the dendrogram containing the results of the cluster
# validation should be displayed. If true the dendrogram will be displayed/saved to a file. If false
# no action will be taken.

# logFileName - Optional character vector containing the name of a file to write the results when
# running trueTree. The results printed are various pieces of timing data, as well as the cophenetic
# correlation analysis. If the file specified already exists the results of the current run will be appended
# to the end of the file. If left to the default of "" all output will be displayed in the console.
	
#trueTree returns a trueTree object that contains the following attributes

# hclust - The hclust object formed from the original dataset

# edges - Table containing AU and BP values for every clade in the order that the clades were
# formed, as well as the standard error for both, the v and c found during the calculations.
# count - Table containing the number of times each clade in the original clustering was formed
# during the bootstrap for each r value used.

# msfit - msfit object created containing the results of the curve fitting used to calculate the au
# values. See pvclusts documentation for more information (Suzuki and Shimodaira, 2011)

# nboot The number of resampled databases generated for each r value

# r Vector containing all the r values used in the multiscale resampling

# storehClust - List containing all the hclust objects generating while running the bootstrap. This
# will be NULL if store was set to false in trueTree's parameters. The list will contain one sublist
# for each r value used, and each sublist has all the hclust objects generated for that r value.
# distance The distance matrix generated from the original dataset

# seed - A vector containing all the seeds for random number generator used to run the bootstrap.
# There will be one seed for each processor used.

# storeChunks - List containing all the resampled datasets objects generating while running the
# bootstrap. This will be NULL if storeChunks was set to false in trueTree's parameters. The list
# will contain one sublist for each r value used, and each sublist has all the resampled datasets
# generated for that r value.	

# storeCop - List containing all the cophenetic correlation coefficents generating while running the
# bootstrap. This will be NULL if storeChunks was set to false in trueTree's parameters. The list
# will contain one sublist for each r value used, and each sublist has all the cophenetic correlation coefficents
# generated for that r value.	

# metaTable - table containing only the metadata from the input file assuming said input file contained any metadata
	
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
	
	#write(copValues)
	copValues <- unlist(copValues)
	copValues <- sort(copValues) #sort the list
	copSize <- length(copValues)
	
	lowerBound = (1-confidenceInterval) / 2 #to get the lower bound subtract the confidence interval from 1 to get the precentage outside the interval and divide by 2 to get the precentage below the interval
	upperBound = 1 - lowerBound #subtract the lowerbound from 100% to get upper end of ranger
	
	write(paste("Original Cophenetic correlation", originalCor, sep=" "), file = logFileName, append=TRUE)
	write(paste("Number of Cophenetic correlation values", length(copValues), sep=" "), file = logFileName, append=TRUE)
	write(paste("Minimum Cophenetic correlation", min(copValues), sep=" "), file = logFileName, append=TRUE)
    write(paste(lowerBound * 100, "% interval Cophenetic correlation ", copValues[round(copSize * lowerBound)], sep=""), file = logFileName, append=TRUE)	
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
			if(.Platform$OS.type == "unix")
		  	{
				plot(pCluster, filename=outputFilename, main = main, height=height, width=width, specialLabels=specialLabels)
		 	}

			else
			{
				plot(pCluster, filename=outputFilename, main = main, height=height, width=width, specialLabels=specialLabels)
			}
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
		
		#create the histrogram
		hist(copValues)
		#abline(v = originalCor, col = "red") #add a line showing where original falls
		
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
	
	return(pCluster)
}

varianceTest <- function(input.file, distMetric = "euclidean" , clustMethod = "average" , input.transposed = TRUE, nboot = 10, runParallel = FALSE,
					numCPUs = 2, clusterType = 'SOCK', cladeChunkIn=NULL, rowSample=FALSE, r=seq(.5,1.4,by=.1), 
					metadata = FALSE, testRuns = 5)
{
# Function for testing how great the AU and BP values can vary for a particular dataset when using
# a specific set of parameters. varianceTest runs trueTree a user specified number of times, and finds
# the largest amount an AU value differs over all those runs and the largest amount a BP differs over
# all those runs.

# testRuns How many times should trueTree be run.

# ... All other parameters modify the behavior of trueTree when it is run and function identically
# to the parameters with the same name for the trueTree function.

		#do first test run
		result <- trueTree(input.file, distMetric = distMetric, clustMethod = clustMethod, input.transposed = input.transposed, nboot = nboot, runParallel = runParallel,
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
			result <- trueTree(input.file, distMetric = distMetric, clustMethod = clustMethod, input.transposed = input.transposed, nboot = nboot, runParallel = runParallel,
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
logFile <- ""
result <- trueTree("fedpapers.tsv", outputFilename = "dummy", nboot=1, distMetric = "euclidean", runParallel = FALSE, input.transposed = TRUE, numCPUs = 2, clusterType = "SOCK", plotOut=TRUE, logFileName = logFile, metadata=FALSE, width = 100000, height = 2000)

