source ( 'pvclust.R' )
source ( 'pvclust-internal.R' )

#library(pvclust)

# Idea for text labeling input
# textlabs is a list of character strings for the text part of the the label
# chunksize is a list of numeric for the number of chunks in each text
# textlabs and chunksize must have same length and are assumed to correspond elementwise
# so that first textlabs has first chunksize number of chunks

myCluster <- function(input.file, filename = NULL, textlabs = NULL , chunksize = NULL ,
					distMetric = "euclidean" , clustMethod = "average" , main = "",
					input.transposed = TRUE, nboot = 100, runParallel = FALSE, clusterNumber = 2, clusterType = 'SOCK', confidenceInterval = .95, seed = NULL, cladeNumber=0, store=FALSE, storeChunks=FALSE, rowSample=FALSE)
{
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
	#being able to run the code. The clusterNumber, and clusterType parameters will also be ignored. If runParallel is set to true, the snow and snowfall librarys must be installed and the clusterNumber,
    #and clusterType parameters will be used. Currently there is one other discrepency that needs to be fixed. When runParallel is false each run will produce a new set of results, and when runParallel is true, each run will
	#produce the same set of results due to parpvclust giving each processor a fixed seed.
	
	#clusterType is the type of cluster to create. If set to 'SOCK' the cluster will use raw sockets, which due to security concerns will not work over a network with some computers. If set to MPI will use the MPI protocal
	#which requries that the computer has a MPI implementation, installed, properly set up and running, as well as having the r package rmpi installed. PVM (Parallel Virtual Machine) and NWS(networkspaces) 
	#are theoretically supported but not tested. See the snow/snowfall documentation for more information about the necessary steps to use those.
	
	#clusterNumber is the number of processers to use in the cluster. If the clusterType is set to SOCK all the processers used will be from the machine running this code. If clusterType is set to MPI the processers used 
	#can come from any machine mpi is set up to access. If clusterNumber is set to a value greater then the number of processers actually available the code will still run, but some processers will end up running multiple instances
	#of the program which will increase the total run time.
	
	#confidence interval is the confidence interval for the cophenetic correlations given as a percent entered in decimal form. The program will give the output the cophenetic correlations at the lower and upper ends of the interval. 
	#The lower bound can be computed as (100% - confidence interval) / 2 and the upper bound as (100% - lower bound). For example if the confidence interval is 95% the lowerbound will be 2.5% and the upper bound 97.5%

	#the seed parameter is used to set the seed used for random number generation. By default seed is set to null which will cause the program to use a random seed. If a value for seed is given the program will instead use that
	#value as the seed for random number generation. When not running in parallel a single value should be entered for seed. When running in parallel a list the same length as clusterNumber should be entered which each item in the list
	#being the seed for one particular processor. If seed is set to null and you want to check what the number generated for use as a seed was, it is stored as the attribute seed of the pvclust object returned
	#For example if you named the pvclust object returned result the new chunks can be accessed as follows: result$seed.
	
	#cladeNumber is a parameter to allow resampling by clades. Normally the resampling happens on a per chunk basis, meaning that the words in the new chunks created from the bootstraping are taken from the list of words in that chunk
	#when resampling is done by clades each chunk inside of one of the clades used in the resampling draws from all the words in the clade when resampling is perfomed.  The cladeNumber specifies how many clades are used. The dendrogram
	#will be broken up into cladeNumber clades using r's cutree function, and these clades are then the basis for the resampling
	
	#store is a boolean that sets pvclusts store parameter. If store is set to true then the result objected returned by the function will contain a store attribute which is a list holding all the hclust objects created during
	#the bootstrapping. The hclust obects will be stored in a list as the attribute store of the of the pvclust object returned
	#For example if you named the pvclust object returned result the new chunks can be accessed as follows: result$store. The format of store is as a list containing a number of sublists with there being one sublist for each
	#r value pvclust used, and each sublist having all the hclust objects created for that r value
	#WARNING: setting store to true will likely cause R to run out of memory if nboot is set to a high value like 100,000. As such this parameter should be set to false when preforming actual validation
	#and only used on smaller test runs to see what kind of dendrograms are being created. 
	
	#storeChunks is a boolean that sets pvclusts store parameter. If store is set to true then the result objected returned by the function will contain a store attribute which is a list holding all the new chunks created during
	#the bootstrapping. The new chunks	will be stored in a list as the attribute storeChunk of the pvclust object returned
	#For example if you named the pvclust object returned result the new chunks can be accessed as follows: result$storeChunks. The format of storeChunks is as a list containing
	#a number of sublists with there being one sublist for each r value pvclust used, and each sublist having all the new chunks created for that r value
	#WARNING: setting store to true will likely cause R to run out of memory if nboot is set to a high value like 100,000. As such this parameter should be set to false when preforming actual validation
	#and only used on smaller test runs to see what kind of dendrograms are being created. 
	
	#rowSample is a boolean that if true allows the use of pvclusts original resampling method, where instead of basing the resampling around the number of times each individual word appeared in each individual clade, it builds
	#the new data table by picking out a new set of words to use for all the chunks keeping the counts for the number of times those words appeared in the orignal text for those chunks. If storeRow is set to true then the 
	#cladeNumber parameter will be ignored. WARNING. It is currently unclear whether this is actually a good method of resampling. By default this parameter is set to false and for now it is advised that it be left to false
	#unless you have a specific reason for wanting to use it.
	
	#This function will return a list containing a modified pvclust object with a number of additional attributes. The storeCop attribute is a list of lists with each sublist containing the cophenetic correlation values
	#for all the hclust objects generated for a particular r value. The first sublist contains the values for the first r value used and so on. There is also a seed attribute and a storeChunks attributes which are described
	#in the documentation for the parameters with the same names. If you want to see what the r values where they are stored as a list in the attribute r
	
	library(stats)
	Sys.time()->startTotal; #holds start time of program so time to run entire program can be calculated
	Sys.time()->startSection; #start timing the bootstraping
	
	#set inital values for certain variables
	storeData <- NULL #st
	
	#change this for the text you'd like to input
	input.data <- read.table(as.character(input.file), header=T,
		comment.char="", row.names=1, sep="\t", quote="")

	#tTable <- ifelse( input.transposed, input.data, t( input.data ) )
	if ( input.transposed ) #if the input
		tTable <- input.data
	else
		tTable <- t( input.data )

    #transpose data so it matches format pvclust expects
    tTable <- t(tTable)

	#Old code to change the labels of the chunks. Hasn't been used in a while and probally doesn't work anymore with all the changes since it was written. Still leaving it in
	#as a starting point in case someone decides it should be reimplemented.
	#if( !is.null(textlabs) && !is.null(chunksize)) {
	#	if(length(textlabs) != length(chunksize)) stop("number of texts and corresponding chunk numbers must match")
	#	else {# check that sum(chunksize) == dim(relFreq)[1] , total number of chunks equals number of rows in relFreq
	#			L <- length(chunksize)
	#			temp <- NULL
	#			for(i in 1:L) {
	#				for(k in 1:chunksize[i]){
	#					temp <- c(temp,paste(textlabs[i],as.character(k),sep=""))
	#		}
	#			}
	#	row.names(relFreq) <- temp
	#}
	#}
	# else 0

	copValues <- numeric(0)

	print("Inital runtime:")
    print(Sys.time()-startSection);
	Sys.time()->startSection; #start timing the bootstraping
	
	if(!runParallel)
	{
		pCluster <- pvclust(tTable, nboot=nboot, method.hclust=clustMethod, method.dist=distMetric, storeCop=TRUE, seed=seed, cladeNumber=cladeNumber, store=store, storeChunks=storeChunks, rowSample=rowSample)
	}

	else
	{
		#ADD MORE ERROR CHECKING
		library(snowfall) #I will write more about next line missing ip argument
		sfInit(parallel=TRUE, cpus=clusterNumber,type=clusterType) #This creates the cluster. If clusterType is SOCK all the processors will be taken from the current computer. If the type is MPI, the processors can be taken
		                                                           #from any computer running an MPI implementation, and the MPI program will handle choosing what processors to use
																   
		if(!sfIsRunning())
		{
			print("Error. Cluster not created successfully. Will use nonparallel version of pvclust instead")
			pCluster <- pvclust(tTable, nboot=nboot, method.hclust=clustMethod, method.dist=distMetric, storeCop=TRUE, seed=seed, cladeNumber=cladeNumber, store=store, storeChunks=storeChunks, rowSample=rowSample)			
		}
		cl <- sfGetCluster()
		sfSource( 'pvclust.R' )
		sfSource( 'pvclust-internal.R' )
		## parallel version of pvclust
		pCluster <- parPvclust(cl,tTable, nboot=nboot, method.hclust=clustMethod, method.dist=distMetric, storeCop=TRUE, normalize=TRUE, seed=seed, cladeNumber=cladeNumber, store=store, storeChunks=storeChunks, rowSample=rowSample)
	}
	
	print("Bootstrap runtime:")
    print(Sys.time()-startSection);
	Sys.time()->startSection; #start timing the analysis of the cophenetic correlations

	#find cophenetic correlation of orignal hclustering	
	orginalCopDist <- cophenetic(pCluster$hclust) #get cophenetic distance matrix for original hclust object
	originalCor <- cor(orginalCopDist, pCluster$distance) #get correlation with original distance matrix

	#the cophenentic correlations are stored in a series of lists which are themselves stored in a list, so merge them all into one list
	listH <- pCluster$storeCop
	for(i in listH)
	{
		for(j in i)
		{	
			copValues <- c(copValues, j)
		}
	}

	copValues <- sort(copValues) #sort the list
	copSize <- length(copValues)
	
	lowerBound = (1-confidenceInterval) / 2 #to get the lower bound subtract the confidence interval from 1 to get the precentage outside the interval and divide by 2 to get the precentage below the interval
	upperBound = 1 - lowerBound #subtract the lowerbound from 100% to get upper end of ranger
	
	print(paste("Original Cophenetic correlation", originalCor), sep=" ")
	print(paste("Number of Cophenetic correlation values", length(copValues)), sep=" ")
	print(paste("Minimum Cophenetic correlation", min(copValues)), sep=" ")
    print(paste(lowerBound * 100, "% interval Copheneticcorrelation ", copValues[round(copSize * lowerBound)]), sep="")	
    print(paste("Median Cophenetic correlation", median(copValues)), sep=" ")
	print(paste(upperBound * 100, "% interval Cophenetic correlation ", copValues[round(copSize * upperBound)]), sep="")	
    print(paste("Maximum Cophenetic correlation", max(copValues)), sep=" ")
	
	print(paste("Mean Cophenetic correlation", mean(copValues)), sep=" ")
    print(paste("Standard Deviation Cophenetic correlation", sd(copValues)), sep=" ")
    
	print("Cophenetic analysis runtime:")
    print(Sys.time()-startSection);
	
	#find range between 2.5 % in and 97.5 % sorted make parameters use order/sort
	
	if(!is.null(filename))
	{
		plot(pCluster, filename)
		dev.off()
	}
	
	else
	{
		plot(pCluster)
	}
	
	print("Total runtime:")
    print(Sys.time()-startTotal);
	
	if(runParallel)
	{
		sfStop() #end the cluster
	}
	
	return(pCluster)
}

#result <- myCluster("inputTestTwoClose.tsv", filename = "inputTestTwoClose.png", nboot=100000, distMetric = "euclidean", runParallel = TRUE, input.transposed = TRUE, clusterNumber = 2, clusterType = "SOCK")
#print("Done 1")
#result <- myCluster("inputTestTwoClose2.tsv", filename = "inputTestTwoClose2.png", nboot=100000, distMetric = "euclidean", runParallel = TRUE, input.transposed = TRUE, clusterNumber = 2, clusterType = "SOCK")
#print("Done 2")
#result <- myCluster("inputTestTwoClose3.tsv", filename = "inputTestTwoClose3.png", nboot=100000, distMetric = "euclidean", runParallel = TRUE, input.transposed = TRUE, clusterNumber = 2, clusterType = "SOCK")
#print("Done 3")
#result <- myCluster("inputTestTwoClose4.tsv", filename = "inputTestTwoClose4.png", nboot=100000, distMetric = "euclidean", runParallel = TRUE, input.transposed = TRUE, clusterNumber = 2, clusterType = "SOCK")
#print("Done 4")
#result <- myCluster("inputTestTwoClose.tsv", filename = "rowinputTestTwoClose.png", nboot=100000, distMetric = "euclidean", runParallel = TRUE, input.transposed = TRUE, clusterNumber = 2, clusterType = "SOCK", rowSample = TRUE)
#print("Done 5")
#result <- myCluster("inputTestTwoClose2.tsv", filename = "rowinputTestTwoClose.png", nboot=100000, distMetric = "euclidean", runParallel = TRUE, input.transposed = TRUE, clusterNumber = 2, clusterType = "SOCK", rowSample = TRUE)
#print("Done 6")
#result <- myCluster("inputTestTwoClose3.tsv", filename = "rowinputTestTwoClose.png", nboot=100000, distMetric = "euclidean", runParallel = TRUE, input.transposed = TRUE, clusterNumber = 2, clusterType = "SOCK", rowSample = TRUE)
#print("Done 7")
#result <- myCluster("inputTestTwoClose4.tsv", filename = "rowinputTestTwoClose.png", nboot=100000, distMetric = "euclidean", runParallel = TRUE, input.transposed = TRUE, clusterNumber = 2, clusterType = "SOCK", rowSample = TRUE)
#print("All Done")