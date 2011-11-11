source ( 'pvclust.R' )
source ( 'pvclust-internal.R' )

#library(pvclust)

# Idea for text labeling input
# textlabs is a list of character strings for the text part of the the label
# chunksize is a list of numeric for the number of chunks in each text
# textlabs and chunksize must have same length and are assumed to correspond elementwise
# so that first textlabs has first chunksize number of chunks

myCluster <- function(input.file , textlabs = NULL , chunksize = NULL ,
					distMetric = "euclidean" , clustMethod = "average" , main = "",
					input.transposed = TRUE, nboot = 100, runParallel = FALSE, clusterNumber = 2, clusterType = 'SOCK')
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
	#being able to run the code. The  clusterNumber, and clusterType parameters will also be ignored. If runParallel is set to true, the snow and snowfall librarys must be installed and the clusterNumber,
    #and clusterType parameters will be used
	
	#clusterType is the type of cluster to create. If set to 'SOCK' the cluster will use raw sockets, which due to security concerns will not work over a network with some computers. If set to MPI will use the MPI protocal
	#which requries that the computer has a MPI implementation, installed, properly set up and running, as well as having the r package rmpi installed. PVM (Parallel Virtual Machine) and NWS(networkspaces) 
	#are theoretically supported but not tested. See the snow/snowfall documentation for more information about the necessary steps to use those.
	
	#clusterNumber is the number of processers to use in the cluster. If the clusterType is set to SOCK all the processers used will be from the machine running this code. If clusterType is set to MPI the processers used 
	#can come from any machine mpi is set up to access. If clusterNumber is set to a value greater then the number of processers actually available the code will still run, but some processers will end up running multiple instances
	#of the program which will increase the total run time.
	
	
	library(stats)
	#change this for the text you'd like to input
	input.data <- read.table(as.character(input.file), header=T,
		comment.char="", row.names=1, sep="\t", quote="")

	#tTable <- ifelse( input.transposed, input.data, t( input.data ) )
	if ( input.transposed )
		tTable <- input.data
	else 
		tTable <- t( input.data )
	
	rowSums <- apply(tTable, 1, sum)
	denoms <- matrix(rep(rowSums, dim(tTable)[2]), byrow=F, ncol=dim(tTable)[2])
	relFreq <- tTable/denoms
	
       
	#find cophenetic corelation of orignal hclustering
    distTable <- dist(relFreq, method = distMetric)
	orignalClust <- hclust(distTable, method=clustMethod) #this can be taken out. get if from pvclust results     	
	orginalCopDist <- cophenetic(orignalClust)
	originalCor <- cor(orginalCopDist, distTable)

    #transpose relFreq so it matches the format pvclust expects
	relFreq <- t(relFreq)
	
	if( !is.null(textlabs) && !is.null(chunksize)) {
		if(length(textlabs) != length(chunksize)) stop("number of texts and corresponding chunk numbers must match")
		else {# check that sum(chunksize) == dim(relFreq)[1] , total number of chunks equals number of rows in relFreq
				L <- length(chunksize)
				temp <- NULL
				for(i in 1:L) {
					for(k in 1:chunksize[i]){
						temp <- c(temp,paste(textlabs[i],as.character(k),sep=""))
			}
				}
		row.names(relFreq) <- temp
	}
	}
	# else 0

	Sys.time()->start;
	
	copValues <- numeric(0)

	if(!runParallel)
	{
		pCluster <- pvclust(relFreq, nboot=nboot, method.hclust=clustMethod, method.dist=distMetric, storeCop=TRUE)
	}

	else
	{
		#ADD ERROR CHECKING
		library(snowfall) #I will write more about next line missing ip argument
		sfInit(parallel=TRUE, cpus=clusterNumber,type=clusterType) #This creates the cluster. If clusterType is SOCK all the processors will be taken from the current computer. If the type is MPI, the processors can be taken
		                                                           #from any computer running an MPI implementation
		cl <- sfGetCluster()
		sfSource( 'pvclust.R' )
		sfSource( 'pvclust-internal.R' )
		## parallel version of pvclust
		pCluster <- parPvclust(cl,relFreq, nboot=nboot, method.hclust=clustMethod, method.dist=distMetric, storeCop=TRUE, normalize=TRUE)
	}
	

	listH <- pCluster$storeCop
	for(i in listH)
	{
		for(j in i)
		{	
			copValues <- c(copValues, j)
		}
	}

	#print(class(pCluster))
	copValues <- sort(copValues)
	copSize <- length(copValues)
	
	upperBound = .95
	lowerBound = .05
	
	copValues <- copValues[(copSize * lowerBound):(copSize * upperBound)]
	
	print(copValues)
	
	print(paste("Original Cophenetic correlation", originalCor), sep=" ")
	print(paste("Number of Cophenetic correlation values", length(copValues)), sep=" ")
	print(paste("Minimum Cophenetic correlation", min(copValues)), sep=" ")
    print(paste("Maximum Cophenetic correlation", max(copValues)), sep=" ")
    print(paste("Mean Cophenetic correlation", mean(copValues)), sep=" ")
    print(paste("Median Cophenetic correlation", median(copValues)), sep=" ")
    print(paste("Standard Deviation Cophenetic correlation", sd(copValues)), sep=" ")
     
	#find range between 2.5 % in and 97.5 % sorted make parameters use order/sort

	## plot dendrogram with p-values
	# dev.control()
	#pdf("Testout.pdf" , onefile = TRUE, width=7.25, height=10)
	plot(pCluster)
	#dev.off()

	#ask.bak <- par()$ask
	#par(ask=TRUE)

	## highlight clusters with high au p-values
	#pvrect(pCluster)

	## print the result of multiscale bootstrap resampling
	#print(pCluster, digits=3)

	## plot diagnostic for curve fitting
	#msplot(pCluster, edges=c(2,4,6,7)) #note if the numbers in edges are higher then the number of actual edges (which is the number of observations minus 1)
								       #this line will not wok.

	#par(ask=ask.bak)

	## Print clusters with high p-values
	#pCluster.pp <- pvpick(pCluster)
	#pCluster.pp
	
	print(Sys.time()-start);
	
	if(runParallel)
	{
		sfStop()
	}
}

myCluster("inputTest.tsv", nboot=10, distMetric = "euclidean", runParallel = FALSE, input.transposed = TRUE, clusterNumber = 3)
