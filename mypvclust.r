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
					input.transposed = TRUE, nboot = 100, runParallel = FALSE)
{
	## List of possible distance metrics
	## METHODS <- c("euclidean", "maximum", "manhattan", "canberra",
	## "binary", "minkowski", "correlation", "uncentered", "abscor")

	## List of possible cluster-distance methods
	## METHODS <- c("ward", "single", "complete", "average", "mcquitty",
	## "median", "centroid")
	
	## List if possible output.type for writing to file
	## OUTPUT.TYPE <- c( "pdf", "svg", "phyloxml" )

	## If want to dump to stdout, do not provide output.file

	## If the input has text names on top, set input.transposed <- FALSE

	##

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
	orignalClust <- hclust(distTable, method=clustMethod)     	
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
		library(snowfall)
		sfInit(parallel=TRUE, cpus=2,type='SOCK')
		cl <- sfGetCluster()
		#cl <- makeCluster(c("localhost", "localhost", "localhost"), type = "SOCK", homogeneous = TRUE)
		sfSource( 'pvclust.R' )
		sfSource( 'pvclust-internal.R' )
		#clusterEvalQ(cl, library(snow))
		#clusterEvalQ(cl, library(stats))
		## parallel version of pvclust
		pCluster <- parPvclust(cl,relFreq, nboot=nboot, method.hclust=clustMethod, method.dist=distMetric, storeCop=TRUE)
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
	print(paste("Original Cophenetic correlation", originalCor), sep=" ")
	print(paste("Minimum Cophenetic correlation", min(copValues)), sep=" ")
      print(paste("Maximum Cophenetic correlation", max(copValues)), sep=" ")
      #print(paste("Average Cophenetic correlation", average(copValues)), sep=" ")	

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

myCluster("danile-azarius.txt", nboot=100000, distMetric = "euclidean", runParallel = FALSE, input.transposed = FALSE)