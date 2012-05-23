currentNode = 1 #global variable to keep track of position while recursively travelling dendrogram during the plotting code

pvclust <- function(data, method.hclust="average",
                    method.dist="correlation", use.cor="pairwise.complete.obs",
                    nboot=1000, r=seq(.5,1.4,by=.1), weight=FALSE, normalize=TRUE, seed=NULL, cladeChunkIn=NULL, rowSample=FALSE, store=FALSE, storeCop=FALSE, storeChunks=FALSE)
  {
  
	if(is.null(seed)) #if no seed was specified use the system time as the seed which should effectively be random
	{
		seed = as.numeric(Sys.time())
	}

	set.seed(seed)
  
    copDistance <- NULL #set to null in case cophenetic correlations aren't beening looked for

    # data: (n,p) matrix, n-samples, p-variables
    n <- nrow(data); p <- ncol(data)

    #normalize data before getting distance matrix
    colSums <- apply(data, 2, sum) #each example/observation/object is one column, so find the sums of the columns
    denoms <- matrix(rep(colSums, dim(data)[1]), byrow=T, ncol=dim(data)[2]) #compute matrix to divide current matrix by to normalize matrix. Each entry in a column is the sum of the column
    relFreq <- data/denoms
	
	# hclust for original data
    METHODS <- c("ward", "single", "complete", "average", "mcquitty",
                 "median", "centroid")
    method.hclust <- METHODS[pmatch(method.hclust, METHODS)]
    distance <- dist.pvclust(relFreq, method=method.dist, use.cor=use.cor)

    data.hclust <- hclust(distance, method=method.hclust)
	
	if(is.null(cladeChunkIn)) #if resampling by chunk
	{
		cladeChunkIn <- 1:ncol(data) #holds which clade each chunk is in
	}
	
	else #there is a desired number of clades to resample by 
	{
		cladeCounts <- matrix(rep(0, max(cladeChunkIn * n)), ncol = max(cladeChunkIn), nrow = n)
		for(i in 1:p) #for each chunk
		{
			curClade <- cladeChunkIn[i] #get which clade the chunk is in
			cladeCounts[,curClade] <- cladeCounts[,curClade] + data[,i] #add counts to total	
		}
		
		#convert to relative frequency
		colSums <- apply(cladeCounts, 2, sum) #each clade is one column, so find the sums of the columns
		denoms <- matrix(rep(colSums, dim(cladeCounts)[1]), byrow=T, ncol=dim(cladeCounts)[2]) #compute matrix to divide current matrix by to normalize matrix. Each entry in a column is the sum of the column
		relFreq <- cladeCounts/denoms
	}
	
	chunkSize <- list() #stores total number of words in the chunk
	
	for(i in 1:ncol(data))
	{
		chunkSize[[i]] <- sum(data[,i]) #total number of words in the chunk is sum of the number of each individual word
	}

    #if finding the cophenetic correlations
    if(storeCop)
    {
        copDistance <- distance
    }

    # multiscale bootstrap
    size <- floor(n*r) #get the size of the samples for each stage of the multiscale bootstraping
    rl <- length(size) #get the number of stages

    if(rl == 1) {
      if(r != 1.0)
        warning("Relative sample size r is set to 1.0. AU p-values are not calculated\n")

      r <- list(1.0)
    }
    else
	{
      r <- as.list(size/n) #recalculate r so R matchs size/n exactly instead of approxiametly
	}
	
    mboot <- lapply(r, boot.hclust, data=data, object.hclust=data.hclust, nboot=nboot,
                    method.dist=method.dist, use.cor=use.cor,
                    method.hclust=method.hclust, store=store, weight=weight, storeCop=storeCop, copDistance=copDistance, normalize=normalize, cladeChunkIn=cladeChunkIn, chunkSize=chunkSize,
					storeChunks=storeChunks, rowSample=rowSample, relFreq = relFreq) #do the actual bootstraping

    result <- pvclust.merge(data=data, object.hclust=data.hclust, mboot=mboot, distance=distance, seed=seed)
    
    return(result) 
  }

#gives the node its proper color
lineColor <- function(x, colorOfNodes)
{
	attr(x, "nodePar") <- list("pch"  = NA, "lab.col"  = colorOfNodes[[currentNode]])
	attr(x, "edgePar") <- list("col"  = colorOfNodes[[currentNode]])
	assign("currentNode",  currentNode + 1, envir = .GlobalEnv)
	
	x #this line is necessary for some reason for the dendrapply function that calls this to work
}

#get the color for a given leaf based on it's label
getColor <- function(label, specialLabels, metaTable = NULL)
{
	if(label %in% specialLabels) #if label is one of the labels to watch out for
	{
		return("gold")
	}

	else if(is.null(metaTable)) #if there is no metadata set the node to black. 
	{
		return("black")
	}
	
	else if(metaTable[label, 3] == "Bacteria")
	{
		return("green") 
	}
	
	else #if the node is not bacteria it must be archaea
	{
		return("red")
	}
}

#generates an list containing the color for every node in the tree
#rules are a node is red if it only contains Archeia, green if it only contains Bacteria, Gold if it was specially selected for highlighting, and blue if it contains mulitple of the previous categories'
#The list is ordered in the order that nodes are visited by dendrapply.
generateLineColorList <- function(x, mergeTableRow, specialLabels, metaTable = NULL)
{
	colorlist <- list()
	
	#color the left half of the clade
	if(x$merge[mergeTableRow,1] < 0) #if the left node is a chunk determine the chunk's color
	{
		leftColor <- getColor(x$labels[-x$merge[mergeTableRow,1]], specialLabels=specialLabels, metaTable = metaTable) #the color of the chunk
		leftList <- list(leftColor) #list of the colors of all the nodes to the left
	}
	
	else #if the left node is a clade recursively run the function on that clade
	{
		result <- generateLineColorList(x, x$merge[mergeTableRow,1], specialLabels=specialLabels, metaTable = metaTable) 
		leftColor <- result$color #the overall color of the subclade
		leftList <- result$colorList #list of the colors of all the nodes to the left
	}
	
	#color the right half of the clade
	if(x$merge[mergeTableRow,2] < 0) #if the right node is a chunk determine the chunk's color
	{
		rightColor <- getColor(x$labels[-x$merge[mergeTableRow,2]], specialLabels=specialLabels, metaTable = metaTable) #the color of the chunk
		rightList <- list(rightColor) #list of the colors of all the nodes to the right
	}
	
	else #if the right node is a clade recursively run the function on that clade
	{
		result <- generateLineColorList(x, x$merge[mergeTableRow,2], specialLabels=specialLabels, metaTable = metaTable) 
		rightColor <- result$color #the overall color of the subclade
		rightList <- result$colorList #list of the colors of all the nodes to the right
	}
	
	if(leftColor == rightColor) #check if the colors of the two subclades of the current clade are the same 
	{
		color <- leftColor #if so use the color they share
	}
	
	else #if the colors are different the subclades have different contents
	{
		color <- "blue" #set the clade to blue to mark it's mixed contents
	}

	#the colors found need to be put together in the proper order. The current clade has one node for each of it's childern which contains a clade instead of just a chunk. 
	#Those nodes need to be given the color of the current clade, but only if they exist. These nodes will appear in the list of colors before all the colors for the nodes in the respective
	#subclades
	
	if(x$merge[mergeTableRow,1] > 0  && x$merge[mergeTableRow,2] > 0) #if both childern are subclades
	{
		colorList <- c(color, leftList, color, rightList) #both nodes in the current clade exist so add them into the color list
	}
	
	else if(x$merge[mergeTableRow,1] > 0) #if the right child is a chunk
	{
		colorList <- c(color, leftList, rightList) #there is only a node for the left clade so add that to the color list
	}
	
	else if(x$merge[mergeTableRow,2] > 0) #if the left child is a chunk
	{
		colorList <- c(leftList, color, rightList) #there is only a node for the right clade so add that to the color list
	}
	
	else #both children are individual chunks
	{
	colorList <- append(leftColor, rightColor)
	}
	
	result <- list(colorList=colorList, color=color)
	return(result)
}
  
 #plots a pvclust object
plot.trueTree <- function(x, outputFilename = NULL, print.pv=TRUE, print.num=TRUE, float=0.01,
                         col.pv=c(2,3,8), cex.pv=0.8, font.pv=NULL,
                         col=NULL, cex=NULL, font=NULL, lty=NULL, lwd=NULL,
                         main=NULL, sub=NULL, xlab=NULL, height=800, width=800, specialLabels=NULL, showBP=FALSE, ...)
{
  if(.Platform$OS.type == "windows")
  {
  	if(!is.null(outputFilename))
  	{
		png(paste(outputFilename, ".png", sep=""), width=width, height=height)
  	}
	
	else
	{
		windows(width=width, height=height)
	}
  }
  
  else if(.Platform$OS.type == "unix")
  {
	if(!is.null(outputFilename))
 	{
		if(width > 32766)
		{
				png(paste(outputFilename, ".png", sep=""), width=32766, height=height, type="Xlib")
		}
		
		else
		{
				png(paste(outputFilename, ".png", sep=""), width=width, height=height, type="Xlib")
		}
	}

	else
	{
		X11(width=floor(width/96), height=floor(height/96), type="Xlib") #X11 specifies window size in inches for some bizare reason
									 #I'm not certain of the correct conversion factor but this is my best
									 #guess 
	}
  }
  
  metaTable <- x$metaTable[[1]] #get metadata out of pvclust object
  
  #The line describing the dendrogram needs to be changed depending on if bp values are being displayed or not
  
  if(!showBP) #if not showing bp values
  {
	  if(is.null(main))
		main <- paste("Cluster dendrogram with AU values (%)", paste("Cluster method: ", x$hclust$method, sep=""), paste("Distance: ", x$hclust$dist.method), sep = "\n")

	  else
		main <- paste(main, "Cluster dendrogram with AU values (%)", paste("Cluster method: ", x$hclust$method, sep=""), paste("Distance: ", x$hclust$dist.method), sep = "\n")
  }
  
  else #if not showing bp values
  {
	  if(is.null(main))
		main <- paste("Cluster dendrogram with AU/BP values (%)", paste("Cluster method: ", x$hclust$method, sep=""), paste("Distance: ", x$hclust$dist.method), sep = "\n")

	  else
		main <- paste(main, "Cluster dendrogram with AU/BP values (%)", paste("Cluster method: ", x$hclust$method, sep=""), paste("Distance: ", x$hclust$dist.method), sep = "\n")
  }

  if(is.null(sub))
    #sub=paste("Cluster method: ", x$hclust$method, sep="")
 
  if(is.null(xlab))
    #xlab=paste("Distance: ", x$hclust$dist.method)  
	
  dend <- as.dendrogram(x$hclust) #convert the hclust object into a dendrogram object
  colorList <- generateLineColorList(x$hclust, dim(x$hclust$merge)[1], specialLabels=specialLabels, metaTable = metaTable) #figure out what color each node should be
  colorList <- c(0, colorList$colorList) #the first node checked be dendrapply doesn't seem to be part of the dendrogram so add a dummy value at the start of the list
  
  assign("currentNode",  1, envir = .GlobalEnv) #currentNode is a global variable to keep track of where in the tree we are
  dend <- dendrapply(dend, lineColor, colorList) #add color to all the nodes in the tree

  #find length of longest chunk name
  maxL <- max( nchar( x$hclust$labels ))
  
  # set margins so there is just enough room for the labels
  # The numbers measure margin size in line units
  # The paramets are the size of the bottom,left,top,right margins
  # On average a margin one line wide seems to have room for about 2.5 characters)
  # so the margin on the bottom is set to the number of lines necessary to display
  # the longest label if there was only 2 characters per line which leave's a decent buffer
  par( mar=c((maxL / 2.0), 2.1, 4.1, 2.1))
	
  plot(dend, main=main, sub=sub, xlab="", col=col, cex=cex,
       font=font, lty=lty, lwd=lwd, ...)
  if(print.pv)
    text(x, col=col.pv, cex=cex.pv, font=font.pv, float=float, print.num=print.num, showBP = showBP)

  if(!is.null(outputFilename)) #if writing to a file close the connection
  {
	dev.off()
  }
}

#this function handles the actual writing of the au and bp labels on the plot
text.trueTree <- function(x, col=c(2,3,8), print.num=TRUE,  float=0.01, cex=NULL, font=NULL, showBP = FALSE, ...)
{
  axes <- hc2axes(x$hclust)
  usr  <- par()$usr; wid <- usr[4] - usr[3]
  au <- as.character(round(x$edges[,"au"]*100))
  bp <- as.character(round(x$edges[,"bp"]*100))
  rn <- as.character(row.names(x$edges))
  au[length(au)] <- "au"
  bp[length(bp)] <- "bp"
  rn[length(rn)] <- "edge #"
  a <- text(x=axes[,1], y=axes[,2] + float * wid, au,
            col=col[1], pos=2, offset=.3, cex=cex, font=font)
  if(showBP)
  {
  a <- text(x=axes[,1], y=axes[,2] + float * wid, bp,
            col=col[2], pos=4, offset=.3, cex=cex, font=font)
  }
  if(print.num)
  {
    a <- text(x=axes[,1], y=axes[,2], rn,
              col=col[3], pos=1, offset=.3, cex=cex, font=font)
  }
}

print.trueTree <- function(x, which=NULL, digits=3, ...)
{
  if(is.null(which)) which <- 1:nrow(x$edges)
  cat("\n")
  cat(paste("Cluster method: ", x$hclust$method, "\n", sep=""))
  cat(paste("Distance      : ", x$hclust$dist.method, "\n\n", sep=""))
  cat("Estimates on edges:\n\n")
  print(round(x$edges[which,], digits=digits))
  cat("\n")
}

summary.trueTree <- function(object, ...){
  class(object) <- "list"
  summary(object, ...)
}

pvrect <- function(x, alpha=0.95, pv="au", type="geq", max.only=TRUE, border=2, ...)
  {
    len <- nrow(x$edges)
    member <- hc2split(x$hclust)$member
    order  <- x$hclust$order
    usr <- par("usr")
    xwd <- usr[2] - usr[1]
    ywd <- usr[4] - usr[3]
    cin <- par()$cin

    ht <- c()
    j <- 1

    if(is.na(pm <- pmatch(type, c("geq", "leq", "gt", "lt"))))
       stop("Invalid type argument: see help(pvrect)")
    
    for(i in (len - 1):1)
      {
        if     (pm==1) wh <- (x$edges[i,pv] >= alpha) # Greater than or EQuals
        else if(pm==2) wh <- (x$edges[i,pv] <= alpha) # Lower than or EQuals
        else if(pm==3) wh <- (x$edges[i,pv] >  alpha) # Greater Than
        else if(pm==4) wh <- (x$edges[i,pv] >  alpha) # Lower Than

        if(wh)
          {
            mi <- member[[i]]
            ma <- match(mi, order)
            
            if(max.only == FALSE || (max.only && sum(match(ma, ht, nomatch=0)) == 0))
              {
                xl <- min(ma)
                xr <- max(ma)
                yt <- x$hclust$height[i]
                yb <- usr[3]
                
                mx <- xwd / length(member) / 3
                my <- ywd / 200
                
                rect(xl - mx, yb + my, xr + mx, yt + my, border=border, shade=NULL, ...)
                
                j <- j + 1
              }
            ht <- c(ht, ma)
          }
      }
  }

msplot <- function(x, edges=NULL, ...) 
  {
    if(is.null(edges)) edges <- 1:length(x$msfit)
    d   <- length(edges)

    mfrow.bak <- par()$mfrow
    on.exit(par(mfrow=mfrow.bak))

    par(mfrow=n2mfrow(d))

    for(i in edges) {
      if(i == 1 || (i %% 10 == 1 && i > 20))
        main <- paste(i, "st edge", sep="")
      else if(i == 2 || (i %% 10 == 2 && i > 20))
        main <- paste(i, "nd edge", sep="")
      else if(i == 3 || (i %% 10 == 3 && i > 20))
        main <- paste(i, "rd edge", sep="")
      else
        main <- paste(i, "th edge", sep="")

      plot(x$msfit[[i]], main=main, ...)
    }
  }

lines.trueTree <- function(x, alpha=0.95, pv="au", type="geq", col=2, lwd=2, ...)
  {
    len <- nrow(x$edges)
    member <- hc2split(x$hclust)$member
    order  <- x$hclust$order
    usr <- par("usr")
    xwd <- usr[2] - usr[1]
    ywd <- usr[4] - usr[3]
    cin <- par()$cin

    ht <- c()
    j <- 1
    
    if(is.na(pm <- pmatch(type, c("geq", "leq", "gt", "lt"))))
       stop("Invalid type argument: see help(lines.pvclust)")
    
    for(i in (len - 1):1)
      {
        if     (pm==1) wh <- (x$edges[i,pv] >= alpha) # Greater than or EQuals
        else if(pm==2) wh <- (x$edges[i,pv] <= alpha) # Lower than or EQuals
        else if(pm==3) wh <- (x$edges[i,pv] >  alpha) # Greater Than
        else if(pm==4) wh <- (x$edges[i,pv] >  alpha) # Lower Than

        if(wh)
          {
            mi <- member[[i]]
            ma <- match(mi, order)
            
            if(sum(match(ma, ht, nomatch=0)) == 0)
              {
                xl <- min(ma)
                xr <- max(ma)
                yt <- x$hclust$height[i]
                yb <- usr[3]
                
                mx <- xwd/length(member)/10
                
                segments(xl-mx, yb, xr+mx, yb, xpd=TRUE, col=col, lwd=lwd, ...)

                j <- j + 1
              }
            ht <- c(ht, ma)
          }
      }
  }

pvpick <- function(x, alpha=0.95, pv="au", type="geq", max.only=TRUE)
  {
    len <- nrow(x$edges)
    member <- hc2split(x$hclust)$member
    order  <- x$hclust$order
    
    ht <- c()
    a  <- list(clusters=list(), edges=c()); j <- 1

    if(is.na(pm <- pmatch(type, c("geq", "leq", "gt", "lt"))))
       stop("Invalid type argument: see help(pickup)")
    
    for(i in (len - 1):1)
      {
        if     (pm==1) wh <- (x$edges[i,pv] >= alpha) # Greater than or Equals
        else if(pm==2) wh <- (x$edges[i,pv] <= alpha) # Lower than or Equals
        else if(pm==3) wh <- (x$edges[i,pv] >  alpha) # Greater Than
        else if(pm==4) wh <- (x$edges[i,pv] >  alpha) # Lower Than

        if(wh)
          {
            mi <- member[[i]]
            ma <- match(mi, order)

            if(max.only == FALSE || (max.only && sum(match(ma, ht, nomatch=0)) == 0))
              {
                a$clusters[[j]] <- x$hclust$labels[mi]
                a$edges <- c(a$edges,i)
                
                j <- j + 1
              }
            ht <- c(ht, ma)
          }
      }
    
    a$edges <- a$edges[length(a$edges):1]
    a$clusters <- a$clusters[length(a$edges):1]

    return(a)
  }

parPvclust <- function(cl, data, method.hclust="average",
                       method.dist="correlation", use.cor="pairwise.complete.obs",
                       nboot=1000, r=seq(.5,1.4,by=.1),
                       weight=FALSE, normalize=TRUE,
                       init.rand=TRUE, seed=NULL, cladeChunkIn=NULL, rowSample=FALSE, store=FALSE, storeCop=FALSE, storeChunks=FALSE)
  {

    if(!(require(snow))) stop("Package snow is required for parPvclust.")

    if((ncl <- length(cl)) < 2 || ncl > nboot) { #if nboot is less then the number of clusters
      warning("Too small value for nboot: non-parallel version is executed.")
	  return(pvclust(data,method.hclust,method.dist,use.cor,nboot,r,weight,normalize,NULL,cladeChunkIn,rowSample,store,storeCop,storeChunks))
    }

    copDistance <- NULL #set to null in case cophenetic correlations aren't beening looked for
    
	#normalize data before getting distance matrix
    colSums <- apply(data, 2, sum) #each example/observation/object is one column, so find the sums of the columns
    denoms <- matrix(rep(colSums, dim(data)[1]), byrow=T, ncol=dim(data)[2]) #compute matrix to divide current matrix by to normalize matrix. Each entry in a column is the sum of the column
    relFreq <- data/denoms
	
	if(init.rand) {
	#give all the processers a unique random seed
      if(is.null(seed)) #if the user didn't supply seeds
	  {
		curTime <- as.numeric(Sys.time())
        seed <- curTime:(curTime+length(cl) - 1) #start the seed at the current time and increment the seed by 1 for each additional processor 
      }
	  else if(length(seed) != length(cl)) #if the user supplied seeds make sure there is exactly one seed per processor
        stop("seed and cl should have the same length.")
      
      # setting random seeds
      parLapply(cl, as.list(seed), set.seed) #give each processor in the cluster a seed equal to it's processer number or equal to the seed the user supplied
    }

    # data: (n,p) matrix, n-samples, p-variables
    n <- nrow(data); p <- ncol(data)

    # hclust for original data
    METHODS <- c("ward", "single", "complete", "average", "mcquitty", 
                 "median", "centroid")
    method.hclust <- METHODS[pmatch(method.hclust, METHODS)]
    distance <- dist.pvclust(relFreq, method=method.dist, use.cor=use.cor)
    data.hclust <- hclust(distance, method=method.hclust)
	
	if(is.null(cladeChunkIn)) #if resampling by chunk
	{
		cladeChunkIn <- 1:ncol(data) #holds which clade each chunk is in
	}
	
	else #there is a desired number of clades to resample by 
	{
		cladeCounts <- matrix(rep(0, max(cladeChunkIn * n)), ncol = max(cladeChunkIn), nrow = n)
		for(i in 1:p) #for each chunk
		{
			curClade <- cladeChunkIn[i] #get which clade the chunk is in
			cladeCounts[,curClade] <- cladeCounts[,curClade] + data[,i] #add counts to total	
		}
		
		#convert to relative frequency
		colSums <- apply(cladeCounts, 2, sum) #each clade is one column, so find the sums of the columns
		denoms <- matrix(rep(colSums, dim(cladeCounts)[1]), byrow=T, ncol=dim(cladeCounts)[2]) #compute matrix to divide current matrix by to normalize matrix. Each entry in a column is the sum of the column
		relFreq <- cladeCounts/denoms
	}
	
	chunkSize <- list() #stores total number of words in the chunk
	
	for(i in 1:ncol(data))
	{
		chunkSize[[i]] <- sum(data[,i]) #total number of words in the chunk is sum of the number of each individual word
	}

	
    #if finding the cophenetic correlations
    if(storeCop)
    {
        copDistance <- distance
    }


    # multiscale bootstrap
    size <- floor(n*r)
    rl <- length(size)
    
    if(rl == 1) {
      if(r != 1.0)
        warning("Relative sample size r is set to 1.0. AU p-values are not calculated\n")
      
      r <- list(1.0)
    }
    else
      r <- as.list(size/n)

    nbl <- as.list(rep(nboot %/% ncl,times=ncl)) # %/% is integer division. Divide nboot up evenly across the processers in the cluster
    
    if((rem <- nboot %% ncl) > 0) #if there are some nboots remaining
    nbl[1:rem] <- lapply(nbl[1:rem], "+", 1) #add 1 nboot to each cluster upto the number of remaining bootstraps

    cat("Multiscale bootstrap... ")
    
    mlist <- parLapply(cl, nbl, pvclust.node,
                       r=r, data=data, object.hclust=data.hclust, method.dist=method.dist,
                       use.cor=use.cor, method.hclust=method.hclust,
                       store=store, weight=weight, storeCop=storeCop, copDistance=copDistance, normalize=normalize, cladeChunkIn=cladeChunkIn, chunkSize=chunkSize,
					   storeChunks=storeChunks, rowSample=rowSample, relFreq = relFreq) #do the bootstraping
    cat("Done.\n")
    
	
    mboot <- mlist[[1]]

    for(i in 2:ncl) { #merge all the data into a single object
      for(j in 1:rl) {
        mboot[[j]]$edges.cnt <- mboot[[j]]$edges.cnt + mlist[[i]][[j]]$edges.cnt
        mboot[[j]]$nboot <- mboot[[j]]$nboot + mlist[[i]][[j]]$nboot
        mboot[[j]]$store <- c(mboot[[j]]$store, mlist[[i]][[j]]$store)
        mboot[[j]]$storeCop <- c(mboot[[j]]$storeCop, mlist[[i]][[j]]$storeCop)
        mboot[[j]]$storeChunks <- c(mboot[[j]]$storeChunks, mlist[[i]][[j]]$storeChunks)
      }
    }

    result <- pvclust.merge( data=data, object.hclust=data.hclust, mboot=mboot, distance=distance, seed=seed)
	
    return(result)
  }

#bp = a list of all the bp values for a particular clade 
msfit <- function(bp, r, nboot) {

  if(length(bp) != length(r))
    stop("bp and r should have the same length")
	
  nboot <- rep(nboot, length=length(bp)) 

  use <- bp > 0 & bp < 1 #find all bp with values between 0 and 1

  p <- se <- c(0,0); names(p) <- names(se) <- c("au", "bp")
  coef <- c(0,0); names(coef) <- c("v", "c")

  a <- list(p=p, se=se, coef=coef, df=0, rss=0, pchi=0); class(a) <- "msfit"

  if(sum(use) < 2) { #are there at least two valid bp values
    # if(mean(bp) < .5) a$p[] <- c(0, 0) else a$p[] <- c(1, 1)
    if(mean(bp) < .5) a$p[] <- c(0, bp[r==1.0]) else a$p[] <- c(1, bp[r==1.0])
    return(a)
  }

  bp <- bp[use]; r <- r[use]; nboot <- nboot[use] #get only the bp that had values greater then 0 and less then 1

  zz <- -qnorm(bp) #find where the bp values lie phi inverse

  vv <- ((1 - bp) * bp) / (dnorm(zz)^2 * nboot)
  a$use <- use; a$r <- r; a$zz <- zz

  X   <- cbind(sqrt(r), 1/sqrt(r)); dimnames(X) <- list(NULL, c("v","c"))
  
  fit <- lsfit(X, zz, 1/vv, intercept=FALSE) #fit the curve
  a$coef <- coef <- fit$coef #get the coefficents
  
  h.au <- c(1, -1); h.bp <- c(1, 1)
  
  z.au <- drop(h.au %*% coef); z.bp <- drop(h.bp %*% coef) #%*% is matrix multiplication  au is v - c bp is v + c

  a$p["au"] <- pnorm(-z.au); a$p["bp"] <- pnorm(-z.bp) #phi
  
  V <- solve(crossprod(X, X/vv))
  vz.au <- drop(h.au %*% V %*% h.au); vz.bp <- drop(h.bp %*% V %*% h.bp)
  a$se["au"] <- dnorm(z.au) * sqrt(vz.au); a$se["bp"] <- dnorm(z.bp) * sqrt(vz.bp)
  a$rss <- sum(fit$residual^2/vv)
  
  if((a$df <- sum(use) - 2) > 0) {
    a$pchi <- pchisq(a$rss, lower.tail=FALSE, df=a$df)
  }
  else a$pchi <- 1.0
  return(a)
}

plot.msfit <- function(x, curve=TRUE, main=NULL, sub=NULL, xlab=NULL, ylab=NULL, ...)
{
  if(is.null(main)) main="Curve fitting for multiscale bootstrap resampling"
  if(is.null(sub))
    {
      sub  <- paste("AU = ", round(x$p["au"], digits=2),
                    ", BP = ", round(x$p["bp"], digits=2),
                    ", v = ", round(x$coef["v"], digits=2),
                    ", c = ", round(x$coef["c"], digits=2),
                    ", pchi = ", round(x$pchi, digits=2))
    }
  if(is.null(xlab)) xlab=expression(sqrt(r))
  if(is.null(ylab)) ylab=expression(z-value)
  
  a <- sqrt(x$r); b <- x$zz
  
  if(!is.null(a) && !is.null(b)) {
    plot(a, b, main=main, sub=sub, xlab=xlab, ylab=ylab, ...)
    if(curve) lines(x, ...)
  }
  else if (!is.null(a)){
    plot(0, 0, main=main, sub=sub, xlab=xlab, ylab=ylab,
         type="n", xaxt="n", yaxt="n", ...)
    a <- text(mean(a), 0, "No fitting")
  }
}

lines.msfit <- function(x, col=2, lty=1, ...) {
  v <- x$coef["v"]; c <- x$coef["c"]
  curve(v * x + c / x, add=TRUE, col=col, lty=lty)
}

summary.msfit <- function(object, digits=3, ...) {
  cat("\nResult of curve fitting for multiscale bootstrap resampling:\n\n")

  cat("Estimated p-values:\n")
  pv <- data.frame(object$p, object$se)
  names(pv) <- c("Estimate", "Std. Error"); row.names(pv) <- c("au", "bp")
  print(pv, digits=digits); cat("\n")

  cat("Estimated coefficients:\n")
  coef <- object$coef
  print(coef, digits=digits); cat("\n")

  cat(paste("Residual sum of squares: ", round(object$rss,digits=digits)),
      ",   p-value: ", round(object$pchi, digits=digits),
      " on ", object$df, " DF\n\n", sep="")

}
seplot <- function(object, type=c("au", "bp"), identify=FALSE,
                   main=NULL, xlab=NULL, ylab=NULL, ...)
  {
    if(!is.na(pm <- pmatch(type[1], c("au", "bp")))) {
      wh <- c("au", "bp")[pm]
      
      if(is.null(main))
        main <- "p-value vs standard error plot"
      if(is.null(xlab))
        xlab <- c("AU p-value", "BP value")[pm]
      if(is.null(ylab))
        ylab <- "Standard Error"
      
      plot(object$edges[,wh], object$edges[,paste("se", wh, sep=".")],
           main=main, xlab=xlab, ylab=ylab, ...)
      if(identify)
        identify(x=object$edges[,wh], y=object$edges[,paste("se", wh, sep=".")],
                 labels=row.names(object$edges))
    }
    else stop("'type' should be \"au\" or \"bp\".")
  }
