hc2axes <- function(x)
{
  A <- x$merge # (n,n-1) matrix
  n <- nrow(A) + 1
  x.axis <- c()
  y.axis <- x$height
  
  x.tmp  <- rep(0,2)
  zz     <- match(1:length(x$order),x$order)

    for(i in 1:(n-1)) {
        ai <- A[i,1]

        if(ai < 0)
          x.tmp[1] <- zz[-ai]
        else
          x.tmp[1] <- x.axis[ai]
        
        ai <- A[i,2]
        
        if(ai < 0)
          x.tmp[2] <- zz[-ai]
        else
          x.tmp[2] <- x.axis[ai]

        x.axis[i] <- mean(x.tmp)
      }
  
  return(data.frame(x.axis=x.axis,y.axis=y.axis))
}

hc2split <- function(x)
  {
    A <- x$merge # (n-1,n) matrix
    n <- nrow(A) + 1
    B <- list()

    for(i in 1:(n-1)){
        ai <- A[i,1]
        
        if(ai < 0)
          B[[i]] <- -ai
        else
          B[[i]] <- B[[ai]]        
        
        ai <- A[i,2]
        
        if(ai < 0)
          B[[i]] <- sort(c(B[[i]],-ai))
        else
          B[[i]] <- sort(c(B[[i]],B[[ai]]))
      }

    CC <- matrix(rep(0,n*(n-1)),nrow=(n-1),ncol=n)
    
    for(i in 1:(n-1)){
        bi <- B[[i]]
        m <- length(bi)
        for(j in 1:m)
          CC[i,bi[j]] <- 1
      }

    split <- list(pattern=apply(CC,1,paste,collapse=""), member=B)
    
    return(split)
  }

pvclust.node <- function(x, r,...) #this does the bootstraping for a single node
  {
#    require(pvclust)
    mboot.node <- lapply(r, boot.hclust, nboot=x, ...) #do the bootstraping once for each r value
    return(mboot.node)
  }

boot.hclust <- function(r, data, object.hclust, method.dist, use.cor,
                        method.hclust, nboot, store, weight=F, storeCop, copDistance, normalize)
{

  n     <- nrow(data) #get the number of rows (each row contains a single word)
  size  <- round(n*r, digits=0) #calculate the number of rows to resample
  if(size == 0)
    stop("invalid scale parameter(r)")
  r <- size/n #recalculate r to match exactly size/n instead of approxiametly 

  pattern   <- hc2split(object.hclust)$pattern #get the pattern (the contents of each clade) from the original hclust object
  edges.cnt <- table(factor(pattern)) - table(factor(pattern)) #create a table that shows the number of times each clade is formed in the bootstraped clusters
  st <- list() #list that will store the hclust objects creatd
  st <- list() #list that will store the hclust objects create
  stc <- list() #list storing cophenetic correlations

  # bootstrap start
  rp <- as.character(round(r,digits=2)); if(r == 1) rp <- paste(rp,".0",sep="")
  cat(paste("Bootstrap (r = ", rp, ")... ", sep=""))
  w0 <- rep(1,n) # equal weight
  na.flag <- 0
  
  for(i in 1:nboot){
    if(weight && r>10) {  ## <- this part should be improved
      w1 <- as.vector(rmultinom(1,size,w0)) # resampled weight
	  if(normalize)
	  {
		print("Code for normalizing resampled weight is not in yet. Look at code for nonresampled weight in boot.hclust for example of how to implement")
	  }
	  
      suppressWarnings(distance <- distw.pvclust(data,w1,method=method.dist,use.cor=use.cor))
    } else {
		bootStrapData <- matrix(data = 0, nrow = nrow(data), ncol = ncol(data), dimnames=list(rownames(data),colnames(data)))
		#print(tcounts)
		#row.names(tcounts) <- rownames(data)
		#print(tcounts)
		#names(tcounts) <- colnames(data)		
		#print(tcounts)
		for(i in 1:ncol(data)){
                                                         
		 wordlist <- rep(rownames(data),data[,i]) #gets all words in the chunk ordered by word (each word appears count times)
		 size <- length(wordlist) * r
         tmpindex <- sample(1:sum(data[,i]), size, replace = TRUE) #size will equal sum(data[,i]) when r = 1
         counts <- table(wordlist[tmpindex]) #returns counts for new sample
		 
		 #print(counts)
		 
		 for(j in 1:nrow(data))
		 {
			bootStrapData[j,i] = counts[row.names(bootStrapData)[j]]
			if(is.na(bootStrapData[j,i])) #if there was no occurances of this word in new sample
			{
				bootStrapData[j,i] = 0
			}
			
		 }
		
		#try resample by clades. Set a cutoff for what clades to use and then sum up counts of all words inside
                #each clade. Resample each chunk in the clade using those total sums of the clade
	  }

      #smpl <- sample(1:n, size, replace=TRUE) #creates a index vector with each element being the index of the row chosen
	  
	  #bootStrapData = data[smpl,] #get the new bootStrap values using the index vector
      #print(bootStrapData)
	  if(normalize)
	  {
	  #normalize the new data
		  colSums <- apply(bootStrapData, 2, sum) #each example/observation/object is one column, so find the sums of the columns
		  denoms <- matrix(rep(colSums, dim(bootStrapData)[1]), byrow=T, ncol=dim(bootStrapData)[2]) #compute matrix to divide current matrix by to normalize matrix. Each entry in a column is the sum of the column
		  bootStrapData <- bootStrapData/denoms
	  }
	  
      suppressWarnings(distance  <- dist.pvclust(bootStrapData,method=method.dist,use.cor=use.cor)) #ceate the new distance matrix
    }
    if(all(is.finite(distance))) { # check if distance is valid
      x.hclust  <- hclust(distance,method=method.hclust) #create the hclust object
      pattern.i <- hc2split(x.hclust)$pattern # get the contents of the clades
      edges.cnt <- edges.cnt + table(factor(pattern.i,  levels=pattern)) #add the clades identical to the clades in the original object to the count of the number of times those clades appeared
    } else {
      x.hclust <- NULL
	  na.flag <- 1
    }

    if(store) #if storing the hclust objects
      st[[i]] <- x.hclust


    if(storeCop)  #if storing cophenetic corelations, find and store the correlation now
    {
        bootCop <- cophenetic(x.hclust) #get the cophenetic distance matrix
        stc[[i]] <- cor(copDistance, bootCop) #find the correlation between the cophenetic distance matrix and the distance matrix from the original data.
    }

  }
  cat("Done.\n")
  # bootstrap done
  
  if(na.flag == 1)
	warning(paste("inappropriate distance matrices are omitted in computation: r = ", r), call.=FALSE)

  boot <- list(edges.cnt=edges.cnt, method.dist=method.dist, use.cor=use.cor,
               method.hclust=method.hclust, nboot=nboot, size=size, r=r, store=st, storeCop=stc) #store all the data in a list
  class(boot) <- "boot.hclust"
  
  return(boot) #return that list to the function call
}

#this function takes all the data computed and converts it into the final results
pvclust.merge <- function(data, object.hclust, mboot, distance){
  
  pattern <- hc2split(object.hclust)$pattern #get the contents of the clades in the original hclust objects
 
  #get the following data out of the list 
  r     <- unlist(lapply(mboot,"[[","r"))
  nboot <- unlist(lapply(mboot,"[[","nboot"))
  store <- lapply(mboot,"[[", "store")
  storeCop <-lapply(mboot,"[[", "storeCop")
  
  rl <- length(mboot) #get the number of r values
  ne <- length(pattern) #get the number of clades
  
  edges.bp <- edges.cnt <- data.frame(matrix(rep(0,ne*rl),nrow=ne,ncol=rl)) #Creates two big table, with each row being a clade, and each column being a r value. A individual entry denotes either the number of times that clade occured or the bp value for that clade for a particular r value.
  row.names(edges.bp) <- pattern
  names(edges.cnt) <- paste("r", 1:rl, sep="")

  for(j in 1:rl) {
    edges.cnt[,j] <- as.vector(mboot[[j]]$edges.cnt) #holds the number of times each clade is formed for this particular r value
	#print(edges.cnt[,j])
    edges.bp[,j]  <- edges.cnt[,j] / nboot[j]  #holds the number of times each clade is formed divided by the total number of bootstraps for this particular r value
	#print(edges.bp[,j])
  }
  
  ms.fitted <- lapply(as.list(1:ne),									  #does the function one per clade that computes bp and au values
                      function(x, edges.bp, r, nboot){
                        msfit(as.vector(t(edges.bp[x,])), r, nboot)},
                      edges.bp, r, nboot)
  class(ms.fitted) <- "mslist"
  
  p    <- lapply(ms.fitted,"[[","p") #p is a list of objects one per clade that hold both the au and bp values of that clade
  se   <- lapply(ms.fitted,"[[","se") #se holds a different set of bp and au values????
  coef <- lapply(ms.fitted,"[[","coef")
  au    <- unlist(lapply(p,"[[","au"))
  bp    <- unlist(lapply(p,"[[","bp"))
  se.au <- unlist(lapply(se,"[[","au"))
  se.bp <- unlist(lapply(se,"[[","bp"))
  v     <- unlist(lapply(coef,"[[","v"))
  cc    <- unlist(lapply(coef,"[[","c"))
  pchi  <- unlist(lapply(ms.fitted,"[[","pchi"))
  
  #store all the data in a single table
  edges.pv <- data.frame(au=au, bp=bp, se.au=se.au, se.bp=se.bp,
                         v=v, c=cc, pchi=pchi)

  row.names(edges.pv) <- row.names(edges.cnt) <- 1:ne #names by default seem to just be numbers 1-n with n being number of clades

  #combine all the data into a list
  result <- list(hclust=object.hclust, edges=edges.pv, count=edges.cnt,
                 msfit=ms.fitted, nboot=nboot, r=r, store=store, storeCop=storeCop, distance=distance)

  class(result) <- "pvclust"
  return(result) #return that list
}

dist.pvclust <- function(x, method="euclidean", use.cor="pairwise.complete.obs")
{
  if(!is.na(pmatch(method,"correlation"))){
    res <- as.dist(1 - cor(x, method="pearson", use=use.cor))
    attr(res,"method") <- "correlation"
    return(res)
  }
  else if(!is.na(pmatch(method,"abscor"))){
    res <- as.dist(1 - abs(cor(x,method="pearson",use=use.cor)))
    attr(res,"method") <- "abscor"
    return(res)
  }
  else if(!is.na(pmatch(method,"uncentered"))){
    if(sum(is.na(x)) > 0){
      x <- na.omit(x)
      warning("Rows including NAs were omitted")
    }
    x  <- as.matrix(x)
    P  <- crossprod(x)
    qq <- matrix(diag(P),ncol=ncol(P))
    Q  <- sqrt(crossprod(qq))
    res <- as.dist(1 - P/Q)
    attr(res,"method") <- "uncentered"
    return(res)
  }
  else
    dist(t(x),method)
}


corw <- function(x,w,
                 use=c("all.obs","complete.obs","pairwise.complete.obs")
                 ) {
  if(is.data.frame(x)) x <- as.matrix(x)
  x <- x[w>0,,drop=F]
  w <- w[w>0]

  n <- nrow(x) # sample size
  m <- ncol(x) # number of variables
  if(missing(w)) w <- rep(1,n)
  r <- matrix(0,m,m,dimnames=list(colnames(x),colnames(x)))
  diag(r) <- 1
  use <- match.arg(use)

  pairu <- F
  if(use=="all.obs") {
    u <- rep(T,n)
  } else if(use=="complete.obs") {
    u <- apply(x,1,function(y) !any(is.na(y)))
  } else if(use=="pairwise.complete.obs") {
    pairu <- T
    ux <- is.finite(x)
  } else stop("unknown use")
  
  for(i in 1+seq(length=m-1)) {
    for(j in seq(length=i-1)) {
      if(pairu) u <- ux[,i] & ux[,j]
      wu <- w[u]; xi <- x[u,i]; xj <- x[u,j]
      ws <- sum(wu)
      if(ws > 1e-8) {
        xi <- xi - sum(wu*xi)/ws
        xj <- xj - sum(wu*xj)/ws
        vxi <- sum(wu*xi*xi)/ws
        vxj <- sum(wu*xj*xj)/ws
        if(min(vxi,vxj) > 1e-8)  {
          vxij <- sum(wu*xi*xj)/ws
          rij <- vxij/sqrt(vxi*vxj)
        } else {
          rij <- 0
        }
      } else {
        rij <- 0
      }
      r[i,j] <- r[j,i] <- rij
    }
  }
  r
}

### calculate distance by weight
distw.pvclust <- function(x,w,method="correlation", use.cor="pairwise.complete.obs")
{
  if(!is.na(pmatch(method,"correlation"))){
    res <- as.dist(1 - corw(x,w, use=use.cor))
    attr(res,"method") <- "correlation"
    return(res)
  }
  else if(!is.na(pmatch(method,"abscor"))){
    res <- as.dist(1 - abs(corw(x,w, use=use.cor)))
    attr(res,"method") <- "abscor"
    return(res)
  }
  stop("wrong method")
}

createRangeList <- function(x)
{
	numChunks <- ncol(x)
	numWords <- nrow(x)
	rangeList <- list()
	
	#set up one sublist for each list
	for(i in 1:numChunks)
	{
		wordRangeStart <- 1
		wordsSeen <- 1
		chunkRangeList <- list()
		for(j in 1:numWords)
		{
			if(x[j,i] != 0)
			{
				wordRange = list(wordIndex = j, start = wordRangeStart, end = (wordRangeStart + x[j,i] - 1))
				chunkRangeList[[wordsSeen]] <- wordRange
				wordsSeen <- wordsSeen + 1
				wordRangeStart <- wordRangeStart + x[j,i]
			}
		}
		rangeList[[i]] <- chunkRangeList
	}
	
	print(x)
	#print(rangeList)
	
	print(findIndexInRangeList(rangeList, 14, 2))
	
	return(rangeList)
}

findIndexInRangeList <- function(rangeList, wordNum, chunk)
{
	left <- 1
	right <- length(rangeList[[chunk]])
	
	while(left <= right)
	{
		m <- floor((left + right) / 2)
		
		if(rangeList[[chunk]][[m]]$end < wordNum)
		{
			left <- m + 1
		}
		
		else if(rangeList[[chunk]][[m]]$start > wordNum)
		{
			right <- m - 1
		}
		
		else
		{
			return(rangeList[[chunk]][[m]]$wordIndex)
		}

	}
	return(-1)
}