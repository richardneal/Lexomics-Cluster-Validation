currentNode = 1

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
plot.easyhclust <- function(x, filename = NULL, float=0.01,
                         col.pv=c(2,3,8), cex.pv=0.8, font.pv=NULL,
                         col=NULL, cex=NULL, font=NULL, lty=NULL, lwd=NULL,
                         main="", sub=NULL, xlab=NULL, height=800, width=800, specialLabels=NULL, ...)
{
  if(.Platform$OS.type == "windows")
  {
  	if(!is.null(filename))
  	{
		png(paste(filename, ".png", sep=""), width=width, height=height)
  	}
	
	else
	{
		windows(width=width, height=height)
	}
  }
  
  else if(.Platform$OS.type == "unix")
  {
	if(!is.null(filename))
 	{
		png(paste(filename, ".png", sep=""), width=width, height=height, type="Xlib")
	}

	else
	{
		X11(width=floor(width/96), height=floor(height/96), type="Xlib") #X11 specifies window size in inches for some bizare reason
									 #I'm not certain of the correct conversion factor but this is my best
									 #guess 
	}
  }
  
  metaTable <- x$metaTable #get metadata out of pvclust object
  
  #The line describing the dendrogram needs to be changed depending on if bp values are being displayed or not

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
	
  if(!is.null(filename)) #if writing to a file close the connection
  {
	dev.off()
  }
}

easyhclust <- function(inputFile, highlightFileName = NULL, height = 800, width = 800, outputFileName = NULL, metadata = FALSE, main = "", input.Transposed = TRUE)
{
	input.data <- read.table(as.character(inputFile), header=T,
		comment.char="", row.names=1, sep="\t", quote="")

	metaTable <- NULL
		
	#if there is metadata
	if(metadata)
	{
		metaTable <- input.data[,-(1:136)] #put metadata into a seperate table
		input.data <- input.data[,1:136] #keep the non metadata in the original table
	}
	
	if(input.Transposed)
	{
		data <- t(input.data)
	}	
	
	colSums <- apply(data, 2, sum) #each example/observation/object is one column, so find the sums of the columns
	denoms <- matrix(rep(colSums, dim(data)[1]), byrow=T, ncol=dim(data)[2]) #compute matrix to divide current matrix by to normalize matrix. Each entry in a column is the sum of the column
	relFreq <- data/denoms
		
	relFreq <- t(relFreq)	

	distance <- dist(relFreq, method="euclidean")
	data.hclust <- hclust(distance, method="average")
	
	easyhclust <- list(hclustObj = data.hclust, metaTable = metaTable)
	class(easyhclust) <- "easyhclust"
	
	specialLabels = NULL

	if(!is.null(highlightFileName)) #if a file of label names to highlight was specified
	{
		specialLabels <- scan(highlightFileName, what = "character", sep = "\n")
	}
	
	plot(easyhclust, filename = outputFileName, height = height, width = width, specialLabels = specialLabels, main = main)
	
	return(easyhclust)
}

result <- easyhclust("Thesis Defense Example.tsv", outputFileName = "ThesisDefense2")