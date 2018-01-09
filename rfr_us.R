rfrus <- function(X, MinParent=1, trees=100, MaxDepth="inf", bagging=.2, replacement=TRUE, FUN=makeA, options=c(ncol(X), round(ncol(X)^.5),1L, 1/ncol(X)), COOB=TRUE, Progress=TRUE){
    #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    # rfr builds a randomer classification forest structure made up of a list
    # of trees.  This forest is randomer because each node is rotated before 
    # being split (as described by Tyler Tomita).  Each tree is made of nodes
    # which are numbered from left to right within a particular level (depth)
    # of a tree.  The loop over nodes when considering splits is made in the
    # same order as the numbering of the nodes. That is, nodes are traversed
    # across a single level from left to right then moves to the left-most node
    # of the next level.
    #
    #  INPUT:
    #
    # X is an n-by-p matrix, where rows represent observations and columns
    # represent features
    #
    # Y is an n-by-1 array of integer class labels. Unique class labels
    # must be contiguous and start from 1 (e.g. [0,1,2] is not okay;
    # neither is [1,3,4])
    #
    # MinParent is an integer specifying the minimum number of observations
    # a node must have in order for an attempt to split to be made.  Lower
    # values may lead to overtraining and increased training time.
    #
    # trees is the number of trees that will be in the forest.
    #
    # MaxDepth is the maximum depth that a tree can grow to.  If set to "inf"
    # then there is no maximum depth.  If set to 0 then a maximum depth is
    # calculated based on the number of classes and number of samples provided.
    #
    # bagging is the percentage of training data to withhold during each
    # training iteration.  If set to 0 then the entire training set is used
    # during every iteration.  The withheld portion of the training data
    # is used to calculate OOB error for the tree.
    #
    # ClassCt is the number of different classes in Y.  It is calculated 
    # in the calling function to prevent recalculation by each forked function 
    # when in parallel.
    #
    # FUN is the function used to create the projection matrix.  The matrix
    # returned by this function should be a p-by-u matrix where p is the
    # number of columns in the input matrix X and u is any integer > 0.
    # u can also vary from node to node.
    #
    # options is a list of inputs to the user provided projection matrix
    # creation function -- FUN.
    #
    # COOB is a boolean that determines whether or not OOB error is calculated.
    # If bagging equals zero then COOB is ignored.  If bagging does not equal 
    # zero and COOB is TRUE then OOB is calculated and printed to the screen.
    #
    # Progress is a boolean.  When true a progress marker is printed to the 
    # screen every time a tree is grown.  This is useful for large input.
    #
    # OUTPUT:
    #
    # A forest construct made up of trees.  This forest can be used to make 
    # predictions on new inputs.  When COOB=TRUE then the output is a list
    # containing $forest and $OOBmat.  $forest is the forest structure and
    # OOBmat is the OOB error for each tree.
    #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


    #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    # Predefine variables to prevent recreation during loops
    #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    forest <- vector("list",trees)
    BV <- NA # vector in case of ties
    BS <- NA # vector in case of ties
    MaxDeltaI <- 0
    nBest <- 1L
    BestIdx <-0L 
    BestVar <-0L 
    BestSplitIdx<-0L 
    BestSplitValue <- 0
    w <- nrow(X)
    p <- ncol(X)
    perBag <- (1-bagging)*w
    Xnode<-double(w) # allocate space to store the current projection
    SortIdx<-integer(w) 
    if(object.size(X) > 1000000){
        OS<-TRUE
    }else{
        OS<-FALSE
    }

    # Calculate the Max Depth and the max number of possible nodes
    if(MaxDepth == "inf"){
        StopNode <- 2L*w #worst case scenario is 2*(w/(minparent/2))-1
        MaxNumNodes <- 2L*w # number of tree nodes for space reservation
    }else{
        if(MaxDepth==0){
            MaxDepth <- ceiling(log2(w))
        }
        StopNode <- 2L^(MaxDepth)
        MaxNumNodes <- 2L^(MaxDepth+1L)  # number of tree nodes for space reservation
    }

    CutPoint <- double(MaxNumNodes)
    Children <- matrix(data = 0L, nrow = MaxNumNodes,ncol = 2L)
    NDepth <- integer(MaxNumNodes)
    matA <- vector("list", MaxNumNodes) 
    Assigned2Node<- vector("list",MaxNumNodes) 
    Assigned2Leaf <- vector("list", MaxNumNodes)
    ind <- double(w)
    min_error <- Inf
    #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    #                            Start tree creation
    #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    for(treeX in 1:trees){
        # intialize values for new tree before processing nodes
        CutPoint[] <- 0
        Children[] <- 0L
        NDepth[]<- 0L #delete this?
        NDepth[1]<-1L
        CurrentNode <- 1L
        NextUnusedNode <- 2L
        NodeStack <- 1L
        highestParent <- 1L
    Assigned2Leaf <- vector("list", MaxNumNodes)
        ind[] <- 0L
        # Determine bagging set 
        # Assigned2Node is the set of row indices of X assigned to current node
        if(bagging != 0){
            if(replacement){
                ind<-sample(1:w, w, replace=TRUE)
                Assigned2Node[[1]] <- ind
            }else{
                ind[1:perBag] <- sample(1:w, perBag, replace = FALSE)
                Assigned2Node[[1]] <- ind[1:perBag]        
            }
        }else{
            Assigned2Node[[1]] <- 1:w        
        }
        # main loop over nodes
        while (CurrentNode < NextUnusedNode && CurrentNode < StopNode){
            # determine working samples for current node.
            NodeRows <- Assigned2Node[CurrentNode] 
            Assigned2Node[[CurrentNode]]<-NA #remove saved indexes
            NdSize <- length(NodeRows[[1L]]) #determine node size

            # create projection matrix (sparseM) by calling the custom function FUN
            sparseM <- FUN(options)
            #isolate objective function
                # if node is impure and large enough then attempt to find good split
                if (NdSize < MinParent || NDepth[CurrentNode]==MaxDepth || NextUnusedNode+1L >= StopNode || NdSize == 1){
Assigned2Leaf[[CurrentNode]] <- NodeRows[[1L]]
                    NodeStack <- NodeStack[-1L]
                    CurrentNode <- NodeStack[1L]
                    if(is.na(CurrentNode)){
                        break
                    }
                    next 
                }
                min_error <- Inf
                    nBest <- 1L
                for(q in unique(sparseM[,2])){
                    #Project input into new space
                    lrows <- which(sparseM[,2]==q)
                    Xnode[1:NdSize] <- X[NodeRows[[1L]],sparseM[lrows,1], drop=FALSE]%*%sparseM[lrows,3, drop=FALSE]
                    #Sort the projection, Xnode, and rearrange Y accordingly
                    SortIdx[1:NdSize] <- sort(Xnode[1:NdSize] )

                    ##################################################################
                    #                    Find Best Split
                    ##################################################################

                    # initialize variables for loop through projection
# run the first iteration outside of the loop
                    sumLeft <- 0
                    sumRight <- sum(SortIdx[1:NdSize])
                        errLeft <- 0
                        errRight <-0 
                       meanLeft <- 0
                       meanRight <- 0
                       errCurr <- 0
                    for (m in 1:(NdSize-1)){
                        sumLeft <- sumLeft + SortIdx[m]
                        sumRight <- sumRight - SortIdx[m]
                        meanLeft <- sumLeft/(m)
                        meanRight <- sumRight/(NdSize-m)
                        errLeft <-sum((SortIdx[1:m]-meanLeft)^2) 
                        errRight <-sum((SortIdx[(m+1):NdSize]-meanRight)^2) 

                        errCurr <- errLeft + errRight
                        # Determine if this split is currently the best option
                        if (errCurr <= min_error){
                            # Save current best DeltaI
                            if (min_error != errCurr){
                                min_error <- errCurr
                                nBest <- 1L
                            }else{
                                # Save all DeltaI equal to current max DeltaI
                                nBest <- nBest + 1L
                            }
                            BV[nBest] <- q
                            BS[nBest] <- SortIdx[m]+SortIdx[m+1]
                        }
                    }

                    sumLeft <- sumLeft +SortIdx[NdSize]
                        meanLeft <- sumLeft/NdSize
errLeft <-sum((SortIdx[1:NdSize]-meanLeft)^2)
                    errCurr <- errLeft

                    if (errCurr <= min_error){
                        # Save current best DeltaI
                        if (min_error != errCurr){
                            min_error <- errCurr
                            nBest <- 1L
                        }else{
                            # Save all DeltaI equal to current max DeltaI
                            nBest <- nBest + 1L
                        }
                        BV[nBest] <- q
                        BS[nBest] <- SortIdx[NdSize]
                    }
                }#end loop through projections.

                # If there were multiple best splits then randomly choose
                # between the best.
                if (nBest > 1L){
                    # Break ties at random
                    BestIdx <- ceiling(runif(1,0,nBest))
                    BestVar <- BV[BestIdx]
                    BestSplitIdx <- BS[BestIdx]/2
                }else{
                    BestVar <- BV[1L]
                    BestSplitIdx <- BS[1L]/2
                }
                                # Recalculate the best projection
                lrows<-which(sparseM[,2L]==BestVar)
                Xnode[1:NdSize]<-X[NodeRows[[1L]],sparseM[lrows,1], drop=FALSE]%*%sparseM[lrows,3, drop=FALSE]
            
            # determine split value as mean of values on either side of split
            BestSplitValue <- BestSplitIdx

            # find which child node each sample will go to and move
            # them accordingly
# changed this from <= to < just in case best split split all values
            MoveLeft <- Xnode[1:NdSize]  < BestSplitValue
            numMove <- sum(MoveLeft)

            if (is.null(numMove)){
                print("numMove is null")
                flush.console()
            }
            if(is.na(numMove)){
                print("numMove is na")
                flush.console()
            }
            #Check to see if a split occured, or if all elements being moved one direction.
            if(numMove!=0L && numMove!=NdSize){
                # Move samples left or right based on split
                Assigned2Node[[NextUnusedNode]] <- NodeRows[[1L]][MoveLeft]
                Assigned2Node[[NextUnusedNode+1L]] <- NodeRows[[1L]][!MoveLeft]

                #highest Parent keeps track of the highest needed matrix and cutpoint
                # this reduces what is stored in the forest structure
                if(CurrentNode>highestParent){
                    highestParent <- CurrentNode
                }
                # Determine children nodes and their attributes
                Children[CurrentNode,1L] <- NextUnusedNode
                Children[CurrentNode,2L] <- NextUnusedNode+1L
                NDepth[NextUnusedNode]=NDepth[CurrentNode]+1L
                NDepth[NextUnusedNode+1L]=NDepth[CurrentNode]+1L
                # Pop the current node off the node stack
                # this allows for a breadth first traversal
Assigned2Leaf[[CurrentNode]] <- NodeRows[[1L]]
                NodeStack <- NodeStack[-1L]
                NodeStack <- c(NextUnusedNode, NextUnusedNode+1L, NodeStack)
                NextUnusedNode <- NextUnusedNode + 2L
                # Store the projection matrix for the best split
                matA[[CurrentNode]] <- as.integer(t(sparseM[which(sparseM[,2]==BestVar),c(1,3)]))
                CutPoint[CurrentNode] <- BestSplitValue
            }else{
                # There wasn't a good split so ignore this node and move to the next
                NodeStack <- NodeStack[-1L]
            }
            # Store ClassProbs for this node.
            # Only really useful for leaf nodes, but could be used instead of recalculating 
            # at each node which is how it is currently.
            CurrentNode <- NodeStack[1L]
            if(is.na(CurrentNode)){
                break
            }
        }
        #If input is large then garbage collect prior to adding onto the forest structure.
        if(OS){
            gc()
        }
        # save current tree structure to the forest
        if(bagging!=0 && COOB){
            forest[[treeX]] <- list("CutPoint"=CutPoint[1:highestParent],"Children"=Children[1L:(NextUnusedNode-1L),,drop=FALSE], "matA"=matA[1L:highestParent], "ALeaf"=Assigned2Leaf[1L:(NextUnusedNode-1L)])
        }else{
            forest[[treeX]] <- list("CutPoint"=CutPoint[1:highestParent],"Children"=Children[1L:(NextUnusedNode-1L),,drop=FALSE], "matA"=matA[1L:highestParent], "ALeaf"=Assigned2Leaf[1L:(NextUnusedNode-1L)])
        }
        if(Progress){
            cat("|")
            flush.console()
        }
    }
    return(forest)
}


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#                      Default option to make projection matrix 
#
# this is the randomer part of random forest. The sparseM 
# matrix is the projection matrix.  The creation of this
# matrix can be changed, but the nrow of sparseM should
# remain p.  The ncol of the sparseM matrix is currently
# set to mtry but this can actually be any integer > 1;
# can even greater than p.
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
makeA <- function(options){
    p <- options[[1L]]
    d <- options[[2L]]
    method <- options[[3L]]
    if(method == 1L){
        rho<-options[[4L]]
        nnzs <- round(p*d*rho)
        sparseM <- matrix(0L, nrow=p, ncol=d)
        sparseM[sample(1L:(p*d),nnzs, replace=F)]<-sample(c(1L,-1L),nnzs,replace=T)
    }
    #The below returns a matrix after removing zero columns in sparseM.
    ind<- which(sparseM!=0,arr.ind=TRUE)
    return(cbind(ind,sparseM[ind]))        
}


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#                     Create Distance Matrix
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
dist <- function(X, Forest, maxDepth=0){
    n <- nrow(X)
    dist <- matrix(0, nrow=n, ncol=n)

    numT <- length(Forest)
    currBin <- integer(n)
    if (maxDepth==0){
    for(j in 1:numT){
        for(i in 1:n){
            currentNode <- 1L
            while(Forest[[j]]$Children[currentNode]!=0L){
                s<-length(Forest[[j]]$matA[[currentNode]])/2
                rotX <-sum(Forest[[j]]$matA[[currentNode]][(1:s)*2]*X[i,Forest[[j]]$matA[[currentNode]][(1:s)*2-1]])
                if(rotX<=Forest[[j]]$CutPoint[currentNode]){
                    currentNode <- Forest[[j]]$Children[currentNode,1L]
                }else{
                    currentNode <- Forest[[j]]$Children[currentNode,2L]
                }
            }
            dist[Forest[[j]]$ALeaf[[currentNode]]] <- dist[Forest[[j]]$ALeaf[[currentNode]]] + 1
                for(z in 1:(i-1)){
                    if(currBin[z] == currentNode){
                        dist[i,z] <- dist[i,z]+1
                        dist[z,i] <- dist[z,i]+1
                    }
                }
            dist[i,i] <- dist[i,i]+1
        }
    }
    }else{
for(j in 1:numT){
        for(i in 1:n){
            currentNode <- 1L
            depth <- 1L
            while(Forest[[j]]$Children[currentNode]!=0L && depth <= maxDepth){
                s<-length(Forest[[j]]$matA[[currentNode]])/2
                rotX <-sum(Forest[[j]]$matA[[currentNode]][(1:s)*2]*X[i,Forest[[j]]$matA[[currentNode]][(1:s)*2-1]])
                if(rotX<=Forest[[j]]$CutPoint[currentNode]){
                    currentNode <- Forest[[j]]$Children[currentNode,1L]
                }else{
                    currentNode <- Forest[[j]]$Children[currentNode,2L]
                }
                depth <- depth+1L
            }
            currBin[i] <- currentNode
            if(i>1){
                for(z in 1:(i-1)){
                    if(currBin[z] == currentNode){
                        dist[i,z] <- dist[i,z]+1
                        dist[z,i] <- dist[z,i]+1
                    }
                }
            }
            dist[i,i] <- dist[i,i]+1
        }
    }
    }
    return(dist)        
}



#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#                    Find Potential Nearest Neighbors Vector
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
distNN <- function(y, X, Forest, maxDepth=0){
    dist <- integer(nrow(X))

    numT <- length(Forest)
    if (maxDepth==0){
        maxDepth <- Inf
    }
    for(j in 1:numT){
            currentNode <- 1L
            depth <- 1L
            while(Forest[[j]]$Children[currentNode]!=0L && depth <= maxDepth){
                s<-length(Forest[[j]]$matA[[currentNode]])/2
                rotX <-sum(Forest[[j]]$matA[[currentNode]][(1:s)*2]*y[Forest[[j]]$matA[[currentNode]][(1:s)*2-1]])
                if(rotX<=Forest[[j]]$CutPoint[currentNode]){
                    currentNode <- Forest[[j]]$Children[currentNode,1L]
                }else{
                    currentNode <- Forest[[j]]$Children[currentNode,2L]
                }
                depth <- depth+1L
            }
            dist[Forest[[j]]$ALeaf[[currentNode]]] <- dist[Forest[[j]]$ALeaf[[currentNode]]] + 1
    }
    return(dist) #this is the similarity vector 
}

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#                    Find Nearest Neighbors from similarity vector
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
distNNk <- function(y, X, sv, k, adder){
    index <- order(sv, decreasing=TRUE)
    simCount <- tabulate(sv)
    multiplier  <- adder
    if(sum(simCount) < multiplier+k){
        remainingNN <- sum(simCount)
        print("Not enough points.  Decrease search depth.")
        flush.console()
        return(NULL)
    }else{
    remainingNN <- multiplier+k
    }
    simLength <- length(simCount)
    NNindex <- NULL
    while (remainingNN >0){
        if (remainingNN >= simCount[simLength]){
            if(simCount[simLength]>0){
            NNindex <- c(NNindex, index[1:simCount[simLength]])
            index <- index[-(1:simCount[simLength])]
            remainingNN <- remainingNN-simCount[simLength]
            }
            simLength <- simLength -1
        }else{
           #NNorder <- order(sqrt(rowSums((y-X[index[1:simCount[simLength]],])^2))) 
           NNorder <- order(sqrt(rowSums(sweep(X[index[1:simCount[simLength]],],2,y)^2))) 
NNindex <- c(NNindex, index[NNorder[1:remainingNN]])
remainingNN = 0
        }
    }
    kNearest <- order(sqrt(rowSums(sweep(X[NNindex,],2,testSamples[1,])^2)))[1:k]

    
return(NNindex[kNearest])
}


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#                    Check K-Means 
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
CheckKmeans <- function(Y, Yp){
    uY <- length(unique(Y))
    classCt <- tabulate(Y, uY)

    class_order <- order(classCt, decreasing=TRUE)
    used_class <-NULL 
    curr_class <- NA
    class_error <- NA
    for(z in 1:uY){
        Cindex <- which(Y==class_order[z])
        subClassCt <- tabulate(Yp[Cindex], uY)
        subClass_order <- order(subClassCt, decreasing=TRUE)
        if(!is.null(used_class)){
            for(m in 1:uY){
                curr_class <- subClass_order[m]
                if(!(curr_class%in%used_class)){
                    break
                }
            }
            used_class <- c(used_class, curr_class)
        }else{
            curr_class <- subClass_order[1]
            used_class <- curr_class
        }

        class_error[z] <- subClassCt[curr_class]/classCt[class_order[z]]
    }
    print(class_error)
    oe <- sum(class_error*classCt[class_order])/length(Y)
    cat("the overall error is: ", oe, "\n")
    flush.console()
}


#############################Swiss Roll Code###################################
swissRoll <- function(n1, n2 = NULL, size = 6, dim3 = FALSE, rand_dist_fun = NULL, ...) {

### If n2 is NULL, then generate a balanced dataset of size 2*n1
     if (is.null(n2)) n2 <- n1
 xdim <- ifelse(dim3, 3, 2) 
  ### GROUP 1 
  # Generate Angles 
  rho <- runif(n1, 0, size*pi)
  # Create Swiss Roll
  g1x1 <- rho*cos(rho)
   g1x2 <- rho*sin(rho)
   
   ### GROUP 2
   # Generate Angles
   rho <- runif(n2, 0, size*pi)
    # Create Inverse Swiss Roll
    g2x1 <- -rho*cos(rho)
    g2x2 <- -rho*sin(rho)
     
     ### Generate the 3rd dimension
     if (dim3) {
          z_range <- range(c(g1x1, g1x2, g2x1, g2x2))
     x3 <- runif(n1 + n2, z_range[1], z_range[2])
      } 
     
     ### If needed random perturbation on the data
     ### please specify the random generation funciton in R to 'rand_dist_fun'
     ### and the corresponding parameters in '...'.
     ### For example, 
     ### rand_dist_fun = rnorm, mean = 0, sd = 0.2
     err <- matrix(0, n1 + n2, xdim)
      if (!is.null(rand_dist_fun)) err <- matrix(rand_dist_fun(xdim*(n1 + n2), ...), n1 + n2, xdim)
      
      ### Output the Swiss Roll dataset
      if (dim3) {
           out <- data.frame(y = c(rep(0:1, c(n1, n2))), x1 = c(g1x1, g2x1) + err[,1], x2 = c(g1x2, g2x2) + err[,2], x3 = x3 + err[,3])
       } else {
            out <- data.frame(y = c(rep(0:1, c(n1, n2))), x1 = c(g1x1, g2x1) + err[,1], x2 = c(g1x2, g2x2) + err[,2])
        }
       out
}
################################################################################
################################################################################

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#                     Hartigan's Method
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
findClusters <- function(nearnessMatrix, numClusters=3, numNearestNeighbors=10){
    q <- rep(0,numClusters)
    clusters <- vector("list", numClusters)
    numSamples <- nrow(nearnessMatrix)
    numNN <- numNearestNeighbors

    randomOrdering <- sample(1:numSamples, numSamples)
  #  randomOrdering <- 1:numSamples


    step <- floor(numSamples/numClusters)
    stepStart <- 1
    stepEnd <- stepStart+step
    for(z in 1:(numClusters-1)){
        clusters[[z]] <- randomOrdering[stepStart:stepEnd]
        stepStart <- stepEnd+1
        stepEnd <- stepStart+step
    }
    clusters[[numClusters]] <- randomOrdering[stepStart:numSamples]

    for(z in 1:numSamples){
        nearnessMatrix[z,z] <- 0
    }

    for(z in 1:numClusters){
        for(m in clusters[[z]]){
            biggestNN <- order(nearnessMatrix[m,], decreasing=TRUE)[1:numNN]
            q[z] <- q[z] + sum(biggestNN%in%clusters[[z]])
        }
    }
    print(paste("initial q", q))

    currQ <- rep(0,numClusters)
    for(p in 1:30){
        for(z in 1:numClusters){
            for(m in clusters[[z]]){
                biggestNN <- order(nearnessMatrix[m,], decreasing=TRUE)[1:numNN]
                for(k in 1:numClusters){
                    currQ[k] <- sum(biggestNN%in%clusters[[k]])
                }
                QOrder <- order(currQ, decreasing=TRUE)
                if(QOrder[1] != z){
                    q[z] <- q[z] - currQ[z]
                    q[QOrder[1]] <- q[QOrder[1]] + currQ[QOrder[1]]
clusters[[z]] <- clusters[[z]][-which(clusters[[z]]==m)]
clusters[[QOrder[1]]] <- c(clusters[[QOrder[1]]], m)
                }
            }
        }
    }

    print(paste("after 1", q))
    return(clusters)
}

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#                     Spectral Cluster
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
specN <- function(distMat, numClust){
    Y <- kmeans(distMat, numClust)$cluster
    return(Y)        
}

#############################################################################
require(compiler)
rfrus <- cmpfun(rfrus)
distNN <- cmpfun(distNN)


createSimilarityMatrix <- function(X, numTrees=100, K=10){
numberSamples <- nrow(X)
similarityMatrix <- matrix(0,nrow= numberSamples, ncol=numberSamples)

forest <- invisible(rfrus(X,trees=numTrees, MinParent=K))

for(z in 1:numberSamples){
    NN1 <- distNN(X[z,], X, forest)
    similarityMatrix[z,] = NN1
    for(q in 1:numberSamples){ #Why did I do this?
        if(NN1[q]==0){
            similarityMatrix[z,q]<-0
        }
    }
}
return(similarityMatrix)
}

