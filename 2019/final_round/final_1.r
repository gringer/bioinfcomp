setwd("/bioinf/voluntary/bioinfcomp19/final_1");
fileName <- "task1.in";

nTests <- scan(fileName, nlines=1);
inData <- scan(fileName, nlines=1, skip=1);
nVerts <- inData[1];
nSegs <- inData[2];
nAdj <- inData[3];

library(igraph);
library(ggraph);
library(tidygraph);

library(Matrix);

segEdges <- matrix(scan(fileName, skip=2, nlines=nSegs), ncol=nSegs);
adjEdges <- matrix(scan(fileName, skip=2+nSegs, nlines=nAdj), ncol=nAdj);
combEdges <- cbind(segEdges, adjEdges);

outLines <- NULL;
## Set up variables
adjMatSeg <- sparseMatrix(i=c(segEdges[1,], segEdges[2,]),
                          j=c(segEdges[2,], segEdges[1,]),
                          x=rep(segEdges[3,], 2));
adjMatAdj <- sparseMatrix(i=c(adjEdges[1,], adjEdges[2,]),
                          j=c(adjEdges[2,], adjEdges[1,]),
                          x=rep(adjEdges[3,], 2));
adjMatIDs <- sparseMatrix(i=c(combEdges[1,], combEdges[2,]), 
                          j=c(combEdges[2,], combEdges[1,]), 
                          x=rep(seq_len(nSegs+nAdj), 2));
outData <- NULL;
## Start the path traversal process by finding an undiscovered segment
whichSegPoss <- which(adjMatSeg > 0, arr.ind = TRUE);
while(nrow(whichSegPoss) > 0){
  cp <- head(whichSegPoss, 1);
  adjMatSeg[cp[1],cp[2]] <- 
    adjMatSeg[cp[1],cp[2]]-1;
  adjMatSeg[cp[2],cp[1]] <- 
    adjMatSeg[cp[2],cp[1]]-1;
  vertexPath <- c(cp);
  edgePath <- adjMatIDs[cp[1],cp[2]];
  ## Find a matching adjacency
  nextAdjVertex <- head(which(adjMatAdj[cp[2],] > 0),1);
  while(length(nextAdjVertex) > 0){
    cp <- cbind(row=cp[2], col=nextAdjVertex);
    vertexPath <- c(vertexPath, nextAdjVertex);
    edgePath <- c(edgePath, adjMatIDs[cp[1],cp[2]]); ## add adjacency to path
    adjMatAdj[cp[1],cp[2]] <- adjMatAdj[cp[1],cp[2]] - 1;
    adjMatAdj[cp[2],cp[1]] <- adjMatAdj[cp[2],cp[1]] - 1;
    ## Find a matching segment
    nextSegVertex <- head(which(adjMatSeg[cp[2],] > 0),1);
    if(length(nextSegVertex) == 0){
      nextAdjVertex <- NULL;
    } else {
      cp <- cbind(row=cp[2], col=nextSegVertex);
      vertexPath <- c(vertexPath, nextSegVertex);
      edgePath <- c(edgePath, adjMatIDs[cp[1],cp[2]]); ## add segment to path
      adjMatSeg[cp[1],cp[2]] <- adjMatSeg[cp[1],cp[2]] - 1;
      adjMatSeg[cp[2],cp[1]] <- adjMatSeg[cp[2],cp[1]] - 1;
      ## Back to the start of the loop
      nextAdjVertex <- head(which(adjMatAdj[cp[2],] > 0),1);
    }
  }
  ## Spit out result
  outData <- c(outData, paste(c(length(edgePath), edgePath), collapse=" "));
  whichSegPoss <- which(adjMatSeg > 0, arr.ind = TRUE);
  break;
}
outLines <- c(outLines, "YES", length(outData), outData);

cat(file="test1.out", outLines, sep="\n");