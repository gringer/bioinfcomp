setwd("/bioinf/voluntary/bioinfcomp19/final_3");
fileName <- "tests/4.txt";

library(Biostrings);

nTests <- scan(fileName, nlines=1);

sM <- nucleotideSubstitutionMatrix(match = 0, mismatch=-1, baseOnly = FALSE);

processed <- 0;
linePos <- 1;
outOrders <- NULL;
## Initial guess: look at pairwise differences within each population
while(processed < nTests){
  cat(processed,"");
  nGens <- scan(fileName, nlines=1, skip=linePos, quiet=TRUE);
  pop <- matrix(ncol = nGens,
    scan(fileName, nlines=nGens, skip=linePos+1, what=character(), quiet=TRUE));
  meanDist <- sapply(seq_len(nGens), function(x){
    p1 <- DNAStringSet(pop[,x]);
    ##mean(stringDist(p1, method="hamming"));
    c1 <- consensusString(p1, threshold=0.25);
    mean(-pairwiseAlignment(p1, DNAString(c1), gapOpening=50, 
                            scoreOnly=TRUE, substitutionMatrix=sM));
  });
  outOrders <- c(outOrders, paste(order(meanDist)-1, collapse=" "));
  processed <- processed + 1;
  linePos <- linePos + 1 + nGens;
}
cat("\n");

mdOrder <- order(meanDist);

write.csv(t(pop[,mdOrder]), "out_ordered.txt");

cat(file=sprintf("output/%s", basename(fileName)), 
    outOrders, sep="\n");

## For problem #4; use condorcet voting to determine ideal rank
## [assuming randomly sampled votes per generation]
orderVotes <- sapply(seq_len(nrow(pop)),function(x){
  hclust(stringDist(pop[x,], method = "hamming"))$order});


orderPlacement <- sapply(seq_len(nrow(orderVotes)), function(x){
  which(orderVotes == x, arr.ind=TRUE)[,1]
  });

image(orderPlacement[order(meanDist), order(meanDist)]);

opOrder <- rev(order(apply(orderPlacement,2,mean)));
image(orderPlacement[opOrder, opOrder]);

pc <- PairCount(orderVotes);
hco <- hclust(dist(pc))$order;

library(viridis);
pdf("out.pdf", width=20, height=20);
image(log(pc[hco,hco]), col=viridis(300));
dummy <- dev.off();

cat(file=sprintf("output/%s", basename(fileName)), 
    paste(order(-apply(orderPlacement,2,mean))-1, collapse=" "), sep="\n");

cat(file=sprintf("output/%s", basename(fileName)), 
    paste(hco-1, collapse=" "), sep="\n");