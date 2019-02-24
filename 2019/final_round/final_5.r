setwd("/bioinf/voluntary/bioinfcomp19/final_5");
fileName <- "tests/1.in";

inData <- scan(fileName, nlines=1);
nchrs <- inData[1];
nsite <- inData[2];
siteDesc <- matrix(scan(fileName, skip=1, nlines=nsite), nrow=2, ncol=nsite);
ninds <- scan(fileName, nlines=1, skip=1+nsite);

minds <- array(dim=c(2, nsite, ninds));
finds <- array(dim=c(2, nsite, ninds));
processed <- 0;
linePos <- 2+nsite;
while(processed < ninds){
  processed <- processed + 1;
  minds[,,processed] <- (scan(fileName, nlines=nsite, sep="/", skip=linePos,
                             quiet=TRUE)-0.5) * 2;
  linePos <- linePos + nsite + 1;
}
processed <- 0;
while(processed < ninds){
  processed <- processed + 1;
  finds[,,processed] <- (scan(fileName, nlines=nsite, sep="/", skip=linePos,
                             quiet=TRUE)-0.5) * 2;
  linePos <- linePos + nsite + 1;
}

#sapply(seq_len(ninds), function(x){
#  sapply(seq_len(ninds), function(y){
#    sum(abs(minds[,,x] + finds[,,y]));
#  })
#});

cat(file=sprintf("output/%s", sub("\\.in$",".out",basename(fileName))), 
    sample(ninds), sep="\n");