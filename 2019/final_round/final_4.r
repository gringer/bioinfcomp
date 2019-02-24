setwd("/bioinf/voluntary/bioinfcomp19/final_4");
fileName <- "tests/input_2_backtran.fa";

library(Biostrings);

DNASeqs <- readDNAStringSet(fileName);

#### Test 1

## remove duplicates
DNASeqs.rmdup <- DNASeqs[1];
for(x in seq_along(DNASeqs)){
  if(x %% 100 == 0){
    print(x);
  }
  xs <- DNASeqs[[x]];
  if(!any(unlist(startIndex(vmatchPattern(xs, DNASeqs.rmdup)))>1)){
    DNASeqs.rmdup[[length(DNASeqs.rmdup)+1]] <- xs;
  }
}

DNASeqs.rmdup <- 
  DNASeqs[!sapply(DNASeqs, function(x){
    any(unlist(startIndex(vmatchPattern(x, DNASeqs)))>1)})];
DNAseq.collapsed <- paste(DNASeqs.rmdup,collapse="");

matchStarts <-
  sapply(DNASeqs, function(x){start(matchPattern(x, DNAseq.collapsed))});

cat(DNAseq.collapsed, 
    paste(matchStarts,"+"), file = "output/out.1.txt", sep="\n");

#### End Test 1

#### Test 2

setwd("/bioinf/voluntary/bioinfcomp19/final_4");
fileName <- "tests/input_2_backtran.fa";

library(Biostrings);

DNASeqs <- readDNAStringSet(fileName);
superGenome <- paste(DNASeqs, collapse="");
matchCounts <-
  sapply(DNASeqs, function(x){length(matchPattern(x, superGenome))});

DNASeqs.rmdup <- DNASeqs[matchCounts == 1];
miniGenome <- paste(DNASeqs.rmdup, collapse="");
matchStarts <-
  sapply(DNASeqs, function(x){min(start(matchPattern(x, miniGenome)))});

cat(miniGenome, 
    paste(matchStarts,"+"), file = "output/out.2.txt", sep="\n");

#### Tests 3-6 (stupid / forward-only, no ambiguities)

setwd("/bioinf/voluntary/bioinfcomp19/final_4");
fileName <- "tests/input_6_backtran.fa";

library(Biostrings);

DNASeqs <- readDNAStringSet(fileName);
superGenome <- paste(DNASeqs, collapse="");
matchCounts <-
  sapply(DNASeqs, function(x){length(matchPattern(x, superGenome))});

DNASeqs.rmdup <- DNASeqs[matchCounts == 1];
miniGenome <- paste(DNASeqs.rmdup, collapse="");
matchStarts <-
  sapply(DNASeqs, function(x){min(start(matchPattern(x, miniGenome)))});

cat(miniGenome, 
    paste(matchStarts,"+"), file = "output/out.6.txt", sep="\n");

#### Tests 3-6 (forward-only with ambiguities)... too difficult 1h away from end

setwd("/bioinf/voluntary/bioinfcomp19/final_4");
fileName <- "tests/input_3_backtranambig.fa";

library(Biostrings);

DNASeqs <- readDNAStringSet(fileName);
superGenome <- DNAString(paste(DNASeqs, collapse=""));
matchCounts <-
  sapply(DNASeqs, function(x){
    length(matchPattern(x, superGenome, fixed=FALSE)) +
      length(matchPattern(reverseComplement(x), superGenome, fixed=FALSE))
  });

DNASeqs.rmdup <- DNASeqs[matchCounts == 1];
miniGenome <- DNAString(paste(DNASeqs.rmdup, collapse=""));
matchStarts <-
  sapply(DNASeqs, function(x){
    min(start(matchPattern(x, miniGenome, fixed=FALSE)),
        start(matchPattern(reverseComplement(x), miniGenome)))});
matchDirs <-
  sapply(DNASeqs, function(x){
    ifelse(start(matchPattern(x, miniGenome, fixed=FALSE)),
        start(matchPattern(reverseComplement(x), miniGenome)))});

cat(miniGenome, 
    paste(matchStarts,"+"), file = "output/out.6.txt", sep="\n");
