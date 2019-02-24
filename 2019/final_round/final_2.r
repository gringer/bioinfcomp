setwd("/bioinf/voluntary/bioinfcomp19/final_2");
fileName <- "tests/input_1.txt";

nStates <- scan(fileName, nlines=1);
markLines <- sapply(strsplit(scan(fileName, skip=1, what=character()),""),
                    as.numeric);

plot(hclust(dist(t(markLines))));


markText <- apply(markLines,1,paste,collapse="");


myStates <- rowSums(markLines);

00000 - 0
00010 - 1 [00110, 00100]
10010 - 2
10100 - 3

outStates <- rep(0, length(myStates));
outStates[markText %in% c("01000")] <- 0;
outStates[markText %in% c("01010","00100","00010",
                          "00110","00111", "01101","00011")] <- 1;
outStates[markText %in% c("10000","10010","11011", "11000",
                          "10001","10100","11100", "11010")] <- 2;
outStates[markText %in% c("01110","01100")] <- 3;

##myStates[myStates == 4] <- 2; -- only needed if there are more counts than states

cat(file=sprintf("output/%s", sub("^input","output",basename(fileName))), 
    outStates, sep="\n");

cat(file=sprintf("output/%s", sub("^input","output",basename(fileName))), 
    myStates, sep="\n");

cat(file=sprintf("output/%s", sub("^input","output",basename(fileName))), 
    sample(nStates, size = nrow(markLines), replace=TRUE)-1, sep="\n");