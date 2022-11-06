# Coll.analysis V 4.0
# Collostructional analysis: Computing the degree of association between words and words/constructions
# Copyright (C) 2022 Stefan Th. Gries (Latest changes in this version: 21 August 2022)

coll.analysis <- function() { # FUNCTION FOR THE FAMILY OF COLLOSTRUCTIONAL ANALYSES
   cat("\nColl.analysis 4.0 was written by Stefan Th. Gries (<https://www.stgries.info/>).\nIt computes all methods belonging to the family of collostructional analysis as developed by\nAnatol Stefanowitsch and Stefan Th. Gries. Thus, it can also be used to compute general\ncollocational strengths of word pairs or distinctive collocates.\n\nThis program is free software; you can redistribute it and/or modify it under the terms of the\nGNU General Public License as published by the Free Software Foundation; either version 2 of\nthe License, or (at your option) any later version.\n   Because the program is licensed free of charge, there is no warranty for the program, to the\nextent permitted by applicable law. Except when otherwise stated in writing the copyright holders\nand/or other parties provide the program 'as is' without warranty of any kind, either expressed\nor implied, including, but not limited to, the implied warranties of merchantability and fitness\nfor a particular purpose. The entire risk as to the quality and performance of the program is\nwith you. Should the program prove defective, you assume the cost of all necessary servicing,\nrepair or correction.\n   In no event unless required by applicable law or agreed to in writing will any copyright holder,\nor any other party who may modify and/or redistribute the program as permitted above, be liable\nto you for damages, including any general, special, incidental or consequential damages arising\nout of the use or inability to use the program (including but not limited to loss of data or\ndata being rendered inaccurate or losses sustained by you or third parties or a failure of the\nprogram to operate with any other programs), even if such holder or other party has been advised\nof the possibility of such damages.\n\nLatest changes in this version: 15 August 2022\n----------------------------\n\nYou should have received this program with a collection of example files and a readme file;\nI recommend that you have a look at them before you execute this program for the first time ...\n\n"); pause()
   cat("\nIf you use the program, you should make sure you read the following papers on my website, because they are behind many of the changes implemented in this version: 2019c, 2022c, d, in progress e\n\n"); pause()
   cat("\nIf you use the program, PLEASE QUOTE IT as follows:\nGries, Stefan Th. 2022. Coll.analysis 4.0. A script for R to compute perform collostructional analyses.\n\n"); pause()

   switch(menu(
      title="\nWhich kind of analysis do you want to perform?",
      choices=c("collocation/collexeme analysis (see <1.csv> for an example)",
                "(multiple) distinctive collocate/collexeme analysis (see <2a.csv>, <2b.csv>, <2c.csv> for examples)",
                "co-varying collexeme analysis (see <3.csv> for an example)")),
      collexemes(),
      dist.collexemes(),
      covar.collexemes())
} # END OF FUNCTION FOR THE FAMILY OF COLLOSTRUCTIONAL ANALYSES



collexemes <- function() { # FUNCTION FOR COLLEXEME ANALYSIS
   options(warn=-1)
   cat("\nC o l l o c a t i o n a l / c o l l e x e m e    a n a l y s i s   . . .\n\nThis kind of analysis computes the degree of attraction and repulsion between\none word or construction and many other words using a user-defined statistic;\nall these statistics are based on 2-by-2 tables, and attraction and repulsion\nare indicated in a separate column in the output.\n")

   # input
   cat("\nWhat is the word/construction you're investigating (without spaces)?\n")
      construction.name <- scan(nmax=1, what=character(), quiet=TRUE)
   cat("\nEnter the size of the corpus (without digit grouping symbols)!\n")
      corpus <- scan(nmax=1, quiet=TRUE)
   cat("\nDo want the results of (two-tailed!) Fisher-Yates exact tests ('yes' or 'no')?\n(I recommend not using this. The present script improves on the old one(s) by being able to avoid the (-)Inf problem of the old one, but (i) this means it requires the package Rmpfr for that, (ii) because that can take a very long time, this script uses multiple threads and, thus, recommends the package doParallel, (iii) even with parallelization and one other 'trick', this can STILL take a long time, (iv) the resulting -log10(pFYE) values can be 0.95 or more correlated with the fast-to-compute LLR/G^2-values or the superfast-to-compute Pearson residuals so why bother with the overkill?, and (v) as I wrote in multiple pubs by now, maybe keeping frequency and association more separate is a better approach (for all but the most exploratory studies!).\n")
      fye.mpfr <- scan(nmax=1, what=character(), quiet=TRUE)
   cat("\nLoad the tab-delimited input file:\n"); pause()
      input.data <- read.table(file.choose(), header=TRUE, sep="\t", quote="", comment.char="")

   # computation
   input.matrix <- matrix(c(input.data[,3], input.data[,2]-input.data[,3]), ncol=2)
      input.matrix <- rbind(input.matrix, c(0, corpus-sum(input.matrix)))
      rownames(input.matrix) <- c(input.data[,1], "___OTHER___")
      colnames(input.matrix) <- c(construction.name, "___OTHER___")
   pearson.residuals <- chisq.test(input.matrix, correct=FALSE)$residuals[,1]

   all.2.by.2.matrices <- apply(
      input.matrix, 1,
      \(af) { matrix(c(af, colSums(input.matrix)-af), byrow=TRUE, ncol=2) },
      simplify=FALSE)

   if (fye.mpfr=="yes") {
      FYE.values <- lapply(all.2.by.2.matrices,
         \(af) fisher.test.mpfr(af))
   }
   glms <- lapply(all.2.by.2.matrices,
      \(af) glm(rbind(af[1,], af[2,]) ~ c(1:2), family=binomial))
   log.odds.ratios <- -sapply(glms, coefficients)[2,]
   log.likelihood.values <- sapply(glms, "[[", "null.deviance")
   mi.scores <- sapply(all.2.by.2.matrices,
      \(af) log2(af[1,1] / chisq.test(af, correct=FALSE)$exp[1,1]))
   delta.p.constr.cues.word <- sapply(all.2.by.2.matrices,
      \(af) af[1,1]/sum(af[,1]) - af[1,2]/sum(af[,2]))
   delta.p.word.cues.constr <- sapply(all.2.by.2.matrices,
      \(af) af[1,1]/sum(af[1,]) - af[2,1]/sum(af[2,]))
   relations <- sapply(pearson.residuals,
      \(af) switch(sign(af)+2, "repulsion", "chance", "attraction"))

   # output
   output.table <- data.frame(WORD=rownames(input.matrix), CONSTRUCTION=input.matrix[,1], OTHER=input.matrix[,2], row.names=NULL)
   output.table <- data.frame(output.table, RELATION=relations, LLR=log.likelihood.values, PEARSONRESID=pearson.residuals,
                              LOGODDSRATIO=log.odds.ratios, MI=mi.scores,
                              DELTAPC2W=delta.p.constr.cues.word, DELTAPW2C=delta.p.word.cues.constr, row.names=NULL)
   if (fye.mpfr=="yes") {
      output.table <- data.frame(output.table,
         # FYE=sapply(FYE.values, formatMpfr, digits=12))
         FYE=sapply(sapply(FYE.values, \(af) -log10(af)), asNumeric))
   }
   output.table <- output.table[-nrow(output.table),]
   colnames(output.table)[2] <- construction.name
   output.table <- output.table[order(output.table$RELATION, -output.table$LOGODDSRATIO),]
   write.table(output.table, file=save.date <- gsub(":", "-", paste0(Sys.time(), ".csv")), sep="\t", row.names=FALSE, col.names=TRUE, quote=FALSE)
   cat("\n\nThe results are in the file called ", save.date, ".\n")

   plot(log2(output.table[,2]), output.table$LOGODDSRATIO, type="n",
      xlab="Logged co-occurrence frequency", ylab="Association (log odds ratio)")
      grid(); abline(h=0, lty=2); abline(v=0, lty=2)
      text(log2(output.table[,2]), output.table$LOGODDSRATIO, output.table$WORD, font=3)
   options(warn=0)
} # END OF FUNCTION FOR COLLEXEME ANALYSIS



dist.collexemes<-function(precbitsexponent=precbitsexponent) { # FUNCTION FOR DISTINCTIVE COLLEXEME ANALYSIS
   options(warn=-1)
   # introduction and first input
   cat("\nD i s t i n c t i v e   c o l l o c a t e / c o l l e x e m e   a n a l y s i s   . . .\n\nThis kind of analysis compares 2+ words or constructions with respect to n words they co-occur with differently frequently.\nYou must first enter whether you have two distinctive categories (e.g., when you look at English ditransitive\nvs. prep. dative) or more (e.g., when you compare English active vs. be-passive vs. get-passive)?\n")
   dists <- menu(title="How many distinctive categories do you have?",
                 choices=c(" 2 alternatives", " 3+ alternatives"))

   if (dists==1) {
      # introduction
      cat("\nColl.analysis accepts two kinds of input for such an analysis of distinctive collexemes:\nOn the one hand, you can use as input a file with a table of all tokens. That is, the first column\ncontains for each co-occurrence item the code for one of the two words/constructions W1/C1 and\nW2/C2 you want to investigate; the second column contains the word co-occurring with W1/C1 and W2/C2\nas listed in the first column.\n\nW/C\tColl_Word\nA\tX\nB\tY\n...\t...\n\nOn the other hand, if you have already down more work, you can also use a text file\nwith the following kind of table (with informative column names!), where the columns 2 and 3\ncontain the co-occurrence frequencies of each word listed in column 1 with/in W/C1 and W/C2.\n\nColl_Word\tFreq_CollWord_&_W/C1\tFreq_CollWord_&_W/C2\nA\t\t...\t\t\t...\nB\t\t...\t\t\t...\n...\t\t...\t\t\t...\n\nWhichever input format you choose, your file must not have decimal points/separators and ideally has no spaces (for the latter, use '_' instead)!\nAlso, don't forget that R's treatment of alphanumeric characters is case-sensitive!\n\n")
      input.dc <- menu(title="Which input format do you want to use?",
                       choices=c("Raw list of all tokens", "Frequency table"))

      cat("\nDo want the results of (two-tailed!) Fisher-Yates exact tests ('yes' or 'no')?\n(I recommend not using this. The present script improves on the old one(s) by being able to avoid the (-)Inf problem of the old one, but (i) this means it requires the package Rmpfr for that, (ii) because that can take a very long time, this script uses multiple threads and, thus, recommends the package doParallel, (iii) even with parallelization and one other 'trick', this can STILL take a long time, (iv) the resulting -log10(pFYE) values can be 0.95 or more correlated with the fast-to-compute LLR/G^2-values or the superfast-to-compute Pearson residuals so why bother with the overkill?, and (v) as I wrote in multiple pubs by now, maybe keeping frequency and association more separate is a better approach (for all but the most exploratory studies!).\n")
      fye.mpfr <- scan(nmax=1, what=character(), quiet=TRUE)

      cat("\nLoad the tab-delimited input file:\n"); pause()
         input.matrix <- read.table(file.choose(), header=TRUE, sep="\t", quote="", comment.char="")

      if (input.dc==1) { # DCA
         interim <- t(table(input.matrix))
         input.matrix <- data.frame(
            as.vector(interim[,1]),
            as.vector(interim[,2]), row.names=rownames(interim))
         colnames(input.matrix) <- colnames(interim)
      } else {
         temp <- colnames(input.matrix)[2:3]
         input.matrix <- data.frame(
            as.vector(input.matrix[,2]),
            as.vector(input.matrix[,3]), row.names=input.matrix[,1])
         colnames(input.matrix) <- temp
      }
      construction1.name <- colnames(input.matrix)[1]; construction2.name <- colnames(input.matrix)[2]

      # computation
      options(warn=-1)
      pearson.residuals <- chisq.test(input.matrix, correct=FALSE)$residuals[,1]

      all.2.by.2.matrices <- apply(
         input.matrix, 1,
         \(af) { matrix(c(af, colSums(input.matrix)-af), byrow=TRUE, ncol=2) },
         simplify=FALSE)

      if (fye.mpfr=="yes") {
         FYE.values <- lapply(all.2.by.2.matrices,
            \(af) fisher.test.mpfr(af))
      }
      glms <- lapply(all.2.by.2.matrices,
         \(af) glm(rbind(af[1,], af[2,]) ~ c(1:2), family=binomial))
      log.odds.ratios <- -sapply(glms, coefficients)[2,]
      log.likelihood.values <- sapply(glms, "[[", "null.deviance")
      mi.scores <- sapply(all.2.by.2.matrices,
         \(af) log2(af[1,1] / chisq.test(af, correct=FALSE)$exp[1,1]))
      delta.p.constr.cues.word <- sapply(all.2.by.2.matrices,
         \(af) af[1,1]/sum(af[,1]) - af[1,2]/sum(af[,2]))
      delta.p.word.cues.constr <- sapply(all.2.by.2.matrices,
         \(af) af[1,1]/sum(af[1,]) - af[2,1]/sum(af[2,]))
      relations <- sapply(pearson.residuals,
         \(af) switch(sign(af)+2, construction2.name, "chance", construction1.name))

      # output
      output.table <- data.frame(WORD=rownames(input.matrix), CONSTRUCTION1=input.matrix[,1], CONSTRUCTION2=input.matrix[,2], row.names=NULL)
      output.table <- data.frame(output.table, PREFERENCE=relations, LLR=log.likelihood.values, PEARSONRESID=pearson.residuals,
                                 LOGODDSRATIO=log.odds.ratios, MI=mi.scores,
                                 DELTAPC2W=delta.p.constr.cues.word, DELTAPW2C=delta.p.word.cues.constr, row.names=NULL)
      if (fye.mpfr=="yes") {
         output.table <- data.frame(output.table,
            # FYE=sapply(FYE.values, formatMpfr, digits=12))
            FYE=sapply(sapply(FYE.values, \(af) -log10(af)), asNumeric))
      }
      colnames(output.table)[2:3] <- c(construction1.name, construction2.name)
      output.table <- output.table[order(output.table$PREFERENCE, -output.table$LOGODDSRATIO),]
      write.table(output.table, file=save.date <- gsub(":", "-", paste0(Sys.time(), ".csv")), sep="\t", row.names=FALSE, col.names=TRUE, quote=FALSE)
      cat("\n\nThe results are in the file called ", save.date, ".\n")

      plot(log2(output.table[,2]+output.table[,3]), output.table$LOGODDSRATIO, type="n",
         xlab="Logged co-occurrence frequency", ylab="Association (log odds ratio)")
         grid(); abline(h=0, lty=2); abline(v=0, lty=2)
         text(log2(output.table[,2]+output.table[,3]), output.table$LOGODDSRATIO, output.table$WORD, font=3)

   } else { # MDCA

      # introduction
      cat("\nFor a multiple distinctive collexeme analysis, this version of coll.analysis expects as input\na file with a table of all tokens and returns collexeme strengths based on Pearson residuals. That is, the first column contains for\neach co-occurrence item the code for one of the X words/constructions W/C\nyou want to investigate; the second column contains the word co-occurring with W/C\nas listed in the first column.\n\nW/C\tColl_Word\nA\tX\nB\tY\nC\tZ\n...\t...\n\nYour file ideally has no spaces (use '_' instead) and don't forget that R's treatment of alphanumeric characters\nis case-sensitive!\n\nChoose the text file with the input data!\t"); pause()

      # input
      cat("\nLoad the tab-delimited input file:\n"); pause()
         input.matrix <- read.table(file.choose(), header=TRUE, sep="\t", quote="", comment.char="")
         names(input.matrix) <- c("W_C", "Coll_Word")
         input.matrix <- table(input.matrix$Coll_Word, input.matrix$W_C)

      # computation
      pearson.residuals <- as.data.frame.matrix(chisq.test(input.matrix, correct=FALSE)$residuals)
      output.table <- data.frame(COLLOCATE=rownames(pearson.residuals), pearson.residuals,
         SUMABSDEV=apply(pearson.residuals, 1, \(af) sum(abs(af))),
         LARGESTPREF=colnames(pearson.residuals)[apply(pearson.residuals, 1, \(af) which.max(af))])
      output.table <- output.table[order(output.table$LARGESTPREF, -output.table$SUMABSDEV),]
      write.table(output.table, file=save.date <- gsub(":", "-", paste0(Sys.time(), ".csv")), sep="\t", row.names=FALSE, col.names=TRUE, quote=FALSE)
      cat("\n\nThe results are in the file called ", save.date, ".\n")
   }
   options(warn=0)
} # END OF FUNCTION FOR DISTINCTIVE COLLEXEME ANALYSIS



covar.collexemes<-function(precbitsexponent=precbitsexponent) { # FUNCTION FOR CO-VARYING COLLEXEME ANALYSIS
   options(warn=-1)
   # introduction
   cat("\nC o v a r y i n g   c o l l e x e m e   a n a l y s i s   . . .\n\n\nThis kind of analysis investigated dependencies within two slots of a single construction.\nThis script so far only implements the so-called item-based analysis since comparative studies\n have shown that the system-based correction may require many days computational time with only minor differences in the results (cf. Stefanowitsch and Gries 2005). However, somewhere down the road I may find \ntime to work on an implementation of this technique so that arbitrarily many additional variables\n(e.g. register, corpora etc.) can be included.\n\nColl.analysis 3.2a requires as input for the item-based co-varying collexeme analysis:\na file with a table of all token instances of the construction C with\nthe two words W1 and W2 occurring in the slots of each instance of C.That is, you need the following kind of input file (with column names!)),\nwhere the number of rows corresponds to the number of construction tokens you have.\n\nWord_Slot1\tWord_Slot2\nA\t\tX\nB\t\tX\n...\t...\n\nYour file must not have decimal points/separators and ideally has no spaces (for the latter, use '_' instead)!\nAlso, don't forget that R's treatment of alphanumeric characters is case-sensitive!\n\n")

   # input
   cat("\nLoad the tab-delimited input file:\n"); pause()
   input.data <- read.table( file.choose(), header=TRUE, sep="\t", quote="", comment.char="")
   output.table <- as.data.frame(table(input.data))
      output.table$FREQOFSLOT1 <- table(input.data$WORD_SLOT1)[output.table[,1]]
      output.table$FREQOFSLOT2 <- table(input.data$WORD_SLOT2)[output.table[,2]]
   construction.freq <- sum(output.table[,3])

   # computation
   all.2.by.2.matrices <- apply(
      output.table[,3:5], 1,
      \(af) { matrix(c(af[1], af[2]-af[1], af[3]-af[1], construction.freq-sum(af)), byrow=TRUE, ncol=2) },
      simplify=FALSE)

   glms <- lapply(all.2.by.2.matrices,
      \(af) glm(rbind(af[1,], af[2,]) ~ c(1:2), family=binomial))
   log.odds.ratios <- -sapply(glms, coefficients)[2,]
   log.likelihood.values <- sapply(glms, "[[", "null.deviance") * sign(log.odds.ratios)
   mi.scores <- sapply(all.2.by.2.matrices,
      \(af) log2(af[1,1] / chisq.test(af, correct=FALSE)$exp[1,1]))
   delta.p.slot1.cues.slot2 <- sapply(all.2.by.2.matrices,
      \(af) af[1,1]/sum(af[1,]) - af[2,1]/sum(af[2,]))
   delta.p.slot2.cues.slot1 <- sapply(all.2.by.2.matrices,
      \(af) af[1,1]/sum(af[,1]) - af[1,2]/sum(af[,2]))
   relations <- sapply(log.odds.ratios,
      \(af) switch(sign(af)+2, "repulsion", "chance", "attraction"))

   # output
   output.table <- data.frame(output.table, RELATION=relations, LLR=log.likelihood.values,
                              LOGODDSRATIO=log.odds.ratios, MI=mi.scores,
                              DELTAP1TO2=delta.p.slot1.cues.slot2, DELTAP2TO1=delta.p.slot2.cues.slot1, row.names=NULL)
      output.table <- output.table[order(output.table$RELATION, -output.table$LOGODDSRATIO),]
   write.table(output.table, file=save.date <- gsub(":", "-", paste0(Sys.time(), ".csv")), sep="\t", row.names=FALSE, col.names=TRUE, quote=FALSE)
   cat("\n\nThe results are in the file called ", save.date, ".\n")
   options(warn=0)
} # END OF FUNCTION FOR CO-VARYING COLLEXEME ANALYSIS

pause <- function() { cat("Press <Enter> to continue ... "); readline(); invisible() }

fisher.test.mpfr <- function (some.2by2.matrix, precBits=1000) {
   options(warn=-1)
   stopifnot("Rmpfr" %in% rownames(installed.packages()))
   library(Rmpfr, warn.conflicts=FALSE, quietly=TRUE, verbose=FALSE)

   # define helper functions
   binomcoeff.mpfr <- function (k.successes, n.trials, precBits=1000) {
      numerator <-   factorial(mpfr(n.trials,             precBits=precBits))
      denominator <- factorial(mpfr(k.successes,          precBits=precBits)) *
                     factorial(mpfr(n.trials-k.successes, precBits=precBits))
      return(numerator/denominator)
   }
   hypergeom.mpfr <- function (a, a.and.b, a.and.c, n, precBits=1000) {
      numerator.1 <- binomcoeff.mpfr(a        , a.and.b        , precBits=precBits)
      numerator.2 <- binomcoeff.mpfr(a.and.c-a, n      -a.and.b, precBits=precBits)
      denominator <- binomcoeff.mpfr(a.and.c  , n              , precBits=precBits)
      return((numerator.1*numerator.2)/denominator)
   }

   range.of.as <- 0:min(rowSums(some.2by2.matrix)[1], colSums(some.2by2.matrix)[1])

   if ("doParallel" %in% rownames(installed.packages())) { # if doParallel is installed
      library(doParallel); registerDoParallel(cl <- makePSOCKcluster(detectCores()-1, outfile=""))
      all.results <- # make all.results the result of doing w/ a cluster
      foreach(i=seq(range.of.as), .packages=c("Rmpfr"), .verbose=FALSE) %dopar% { cat(".")
         curr.matrix <- matrix(c(
            a <- range.of.as[i],                              # a
            b <- rowSums(some.2by2.matrix)[1]-range.of.as[i], # b
            c <- colSums(some.2by2.matrix)[1]-range.of.as[i], # c
            sum(some.2by2.matrix)-sum(a,b,c)),                # d
            byrow=TRUE, ncol=2)
         # check if the 'normal' way works, because it's so much faster
         hypergeom.quick <- dhyper(curr.matrix[1,1], rowSums(curr.matrix)[1], rowSums(curr.matrix)[2], colSums(curr.matrix)[1])
         if (hypergeom.quick==0) {
            output <- c(
               "chi.sq.values"=chisq.test(curr.matrix, correct=FALSE)$statistic,
               "pointwise.hypgeom.p"=hypergeom.mpfr(
               curr.matrix[1,1],        # a
               rowSums(curr.matrix)[1], # a+b
               colSums(curr.matrix)[1], # a+c
               sum(curr.matrix),        # n
               precBits=precBits))      # precision
            return(output)
         } else {
            output <- c(
               "chi.sq.values"=chisq.test(curr.matrix, correct=FALSE)$statistic,
               "pointwise.hypgeom.p"=mpfr(hypergeom.quick, precBits=precBits))
            return(output)
         }
      }; stopCluster(cl)
      chi.sq.values <- sapply(all.results, "[[", 1)
      pointwise.hypgeom.p <- sapply(all.results, "[[", 2)
         names(pointwise.hypgeom.p) <- names(chi.sq.values) <- range.of.as
   } else { # if doParallel is not installed
      # set up collectors for the relevant scenario
      pointwise.hypgeom.p <- vector(mode="list", length=length(range.of.as))
      chi.sq.values <- rep(NA, length(range.of.as))
         names(pointwise.hypgeom.p) <- names(chi.sq.values) <- range.of.as
      for (i in seq(range.of.as)) {
         curr.matrix <- matrix(c(
            a <- range.of.as[i],                              # a
            b <- rowSums(some.2by2.matrix)[1]-range.of.as[i], # b
            c <- colSums(some.2by2.matrix)[1]-range.of.as[i], # c
            sum(some.2by2.matrix)-sum(a,b,c)),                # d
            byrow=TRUE, ncol=2)
         chi.sq.values[i] <- chisq.test(curr.matrix, correct=FALSE)$statistic
         # check if the 'normal' way works, because it's so much faster
         hypergeom.quick <- dhyper(curr.matrix[1,1], rowSums(curr.matrix)[1], rowSums(curr.matrix)[2], colSums(curr.matrix)[1])
         if (hypergeom.quick==0) {
            pointwise.hypgeom.p[i] <- hypergeom.mpfr(
               curr.matrix[1,1],        # a
               rowSums(curr.matrix)[1], # a+b
               colSums(curr.matrix)[1], # a+c
               sum(curr.matrix),        # n
               precBits=precBits)       # precision
         } else {
            pointwise.hypgeom.p[i] <- mpfr(hypergeom.quick, precBits=precBits)
         }
      }
   }

   # compute result
   pickito <- which(chi.sq.values >= chi.sq.values[which(names(pointwise.hypgeom.p) == some.2by2.matrix[1,1])])
   return(sum(mpfr(pointwise.hypgeom.p)[pickito]))
   options(warn=0)
}
