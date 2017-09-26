.countGuitarDensity <- function(peak,GuitarCoords,maximalFeatureAmbiguity) {
  
  # count overlaps
  n <- countOverlaps(GuitarCoords,peak)
  
  # normalize by overlaps
  tx <- as.character(mcols(GuitarCoords)[[1]])
  txlist <- split(x=GuitarCoords,f=tx)
  peak_overlap <- countOverlaps(peak,txlist)
  
  # split features based on ambiguity
  z <- matrix(rep(n,maximalFeatureAmbiguity),ncol=maximalFeatureAmbiguity)
  for (i in 1:maximalFeatureAmbiguity) {
    temp <- countOverlaps(GuitarCoords,peak[peak_overlap==i])
    z[,i] = temp/i
  }
  n <- rowSums(z)
  q <- data.frame(mcols(GuitarCoords),count=as.numeric(n))
  return(q)
  
}
.combineListOfDataFrame <- function(t){
  if (length(t)==1) {
    return(t[[1]]) 
  } else {
    txt <- "rbind(t[[1]]"
    for (i in 2:length(t)) {
      txt <- paste(txt,",t[[",i,"]]",sep="")
    }
    txt <- paste(txt,")",sep="")
    c <- parse(text=txt)
    newframe <- eval(c)
    return(newframe)
  } 
}

.countDensityCI <- function(ct1, adjust, CI_interval){
  
  den_ct1 <- data.frame()
  
  for (l in 1:length(levels(ct1$Feature))){

    id <- which(ct1$Feature == levels(ct1$Feature)[l])
    ct_tmp <- ct1[id,]

    xx <- ct_tmp$pos
    ww <- ct_tmp$weight
    fit1 <- density(xx, adjust = adjust, weights = ww) 

    fit2 <- replicate(10000, { sid <- 1:length(xx);idd <- sample(sid, replace=TRUE);x <- xx[idd];wwx <- ww[idd]; 
    density(x, adjust = adjust, weights = wwx)$y}) # , from=min(fit1$x), to=max(fit1$x)
    fit3 <- apply(fit2, 1, quantile, CI_interval ) 
    tmp <-  data.frame( x = fit1$x, y = fit1$y *5, Feature = rep(levels(ct1$Feature)[l], times = length(fit1$y)), confidence_down = fit3[1,] *5, confidence_up = fit3[2,] *5)
    den_ct1 <- rbind(tmp, den_ct1)
  }
 
   return(den_ct1)
}


Guitar2Plot <- function(gfeatures,
                       GuitarCoordsFromTxDb=NA,
                       txdb=NA,
                       genome=NA,
                       noBins=10,
                       transcript='both',
                       saveToPDFprefix=NA,
                       returnCount=FALSE,
                       includeNeighborDNA=FALSE,
                       maximalFeatureAmbiguity=5,
                       rescaleComponent=TRUE,
                       CI_interval=c(0.025,0.975),
                       adjust=1){

  # make sure the Guitar coordinates are available
  suppressWarnings(
    if (is.na(GuitarCoordsFromTxDb)&is.na(txdb)&is.na(genome)) {
      stop("Must provide at least one of the three: GuitarCoords, txdb or genome")
    } 
  )
  
  if ( suppressWarnings(is.na(GuitarCoordsFromTxDb)) ) {
    if (suppressWarnings(is.na(txdb))) {
      print("Downloading Transcriptome Information from UCSC ...")
      txdb <- suppressMessages(makeTxDbFromUCSC(genome=genome))
      print("Making Guitar Coordinates ...")
      GuitarCoordsFromTxDb <- suppressMessages(makeGuitarCoordsFromTxDb(txdb))
      GuitarCoords <- GuitarCoordsFromTxDb
    } else {
      print("Making Guitar Coordinates from provided TranscriptDb Object ...")
      GuitarCoordsFromTxDb <- makeGuitarCoordsFromTxDb(txdb, noBins=noBins)
      GuitarCoords <- GuitarCoordsFromTxDb
    }
  } else {
    print("Using provided Guitar Coordinates")
    GuitarCoords <- GuitarCoordsFromTxDb
  }
  
  # Generate named List
  noGroup <- length(gfeatures)
  group_names <- names(gfeatures)
  m <- gfeatures
  if (is.null(group_names)) {
    group_names <- paste("item",1:noGroup)
  }
  
  print("resolving ambiguious features ...")
  for (i in 1:noGroup) {
    temp = .countGuitarDensity(
      gfeatures[[i]],
      GuitarCoords,
      maximalFeatureAmbiguity)
    temp = cbind(temp,Feature=group_names[i])
    m[[i]] =temp
  }
  ct=.combineListOfDataFrame(m)
  ct[[4]] <- as.character(ct[[4]])
  ct <- ct[ct$count>0,]  # remove redundant info
  
  # # organize output
    .makeFigure_nofill(ct,
                       GuitarCoordsFromTxDb,
                       includeNeighborDNA,
                       rescaleComponent,
                       CI_interval,
                       transcript,
                       saveToPDFprefix,
                       adjust)

  # return the result
  if (returnCount) {return(ct)}
}


.makeFigure_nofill <- function(ct,
                               GuitarCoordsFromTxDb,
                               includeNeighborDNA,
                               rescaleComponent, CI_interval, transcript,
                               saveToPDFprefix,adjust=adjust) {

  # extract information of mRNA and lncRNA
  ct$weight <- ct$count # as numeric
  ct1 <- ct[ct$category=="mRNA",] # mRNA
  ct2 <- ct[ct$category=="lncRNA",] # lncRNA
  
  # save(ct1,ct2, includeNeighborDNA,rescaleComponent,saveToPDFprefix,GuitarCoordsFromTxDb, file = "data.Rdata")
  d <- mcols(GuitarCoordsFromTxDb)
  
  # disable notes
  pos=Feature=weight=NULL
  
  # plot
  if (includeNeighborDNA) {
    if (rescaleComponent==FALSE) {
      
      # normalize feature
      featureSet <- as.character(unique(ct$Feature))
      for (i in 1:length(featureSet)) {
        id <- (ct1$Feature==featureSet[i])
        ct1$weight[id] <- ct1$weight[id]/sum(ct1$weight[id])
        
        id <- (ct2$Feature==featureSet[i])
        ct2$weight[id] <- ct2$weight[id]/sum(ct2$weight[id])
      }
      
      # adjust position mRNA
      pos_adjust <- match(ct1$comp,c("Front","UTR5","CDS","UTR3","Back"))-1
      ct1$pos <- ct1$pos + pos_adjust
      
      den_ct1 <- .countDensityCI(ct1, CI_interval)
      p1 <- 
        ggplot(den_ct1, aes(x=x, y = y, group = Feature))  +
        scale_x_continuous(minor_breaks = seq(0, 5, 1)) + 
        ggtitle("Distribution on mRNA") +
        theme(axis.ticks = element_blank(), axis.text.x = element_blank()) + 
        xlab("") + 
        ylab("Frequency") +
        geom_ribbon(aes(ymin =  rep(0, length(y)), ymax =  y, colour = factor(Feature), fill = factor(Feature)), alpha = 0.2 ) +
        geom_ribbon(aes(ymin =  confidence_down, ymax =  confidence_up, colour = factor(Feature), fill = factor(Feature)), alpha = 0.3 ) +
        geom_line(aes(y = y, colour = factor(Feature)), alpha = 0.8, size = 1) +
        annotate("text", x = 1.5, y = -0.8, label = "5'UTR") +
        annotate("text", x = 2.5, y = -0.8, label = "CDS") +
        annotate("text", x = 0.5, y = -0.8, label = "Promoter (1kb)") +
        annotate("text", x = 4.5, y = -0.8, label = "Tail (1kb)") +
        annotate("text", x = 3.5, y = -0.8, label = "3'UTR")  + 
        geom_vline(xintercept=1:4, linetype="dotted") + 
        annotate("rect", xmin = 1, xmax = 2, ymin = -0.22, ymax = -0.18, alpha = .99, colour = "black")+
        annotate("rect", xmin = 3, xmax = 4, ymin = -0.22, ymax = -0.18, alpha = .99, colour = "black")+
        annotate("rect", xmin = 2, xmax = 3, ymin = -0.26, ymax = -0.14, alpha = .2, colour = "black") +
        xlim(0,5) + 
        theme(legend.position="bottom")
      
      # adjust position lncRNA
      pos_adjust <- match(ct2$comp,c("Front","lncRNA","Back"))-1
      ct2$pos <- ct2$pos + pos_adjust

      den_ct2 <- .countDensityCI(ct2, adjust, CI_interval)
      p2 <- 
        ggplot(den_ct2, aes(x=x, y = y, group = Feature))  +
        ggtitle("Distribution on lncRNA")  +
        scale_x_continuous(minor_breaks = seq(6, 9, 1)) +
        theme(axis.ticks = element_blank(), axis.text.x = element_blank()) + 
        xlab("") + 
        ylab("Frequency") +
        geom_ribbon(aes(ymin =  rep(0, length(y)), ymax =  y, colour = factor(Feature), fill = factor(Feature)), alpha = 0.2 ) +
        geom_ribbon(aes(ymin =  confidence_down, ymax =  confidence_up, colour = factor(Feature), fill = factor(Feature)), alpha = 0.3 ) +
        geom_line(aes(y = y, colour = factor(Feature)), alpha = 0.8, size = 1) +
        annotate("text", x = 1.5, y = -0.6, label = "lncRNA") +
        annotate("text", x = 0.5, y = -0.6, label = "Promoter (1kb)") +
        annotate("text", x = 2.5, y = -0.6, label = "Tail (1kb)")  + 
        geom_vline(xintercept=1:2, linetype="dotted") + 
        annotate("rect", xmin = 1, xmax = 2, ymin = -0.22, ymax = -0.18, alpha = .99, colour = "black")+
        xlim(0,3) + 
        theme(legend.position="bottom")
    } else {
      
      if (rescaleComponent) {
        temp <- unique(d[,c(1,3,4,5)])
        id1 <- which(match(temp$category,"mRNA") >0 )
        temp <- temp[id1,]
        temp <- matrix(temp$interval,ncol=5)
        temp <- temp/rowSums(temp)
        temp <- colSums(temp)
        temp <- temp/sum(temp)
        weight <- temp
        names(weight) <- c("Promoter","5'UTR","CDS","3'UTR","Tail")
        w1 <- weight
        
        temp <- unique(d[,c(1,3,4,5)])
        id1 <- which(match(temp$category,"lncRNA") >0 )
        temp <- temp[id1,]
        temp <- matrix(temp$interval,ncol=3)
        temp <- temp/rowSums(temp)
        temp <- colSums(temp)
        temp <- temp/sum(temp)
        weight <- temp
        names(weight) <- c("Promoter","lncRNA","Tail")
        w2 <- weight
        
        x1 <- cumsum(w1)
        x2 <- cumsum(w2)
      }

      # adjust position
      id <- match(ct1$comp,c("Front","UTR5","CDS","UTR3","Back"))
      ct1$count <- ct1$count*w1[id]
      ct1$pos <- ct1$pos*w1[id] + c(0,x1)[id]
      
      id <- match(ct2$comp,c("Front","lncRNA","Back"))
      ct2$count <- ct2$count*w2[id]
      ct2$pos <- ct2$pos*w2[id] + c(0,x2)[id]

      # normalization by length of components in mRNA
      # calculate relative length of each components
      temp <- unique(d[,c(1,3,4,5)])
      id1 <- which(match(temp$comp,c("Front","Back")) >0 )
      temp <- temp[-id1,] # remove DNA
      id1 <- which(match(temp$category,"mRNA") >0 )
      temp <- temp[id1,]
      temp <- matrix(temp$interval,ncol=3)
      temp <- temp/rowSums(temp)
      temp <- colSums(temp)
      temp <-temp/sum(temp)
      weight <- temp
      names(weight) <- c("5'UTR","CDS","3'UTR")
      
      
      # density
      cds_id <- which(ct1$comp=="CDS")
      utr3_id <- which(ct1$comp=="UTR3")
      utr5_id <- which(ct1$comp=="UTR5")
      ct1$count[utr5_id] <- ct1$count[utr5_id]*weight["5'UTR"]
      ct1$count[cds_id] <- ct1$count[cds_id]*weight["CDS"]
      ct1$count[utr3_id] <- ct1$count[utr3_id]*weight["3'UTR"]
      
      # re-normalization
      featureSet <- as.character(unique(ct$Feature))
      for (i in 1:length(featureSet)) {
        id <- (ct1$Feature==featureSet[i])
        ct1$weight[id] <- ct1$count[id]/sum(ct1$count[id])
      }
      
      # stratch
      x <- cumsum(weight)
      ct1$pos[utr5_id] <- ct1$pos[utr5_id]*weight["5'UTR"] + 0
      ct1$pos[cds_id] <- ct1$pos[cds_id]*weight["CDS"] + x[1]
      ct1$pos[utr3_id] <- ct1$pos[utr3_id]*weight["3'UTR"] + x[2]
      
      # normalize
      featureSet <- as.character(unique(ct$Feature))
      for (i in 1:length(featureSet)) {
        id <- (ct1$Feature==featureSet[i])
        ct1$weight[id] <- ct1$count[id]/sum(ct1$count[id])
        id <- (ct2$Feature==featureSet[i])
        ct2$weight[id] <- ct2$count[id]/sum(ct2$count[id])
      }

      den_ct1 <- .countDensityCI(ct1, adjust, CI_interval)

      p1 <- 
        ggplot(den_ct1, aes(x=x, y = y, group = Feature))  +
        scale_x_continuous(minor_breaks = seq(0, 5, 1)) + 
        ggtitle("Distribution on mRNA") +
        theme(axis.ticks = element_blank(), axis.text.x = element_blank()) + 
        xlab("") + 
        ylab("Frequency") +
        geom_ribbon(aes(ymin =  rep(0, length(y)), ymax =  y, colour = factor(Feature), fill = factor(Feature)), alpha = 0.2 ) +
        geom_ribbon(aes(ymin =  confidence_down, ymax =  confidence_up, colour = factor(Feature), fill = factor(Feature)), alpha = 0.3 ) +
        geom_line(aes(y = y, colour = factor(Feature)), alpha = 0.8, size = 1) +
        annotate("text", x = sum(x1[1:2])/2, y = -0.6, label = "5'UTR") +
        annotate("text", x = sum(x1[2:3])/2, y = -0.6, label = "CDS") +
        annotate("text", x = x1[1]/2, y = -0.6, label = "Promoter (1kb)") +
        annotate("text", x = sum(x1[4:5])/2, y = -0.6, label = "Tail (1kb)") +
        annotate("text", x = sum(x1[3:4])/2, y = -0.6, label = "3'UTR")  + 
        geom_vline(xintercept=x1[1:4], linetype="dotted") + 
        annotate("rect", xmin = x1[1], xmax = x1[2], ymin = -0.22, ymax = -0.18, alpha = .99, colour = "black")+
        annotate("rect", xmin = x1[3], xmax = x1[4], ymin = -0.22, ymax = -0.18, alpha = .99, colour = "black")+
        annotate("rect", xmin = x1[2], xmax = x1[3], ymin = -0.26, ymax = -0.14, alpha = .2, colour = "black") +
        xlim(0,1) +
        theme(legend.position="bottom")
      
      den_ct2 <- .countDensityCI(ct2, adjust, CI_interval)
      
      p2 <- 
        ggplot(den_ct2, aes(x=x, y = y, group = Feature))  +
        ggtitle("Distribution on lncRNA")  +
        scale_x_continuous(minor_breaks = seq(6, 9, 1)) +
        theme(axis.ticks = element_blank(), axis.text.x = element_blank()) + 
        xlab("") + 
        ylab("Frequency") +
        geom_ribbon(aes(ymin =  rep(0, length(y)), ymax =  y, colour = factor(Feature), fill = factor(Feature)), alpha = 0.2 ) +
        geom_ribbon(aes(ymin =  confidence_down, ymax =  confidence_up, colour = factor(Feature), fill = factor(Feature)), alpha = 0.3 ) +
        geom_line(aes(y = y, colour = factor(Feature)), alpha = 0.8, size = 1) +
        annotate("text", x = sum(x2[1:2])/2, y = -0.6, label = "lncRNA") +
        annotate("text", x = x2[1]/2, y = -0.6, label = 'Promoter (1kb)') +
        annotate("text", x = sum(x2[2:3])/2, y = -0.6, label = "Tail (1kb)")  +
        geom_vline(xintercept=x2[1:2], linetype="dotted") + 
        annotate("rect", xmin = x2[1], xmax = x2[2], ymin = -0.22, ymax = -0.18, alpha = .99, colour = "black")+
        xlim(0,1) +
        theme(legend.position="bottom")
    }
    
  } 
  
  if (includeNeighborDNA==FALSE) {
    
    # remove DNA
    id1 <- which(match(ct1$comp,c("Front","Back")) >0 )
    ct1 <- ct1[-id1,]
    id2 <- which(match(ct2$comp,c("Front","Back")) >0 )
    ct2 <- ct2[-id2,]  
    
    # normalize feature
    featureSet <- as.character(unique(ct$Feature))
    for (i in 1:length(featureSet)) {
      id <- (ct1$Feature==featureSet[i])
      ct1$weight[id] <- ct1$weight[id]/sum(ct1$weight[id])
      
      id <- (ct2$Feature==featureSet[i])
      ct2$weight[id] <- ct2$weight[id]/sum(ct2$weight[id])
    }
    
    den_ct2 <- .countDensityCI(ct2, adjust, CI_interval)

    p2 <- 
      ggplot(den_ct2, aes(x=x, y = y, group = Feature))  +
      ggtitle("Distribution on lncRNA")  +
      scale_x_continuous(minor_breaks = seq(6, 9, 1)) +
      theme(axis.ticks = element_blank(), axis.text.x = element_blank()) + 
      xlab("") + 
      ylab("Frequency") +
      geom_ribbon(aes(ymin =  rep(0, length(y)), ymax =  y, colour = factor(Feature), fill = factor(Feature)), alpha = 0.2 ) +
      geom_ribbon(aes(ymin =  confidence_down, ymax =  confidence_up, colour = factor(Feature), fill = factor(Feature)), alpha = 0.3 ) +
      geom_line(aes(y = y, colour = factor(Feature)), alpha = 0.8, size = 1) +
      annotate("text", x = 0.5, y = -0.6, label = "lncRNA")+
      annotate("rect", xmin = 0, xmax = 1, ymin = -0.22, ymax = -0.18, alpha = .99, colour = "black")+
      theme(legend.position="bottom")
    
    if (rescaleComponent) {
      
      # normalization by length of components in mRNA
      # calculate relative length of each components
      temp <- unique(d[,c(1,3,4,5)])
      id1 <- which(match(temp$comp,c("Front","Back")) >0 )
      temp <- temp[-id1,] # remove DNA
      id1 <- which(match(temp$category,"mRNA") >0 )
      temp <- temp[id1,]
      temp <- matrix(temp$interval,ncol=3)
      temp <- temp/rowSums(temp)
      temp <- colSums(temp)
      temp <-temp/sum(temp)
      weight <- temp
      names(weight) <- c("5'UTR","CDS","3'UTR")
      
      
      # density
      cds_id <- which(ct1$comp=="CDS")
      utr3_id <- which(ct1$comp=="UTR3")
      utr5_id <- which(ct1$comp=="UTR5")
      ct1$count[utr5_id] <- ct1$count[utr5_id]*weight["5'UTR"]
      ct1$count[cds_id] <- ct1$count[cds_id]*weight["CDS"]
      ct1$count[utr3_id] <- ct1$count[utr3_id]*weight["3'UTR"]
      
      # re-normalization
      featureSet <- as.character(unique(ct$Feature))
      for (i in 1:length(featureSet)) {
        id <- (ct1$Feature==featureSet[i])
        ct1$weight[id] <- ct1$count[id]/sum(ct1$count[id])
      }
      
      # stratch
      x <- cumsum(weight)
      ct1$pos[utr5_id] <- ct1$pos[utr5_id]*weight["5'UTR"] + 0
      ct1$pos[cds_id] <- ct1$pos[cds_id]*weight["CDS"] + x[1]
      ct1$pos[utr3_id] <- ct1$pos[utr3_id]*weight["3'UTR"] + x[2]
      
      den_ct1 <- .countDensityCI(ct1, adjust, CI_interval)

      p1 <- 
        ggplot(den_ct1, aes(x=x, y = y, group = Feature))  +
        scale_x_continuous(minor_breaks = seq(0, 5, 1)) + 
        ggtitle("Distribution on mRNA") +
        theme(axis.ticks = element_blank(), axis.text.x = element_blank()) + 
        xlab("") + 
        ylab("Frequency") +
        geom_ribbon(aes(ymin =  rep(0, length(y)), ymax =  y, colour = factor(Feature), fill = factor(Feature)), alpha = 0.2 ) +
        geom_ribbon(aes(ymin =  confidence_down, ymax =  confidence_up, colour = factor(Feature), fill = factor(Feature)), alpha = 0.3 ) +
        geom_line(aes(y = y, colour = factor(Feature)), alpha = 0.8, size = 1) +
        annotate("text", x = x[1]/2, y = -0.6, label = "5'UTR") +
        annotate("text", x = x[1] + weight[2]/2, y = -0.6, label = "CDS") +
        annotate("text", x = x[2] + weight[3]/2, y = -0.6, label = "3'UTR") + 
        theme(legend.position="bottom") +
        geom_vline(xintercept= x[1:2], linetype="dotted") +
        annotate("rect", xmin = 0, xmax = x[1], ymin = -0.22, ymax = -0.18, alpha = .99, colour = "black")+
        annotate("rect", xmin = x[2], xmax = 1, ymin = -0.22, ymax = -0.18, alpha = .99, colour = "black")+
        annotate("rect", xmin = x[1], xmax = x[2], ymin = -0.26, ymax = -0.14, alpha = .2, colour = "black")
      
      
    } else {
      
      pos_adjust <- match(ct1$comp,c("UTR5","CDS","UTR3"))-1
      ct1$pos <- ct1$pos + pos_adjust
      
      den_ct1 <- .countDensityCI(ct1, adjust, CI_interval)
      
      p1 <- 
        ggplot(den_ct1, aes(x=x, y = y, group = Feature))  +
        scale_x_continuous(minor_breaks = seq(0, 5, 1)) + 
        ggtitle("Distribution on mRNA") +
        theme(axis.ticks = element_blank(), axis.text.x = element_blank()) + 
        xlab("") + 
        ylab("Frequency") +
        geom_ribbon(aes(ymin =  rep(0, length(y)), ymax =  y, colour = factor(Feature), fill = factor(Feature)), alpha = 0.2 ) +
        geom_ribbon(aes(ymin =  confidence_down, ymax =  confidence_up, colour = factor(Feature), fill = factor(Feature)), alpha = 0.3 ) +
        geom_line(aes(y = y, colour = factor(Feature)), alpha = 0.8, size = 1) +
        annotate("text", x = 0.5, y = -0.6, label = "5'UTR") +
        annotate("text", x = 1.5, y = -0.6, label = "CDS") +
        annotate("text", x = 2.5, y = -0.6, label = "3'UTR") + 
        geom_vline(xintercept=1:2, linetype="dotted") + 
        theme(legend.position="bottom") +
        annotate("rect", xmin = 0, xmax = 1, ymin = -0.22, ymax = -0.18, alpha = .99, colour = "black")+
        annotate("rect", xmin = 2, xmax = 3, ymin = -0.22, ymax = -0.18, alpha = .99, colour = "black")+
        annotate("rect", xmin = 1, xmax = 2, ymin = -0.26, ymax = -0.14, alpha = .2, colour = "black")
    }
  }
  
  suppressWarnings( 
    if (is.na(saveToPDFprefix)) {
      # return the result
      print("no figure saved ...")
    }  else {
      f1 <- paste(saveToPDFprefix,"_",transcript,"_Guitar.pdf",sep="")
      
      pdf(file=f1,width=8, height=4)
      if (transcript == 'both') {
        .multiplot(p1, p2, cols=2)
      }
      
      if (transcript == 'mRNA') {
        .multiplot(p1, cols=1)
      }
      
      if (transcript == 'lncRNA') {
        .multiplot(p2, cols=1)
      }
      dev.off()
      print(paste("Figures saved into",f1,"...", sep=" "))
    }
  )
  
  if (transcript == 'both') {
    suppressWarnings( .multiplot(p1, p2, cols=2))
  }

  if (transcript == 'mRNA') {
    suppressWarnings( .multiplot(p1, cols=1))
  }
  
  if (transcript == 'lncRNA') {
    suppressWarnings( .multiplot(p2, cols=1))
  }
  
}
