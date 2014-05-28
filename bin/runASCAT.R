##First read in the arguments listed at the command line
args=(commandArgs(TRUE))

##args is now a list of character vectors
## First check to see if arguments are passed.
if(length(args)==0){
  stop("No arguments supplied.")
} else{
	ascat_lib = args[1]
  SNP_pos_with_names = args[2]
  GCcorrect_file = args[3]
	tumour_sample = args[4]
	tumour_count_file = args[5]
	normal_sample = args[6]
	normal_count_file = args[7]
	gender = args[8]
  rdat_out = args[9]
  purity = as.numeric(args[10])
  ploidy = as.numeric(args[11])
  refchrs = args[12]
 
  refCN = ifelse(args[13]=="NA",NA,as.numeric(args[13])) 
  if(is.na(refchrs)||refchrs=="NA") {
  	refchrs = NA
	}

}

## checkpointing: if the RData file exists, only rerun the last step of ASCAT
if(length(dir(pattern=rdat_out))==0) {

  ## if the counts exist, skip this step - this can be removed later, but is easier for testing
  if(length(dir(pattern=tumour_count_file))==0) { 
		stop("Tumour count file missing or inaccessible")
  }

  ## if counts of normals exist, this can be skipped.. Again, easier for testing
  if(length(dir(pattern=normal_count_file))==0) {
  	stop("Normal count file missing or inaccessible")
  }

  ## create ASCAT input from read counts
  ## if the LogR file exists, skip this step - this can be removed later, but is easier for testing
  if(length(dir(pattern=paste(tumour_sample, ".tumour.LogR.txt",sep="")))==0) { 

    tumorcounts = read.table(tumour_count_file,sep="\t")
    normalcounts = read.table(normal_count_file,sep="\t")

    SNPpos = matrix(nrow = dim(normalcounts)[1],ncol = 2)
    rownames(SNPpos) = paste("snp",1:dim(SNPpos)[1],sep="")
    colnames(SNPpos) = c("Chr","Position")
    SNPpos[,1] = as.vector(normalcounts[,1])
    SNPpos[,2] = normalcounts[,2]

    ## This gives every SNP a name (to be able to do GC correction)				
    SNPposWithNames = read.table(SNP_pos_with_names,sep="\t",header=T,row.names=1)

    ctrans = 1:24
    names(ctrans)=c(1:22,"X","Y")
    newnames = ctrans[as.vector(SNPposWithNames[,1])]*1000000000+SNPposWithNames[,2]
    newnamesSEQ = ctrans[as.vector(SNPpos[,1])]*1000000000+as.numeric(SNPpos[,2])

    namestrans = rownames(SNPposWithNames)
    names(namestrans)=newnames

    rownames(SNPpos) = newnamesSEQ
    SNPpos = SNPpos[rownames(SNPpos)%in%newnames,]
    rownames(SNPpos)=namestrans[rownames(SNPpos)]

    rownames(tumorcounts) = newnamesSEQ
    tumorcounts = tumorcounts[rownames(tumorcounts)%in%newnames,]
    rownames(tumorcounts)=namestrans[rownames(tumorcounts)]
    rownames(normalcounts) = newnamesSEQ
    normalcounts = normalcounts[rownames(normalcounts)%in%newnames,]
    rownames(normalcounts)=namestrans[rownames(normalcounts)]

    Tumor_BAF = matrix(nrow = dim(normalcounts)[1],ncol = 1)
    rownames(Tumor_BAF) = rownames(SNPpos)
    colnames(Tumor_BAF) = tumour_sample
    acgt = tumorcounts[,c(3:6)]
    acgts = t(apply(acgt,1,sort))
    Tumor_BAF[,1] = acgts[,4]/(acgts[,3]+acgts[,4])
    Tumor_BAF[,1] = ifelse(runif(length(Tumor_BAF[,1]))<0.5,Tumor_BAF[,1],1-Tumor_BAF[,1])
    Tumor_BAF[is.nan(Tumor_BAF)]=NA

    Germline_BAF = matrix(nrow = dim(normalcounts)[1],ncol = 1)
    rownames(Germline_BAF) = rownames(SNPpos)
    colnames(Germline_BAF) = tumour_sample
    acgt = normalcounts[,c(3:6)]
    acgts = t(apply(acgt,1,sort))
    Germline_BAF[,1] = acgts[,4]/(acgts[,3]+acgts[,4])
    Germline_BAF[,1] = ifelse(runif(length(Germline_BAF[,1]))<0.5,Germline_BAF[,1],1-Germline_BAF[,1])
    Germline_BAF[is.nan(Germline_BAF)]=NA


    Tumor_LogR = matrix(nrow = dim(normalcounts)[1],ncol = 1)
    Germline_LogR = matrix(nrow = dim(normalcounts)[1],ncol = 1)
    rownames(Tumor_LogR) = rownames(SNPpos)
    colnames(Tumor_LogR) = tumour_sample
    rownames(Germline_LogR) = rownames(SNPpos)
    colnames(Germline_LogR) = tumour_sample
    Tumor_LogR[,1] = log(tumorcounts[,7]/normalcounts[,7],2)
    Germline_LogR[,1] = 0
    Tumor_LogR[is.infinite(Tumor_LogR)]=NA
    if(gender=="XY") {
      Tumor_LogR[SNPpos[,1]=="X",1] = Tumor_LogR[SNPpos[,1]=="X",1]-1
      Germline_LogR[SNPpos[,1]=="X",1] = Germline_LogR[SNPpos[,1]=="X",1]-1
    }
    Tumor_LogR[,1] = Tumor_LogR[,1] - median(Tumor_LogR[,1],na.rm=T)
    # set regions with 0 reads in tumor and normal to a LogR of 0.
    Tumor_LogR[is.na(Tumor_LogR[,1]),1] = 0

    # limit the number of digits:
    Tumor_LogR = round(Tumor_LogR,4)
    Tumor_BAF = round(Tumor_BAF,4)
    Germline_LogR = round(Germline_LogR,4)
    Germline_BAF = round(Germline_BAF,4)

    # write output to files
    write.table(cbind(SNPpos,Tumor_LogR),paste(tumour_sample,".tumour.LogR.txt",sep=""),sep="\t",row.names=T,col.names=NA,quote=F)
    write.table(cbind(SNPpos,Tumor_BAF),paste(tumour_sample,".tumour.BAF.txt",sep=""),sep="\t",row.names=T,col.names=NA,quote=F)
    write.table(cbind(SNPpos,Germline_LogR),paste(tumour_sample,".normal.LogR.txt",sep=""),sep="\t",row.names=T,col.names=NA,quote=F)
    write.table(cbind(SNPpos,Germline_BAF),paste(tumour_sample,".normal.BAF.txt",sep=""),sep="\t",row.names=T,col.names=NA,quote=F)
  }



  ## run ASCAT
  source(paste(ascat_lib,"ascat.R",sep="/"))

  ascat.bc = ascat.loadData(paste(tumour_sample,".tumour.LogR.txt",sep=""),paste(tumour_sample,".tumour.BAF.txt",sep=""),
                            paste(tumour_sample,".normal.LogR.txt",sep=""),paste(tumour_sample,".normal.BAF.txt",sep=""),gender=gender)

  ascat.bc = ascat.GCcorrect(ascat.bc, GCcorrect_file)

  ascat.plotRawData(ascat.bc)

  ascat.bc = ascat.aspcf(ascat.bc)

  ascat.plotSegmentedData(ascat.bc)

## if RData file exists, read it in
} else {
	
  load(rdat_out)

  ## this can be removed later, but for now, reload ASCAT code before rerunning..
  source(paste(ascat_lib,"ascat.R",sep="/"))

}

ascat.output = ascat.runAscat(ascat.bc, gamma = 1, rho_manual = purity, psi_manual = ploidy)

if(!is.na(refchrs)) {

  ## platform-specific parameter (1 for NGS, 0.55 for SNP6)
  gamma = 1

  dchr = strsplit(refchrs,",")[[1]]

  hetsnps = ascat.bc$SNPpos[!is.na(ascat.bc$Germline_BAF[,1]) & 
    ascat.bc$Germline_BAF[,1] >= 0.3 & ascat.bc$Germline_BAF[,1] <= 0.7,]
  hetIDs=which(rownames(ascat.bc$SNPpos)%in%rownames(hetsnps))
  hetIDs2=hetIDs[ascat.bc$Tumor_BAF_segmented[[1]]==0.5]
  selIDs = hetIDs2[as.vector(ascat.bc$SNPpos[hetIDs2,1])%in%dchr]

  rho = ascat.output$aberrantcellfraction[1]
  psi = ascat.output$psi[1]
  avCN = (2*rho-2+(2*(1-rho)+rho*psi)*2^(mean(ascat.bc$Tumor_LogR[selIDs,1],na.rm=T)/gamma))/rho

  ## update ploidy and rerun ASCAT
  ploidy = ploidy / avCN * refCN
  ascat.output = ascat.runAscat(ascat.bc, gamma = 1, rho_manual = purity, psi_manual = ploidy)

  ## update ploidy once more (fine-tuning)
  rho = ascat.output$aberrantcellfraction[1]
  psi = ascat.output$psi[1]
  avCN = (2*rho-2+(2*(1-rho)+rho*psi)*2^(mean(ascat.bc$Tumor_LogR[selIDs,1],na.rm=T)/gamma))/rho
  ploidy = ploidy / avCN * refCN
  ascat.output = ascat.runAscat(ascat.bc, gamma = 1, rho_manual = purity, psi_manual = ploidy)

}



## make output files
if(!is.null(ascat.output$nA)) {

  ## make output for Catherine's visualization tool
  gCN = matrix(nrow = dim(ascat.bc$Tumor_LogR)[1], ncol = 9)
  rownames(gCN) = rownames(ascat.bc$Tumor_LogR)
  colnames(gCN) = c("Chromosome","Position","Log R", "segmented LogR", "BAF", "segmented BAF", "Copy number", "Minor allele", "Raw copy number")
  gCN[,1]=as.vector(ascat.bc$SNPpos[,1])
  gCN[,2]=ascat.bc$SNPpos[,2]
  gCN[gCN[,1]=="X",1]=23
  gCN[,3]=ascat.bc$Tumor_LogR[,1]
  gCN[,4]=ascat.bc$Tumor_LogR_segmented[,1]
  gCN[,5]=ascat.bc$Tumor_BAF[,1]
  gCN[rownames(ascat.bc$Tumor_BAF_segmented[[1]]),6]=ascat.bc$Tumor_BAF_segmented[[1]][,1]
  seg = ascat.output$segments

  majorallele = NULL
  minorallele = NULL
  for (i in 1:dim(seg)[1]) {
    start = which(ascat.bc$SNPpos[,1]==seg[i,2] & ascat.bc$SNPpos[,2]==seg[i,3])[1]
    end = which(ascat.bc$SNPpos[,1]==seg[i,2] & ascat.bc$SNPpos[,2]==seg[i,4])[1]
    majorallele=c(majorallele,rep(seg[i,5],end-start+1))
    minorallele=c(minorallele,rep(seg[i,6],end-start+1))
  }

  gCN[,7]=majorallele+minorallele
  gCN[,8]=minorallele

  rho = ascat.output$aberrantcellfraction[1]
  psi = ascat.output$psi[1]

  gCN[,9]=(2*rho - 2 + (2*(1-rho)+rho*psi)*2^(ascat.bc$Tumor_LogR[,1]/0.55))/rho

  write.table(gCN,paste("CopyNumber",colnames(ascat.output$nA)[1],".txt",sep=""),sep="\t",row.names=T,col.names=NA,quote=F)

  ## make input for Caveman:
  cavemanSegs = cbind(seg[,2],
          seg[,3],
          seg[,4],2,1,
          seg[,5]+seg[,6],
          seg[,6])

  cavemanSegs = rbind(cavemanSegs,c("Y",2649553,59033810,0,0,0,0))

  rownames(cavemanSegs) = 1:dim(cavemanSegs)[1]

  write.table(cavemanSegs,paste("CopyNumberCaveman",colnames(ascat.output$nA)[1],".csv",sep=""),row.names=T,col.names=F,sep=",",quote=F)

  normalContamination = 2*(1-rho)/(2*(1-rho)+rho*ascat.output$ploidy[1])
  
  ss = matrix(ncol=1,nrow=4)
  rownames(ss) = c("NormalContamination","Ploidy","rho","psi")
  ss[,1] = c(normalContamination,ascat.output$ploidy,rho,psi)
    
  write.table(ss,paste("SampleStatistics",colnames(ascat.output$nA)[1],".csv",sep=""),row.names=T,col.names=F,quote=F)

}

## remove purity, ploidy and diploidchromosomes from RData file, otherwise we can't reiterate
rm(purity)
rm(ploidy)
rm(refchrs)
rm(refCN)

## save file
save.image(rdat_out)

q(save="no")


