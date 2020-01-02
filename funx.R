#################################################################
#################################################################
# This script file contains all auxiliary functions ...
# ...called by "Figure_5.Rmd", "Figure_S3.Rmd" and "Figure_4.Rmd"
#################################################################
#################################################################


#############################################
# FUNCTION 1: CLEAN UP THE CSV IMPORT
# function with 
# input: csv file
# output: clean version
#############################################

Cols_AllMissing <- function(df){ # helper function
  as.vector(which(colSums(is.na(df)) == nrow(df)))
}

clean_up <- function(data_gen){
  
  
  l<-which(data_gen$chr=="spot")
  data_gen$gene[l: (l+4)]<-data_gen$chr[l:(l+4)]
  
  ## remove tumor label, chr and position
  datG<-data_gen %>% select(-c("tumor", "chr", "position")) 
  
  ## Flip the array to get gene/columns and sample/rows
  datG<-t(datG)
  
  ## make the gene name the column name
  colnames(datG)<-datG[1,]
  datG<-datG[-1,]
  
  ## convert to a data frame
  datG<-as.data.frame(datG)
  
  ## create a new column with the sample names
  names<-rownames(datG)
  datG$sample<-names
  datG$sample<-as.factor(datG$sample)
  
  ## remove row names
  rownames(datG)<-NULL
  
  
  ## Make the VAF numeric entries
  for(k in 1:(which(colnames(datG)=="spot")-1)){
    datG[,k]<-as.numeric(as.character(datG[,k]))
  }
  
  
  datG<-select(datG, -Cols_AllMissing(datG))
  
  #d
  
  datG <- datG %>% mutate(spot = paste(slide, spot, sep = '_'))
  
  
   ###################################
   ## More processing: Singletons
   ###################################  
  
  for(j in 1:(which(colnames(datG)=="spot")-1)){

    ## Mutations only present in one spot -> evaluate whether to keep it
    vec.h<-which(datG[,j]>0.05)
  

    if(length(vec.h)==1 & any(datG[vec.h,j]<.06)){
      datG[vec.h,j]<-0

     # print(c(colnames(datG)[j],datG$spot[vec.h]))
    }

    ## Missed public mutation
    vec.h<-which(datG[,j]<.05 | is.na(datG[,j]))

    if(length(vec.h)<=0.1*dim(datG)[1]){

      datG[vec.h,j]<-mean(datG[-vec.h,j])
      #print(c(colnames(datG)[j],datG$spot[vec.h]))
    }
  }
  
    

  #########
  #########
  
  return(datG)
  
  
} 



#############################################
# FUNCTION 2: CLEAN UP THE INDIVIDUAL SLIDES
# input: DATA frame
# output: clean version
#############################################

clean_slide <- function(datG_S1){
  ## clean up the names
  b<-as.character(datG_S1$sample)
  for (k in 1: length(b)){
    
    b[k]<-sub('O.*', '', b[k])
    b[k]<-sub('\\..*', '', b[k])
    b[k]<-sub('_.*', '', b[k])
    b[k]<-sub('E', 'e', b[k])
    
  }
  
  if(length(b)!=length(unique(b))){
    print('ATTENTION: stop - in slide extraction') # that would not be good as we would have missed a "sub-spot"
  }
  
  ## update the names
  datG_S1$sample<-b
  
  datG_S1<-select(datG_S1, -Cols_AllMissing(datG_S1))
  
  return(datG_S1)

}



#############################################
# FUNCTION 3: Distances
# input: DATA frame
# output: clean version
#############################################

clean_dist<-function(distM){
  
## Create a vector that contains the names of the samples
a<-distM$spot
for (k in 1: length(a)){
  
  a[k]<-sub('E', 'e', a[k])
  a[k]<-sub('M.*', '', a[k])
  a[k]<-sub('h.*', '', a[k])
  a[k]<-sub('_.*', '', a[k])
  
  
}

if(length(a)!=length(unique(a))){
  print('stop') # that would not be good as we would have missed a "sub-spot"
}

distM$spot<-a
return(distM)
}

#############################################
# FUNCTION 4: Distance matrices
# input: DATA frame for a single section
# output: data frame containing both physical
# and genetic distance matrix
#############################################

dist_mat <-function(datG_S1, distG_S1, thr){

# identify the genetic information and create a distance matrix
ind<-which(colnames(datG_S1)=="spot")-1
# apply a threshold to call the genotypes

if(thr>0){ datG_S1[,1:ind]<-1*(datG_S1[,1:ind]>thr) }

d_S1 <- dist(datG_S1[,1:ind], method = "manhattan") # distance matrix
dist_gen<-as.vector(d_S1)
distG_S1<-distG_S1[,-1]
dist_space<-as.vector(as.matrix(distG_S1))
dist_space <- dist_space[!is.na(dist_space)]

dist<-data.frame(space=dist_space, gen=dist_gen)

return(dist)

} 



#############################################
# FUNCTION 5: Master cleaning function
# input: DATA from entire tumor
# output: processed genetic and spatial information
# and genetic distance matrix
#############################################


tumor_process <- function(tumor){
  
  
  #####################################
  ## Which tumor to analyze
  #####################################
  
  
  
  if(tumor=="E"){
  
  #####################################
  ## TUMOR E
  #####################################
  X<-read.csv(file="./tumor_E/tumor_E_gen.csv", stringsAsFactors=FALSE)
  X1<-read.csv(file="./tumor_E/tumor_E_distance_A3.csv", stringsAsFactors = FALSE)
  X2<-read.csv(file="./tumor_E/tumor_E_distance_A9.csv", stringsAsFactors = FALSE)
  ind<-c("A3", "A9") # slide names
  
  } else if(tumor=="M"){
  
  #####################################
  ## TUMOR M
  #####################################
  X<-read.csv(file="./tumor_M/tumor_M_gen.csv", stringsAsFactors=FALSE)
  X1<-read.csv(file="./tumor_M/tumor_M_distance_A3.csv", stringsAsFactors = FALSE)
  X2<-read.csv(file="./tumor_M/tumor_M_distance_A4.csv", stringsAsFactors = FALSE)
  ind<-c("A3", "A4") # slide names
  
  
  } else if(tumor=="D"){
    
  #####################################
  ## TUMOR D
  #####################################
  X<-read.csv(file="./tumor_D/tumor_D_gen.csv", stringsAsFactors=FALSE)
  X1<-read.csv(file="./tumor_D/tumor_D_distance_A4.csv", stringsAsFactors = FALSE)
  X2<-read.csv(file="./tumor_D/tumor_D_distance_A3.csv", stringsAsFactors = FALSE)
  ind<-c("A4", "A3") # slide names
  
  } else if(tumor=="K"){
    
  #####################################
  ## TUMOR Kad
  #####################################
  X<-read.csv(file="./tumor_Kad/tumor_Kad_gen.csv", stringsAsFactors=FALSE)
  X1<-read.csv(file="./tumor_Kad/tumor_Kad_distance_A21.csv", stringsAsFactors = FALSE)
  X2<-read.csv(file="./tumor_Kad/tumor_Kad_distance_A22.csv", stringsAsFactors = FALSE)
  ind<-c("A21", "A22") # slide names
  
  
  } else if(tumor=="T"){
    
  #####################################
  ## TUMOR T
  #####################################
  X<-read.csv(file="./tumor_T/tumor_T_gen.csv", stringsAsFactors=FALSE)
  X1<-read.csv(file="./tumor_T/tumor_T_distance_A7.csv", stringsAsFactors = FALSE)
  X2<-read.csv(file="./tumor_T/tumor_T_distance_A10.csv", stringsAsFactors = FALSE)
  ind<-c("A7", "A10") # slide names
  
  } else if(tumor=="R"){
  
  #####################################
  ## TUMOR Knew
  #####################################
  X<-read.csv(file="./tumor_Knew/tumor_Knew_gen.csv", stringsAsFactors=FALSE)
  X1<-read.csv(file="./tumor_Knew/tumor_Knew_distance_A7.csv", stringsAsFactors = FALSE)
  X2<-read.csv(file="./tumor_Knew/tumor_Knew_distance_A7.csv", stringsAsFactors = FALSE)
  ind<-c("A7","A7") # slide names
  
  
  } else if(tumor=="W"){
    
  #####################################
  ## TUMOR W
  #####################################
  X<-read.csv(file="./tumor_W/tumor_W_gen.csv", stringsAsFactors=FALSE)
  X1<-read.csv(file="./tumor_W/tumor_W_distance_A5.csv", stringsAsFactors = FALSE)
  X2<-read.csv(file="./tumor_W/tumor_W_distance_A6.csv", stringsAsFactors = FALSE)
  ind<-c("A5","A6") # slide names
  
  } else if(tumor=="H"){
  #####################################
  ## TUMOR H
  #####################################
  X<-read.csv(file="./tumor_H/tumor_H_gen.csv", stringsAsFactors=FALSE)
  X1<-read.csv(file="./tumor_H/tumor_H_distance_A10.csv", stringsAsFactors = FALSE)
  X2<-read.csv(file="./tumor_H/tumor_H_distance_A11.csv", stringsAsFactors = FALSE)
  ind<-c("A10","A11") # slide names
  
  } else if(tumor=="F"){
  
  #####################################
  ## TUMOR F
  #####################################
  X<-read.csv(file="./tumor_F/tumor_F_gen.csv", stringsAsFactors=FALSE)
  X1<-read.csv(file="./tumor_F/tumor_F_distance_A5.csv", stringsAsFactors = FALSE)
  X2<-read.csv(file="./tumor_F/tumor_F_distance_A6.csv", stringsAsFactors = FALSE)
  ind<-c("A5","A6") # slide names
  
  
  } else if(tumor=="C"){
  #####################################
  ## TUMOR C
  #####################################
  X<-read.csv(file="./tumor_C/tumor_C_gen.csv", stringsAsFactors=FALSE)
  X1<-read.csv(file="./tumor_C/tumor_C_distance_A9.csv", stringsAsFactors = FALSE)
  X2<-read.csv(file="./tumor_C/tumor_C_distance_A26.csv", stringsAsFactors = FALSE)
  ind<-c("A9","A26") # slide names
  
  } else if(tumor=="U"){
  #####################################
  ## TUMOR U
  #####################################
  X<-read.csv(file="./tumor_U/tumor_U_gen.csv", stringsAsFactors=FALSE)
  X1<-read.csv(file="./tumor_U/tumor_U_distance_A4.csv", stringsAsFactors = FALSE)
  X2<-read.csv(file="./tumor_U/tumor_U_distance_A5.csv", stringsAsFactors = FALSE)
  ind<-c("A4","A5") # slide names
  
  
  } else if(tumor=="J"){
    
  #####################################
  ## TUMOR J
  #####################################
  X<-read.csv(file="./tumor_J/tumor_J_gen.csv", stringsAsFactors=FALSE)
  X1<-read.csv(file="./tumor_J/tumor_J_distance_B6.csv", stringsAsFactors = FALSE)
  X2<-read.csv(file="./tumor_J/tumor_J_distance_B8.csv", stringsAsFactors = FALSE)
  ind<-c("B6","B8") # slide names
  
  
  }
  
  
  
  #####################################
  ## Cleaning necessary?
  #####################################
  
  
  if(length(intersect(colnames(X1), colnames(X)))==dim(X1)[1] & length(intersect(colnames(X2), colnames(X)))==dim(X2)[1]){
    
    clean=0; # do not need to clean the names
    
  } else {
    
    clean=1; # need to clean the names
  }
  
  
  #####################################
  ## Prep: physical distance
  #####################################
  
  if(clean==1){
    ## Physical distance S1
    distG_S1<-clean_dist(X1)
    
    ## Physical distance S2
    distG_S2<-clean_dist(X2)
  } else{
    
    
    distG_S1<-X1
    distG_S2<-X2
    
  }
  
  #####################################
  # Prep: sequencing data
  #####################################
  
  ## import the spreadsheet
  data_gen<-X
  
  ## clean up
  datG<-clean_up(data_gen)
  
  
  # Extract the first slide
  datG_S1<-datG %>% filter(datG$slide==ind[1])
  
  # Extract the second slide
  datG_S2<-datG %>% filter(datG$slide==ind[2])
  
  # clean the slides
  if(clean==1){
    datG_S1<-clean_slide(datG_S1)
    datG_S2<-clean_slide(datG_S2)
  } else{
    
    datG_S1$sample<-as.character(datG_S1$sample)
    datG_S2$sample<-as.character(datG_S2$sample)
  }
  
  # order sample names between gen and dist 
  datG_S1<-datG_S1[match(distG_S1$spot, datG_S1$sample),]
  datG_S2<-datG_S2[match(distG_S2$spot, datG_S2$sample),]
  
  out<-list(datG, distG_S1, distG_S2,datG_S1,datG_S2 )
  
  
}




#############################################
#############################################
# FUNCTION 6: MP Tree function
#############################################
#############################################




MP_tree_function <- function(thr, tumor.summary, nam, n.BS){
  
  ind=which(tumor.summary$name==nam)
  reroot=as.character(tumor.summary$reroot[ind])
  
  ##################################
  # Process the samples
  ##################################
  
  Xout<-tumor_process(nam)
  
  # Extract the samples
  datG<-Xout[[1]]
  distG_S1<-Xout[[2]]
  distG_S2<-Xout[[3]]
  datG_S1<-Xout[[4]]
  datG_S2<-Xout[[5]]
  
  # Convert to a phylo object
  ind<-which(colnames(datG)=="spot")
  dat<-datG[,1:ind]
  rownames(dat)<-dat$spot
  dat$spot<-NULL
  
  # remove issues with homoplasia
  dat<-homoplasia_clean(dat, nam)
  
  dat<-(dat>thr)
  dat1<-as.phyDat(dat, type="USER", levels = c(TRUE, FALSE))
  
  # Extract the unique genotypes
  ind<-rownames(unique(dat))
  namz<-data.frame(original=rownames(dat),
                   group=character(dim(dat)[1]),
                   subclone=numeric(dim(dat)[1]),
                   stringsAsFactors = FALSE)
  
  
  for(i in 1:length(ind)){
    
    a<- which(apply(dat, 1, function(x) identical(dat[ind[i],], x)))
    namz$group[a]<-ind[i]
    namz$subclone[a]<-i
    
  }
  
  ##################################
  # MP tree
  ##################################
  
  # Ratchet algorithm for maximum parsimony tree
  treeMP <- pratchet(dat1, trace=0)
  
  # If tree is not binary, convert to binary
  treeMP<-multi2di(treeMP)
  
  # Extract branch length
  treeMP <- acctran(treeMP, dat1)
  
  # Reroot the tree
  if(reroot!="none"){treeMP<-root(treeMP, outgroup=reroot, resolve.root = TRUE)}
  
  # The unique genotypes in proper order
  namz=namz[match(treeMP$tip.label,namz$original),] 
  
  
  ##################################
  ## Bootstrapping
  ##################################
  
  # Bootstrap trees
  BStrees <- bootstrap.phyDat(dat1, pratchet, bs = n.BS, multicore=TRUE)
  
  # Resolve potential non-binary multifurcations
  BStrees<-multi2di(BStrees)
  
  # Reroot the trees
  if(reroot!="none"){BStrees<-root(BStrees, outgroup=reroot, resolve.root = TRUE)}
  
  
  ##################################
  ## Plot the trees
  ##################################
  
  #palette(rainbow(length(unique(namz$group))))
  #plotBS(treeMP, BStrees, "phylogram", bs.col="black", main=paste("Phylogeny of tumor ", nam, sep = ""), cex=0.75,  p=0.00001, tip.color=namz$subclone)
  
  
  ################################
  # MUTATION PANEL HEATMAP
  ################################
  
  #heatmap.2(t(dat*1), distfun = dist, hclustfun = hclust,  scale = "none", col = cm.colors(2), key=TRUE, tracecol=NA, main = "Mutation panel", key.title="Mutation Key", revC = TRUE, cexRow=.5)
  
  ################################
  ################################
  # LIST OF SUBCLONES
  ################################
  ################################
  
  out<-list(namz, treeMP, BStrees, dat)
  
}






#############################################
#############################################
# FUNCTION 7: Avoid homoplasia in the dataset
#############################################
#############################################

homoplasia_clean <- function(dat, nam){
  
  if(nam=="D"){
    
    ## A4_3, remove spot
    h<-which(rownames(dat)=="A4_3")
    dat<-dat[-h,]
    
    # A3_38: remove the two mutations
    h<-which(rownames(dat)=="A3_38")
    ik<-which(colnames(dat)=="TP53_2")
    jk<-which(colnames(dat)=="RC3H1")
    dat[h,c(ik,jk) ]<-0
    
    # A3_37: add two mutations
    h<-which(rownames(dat)=="A3_37")
    ik<-which(colnames(dat)=="CD109")
    jk<-which(colnames(dat)=="CCKBR")
    dat[h,c(ik,jk) ]<-1
    
    dat[is.na(dat)]<-0
    
  }
  
  
  if(nam=="E"){
    # decided to remove A3_13
    h<-which(rownames(dat)=="A3_13")
    ik<-which(dat[h,]>.05)
    jk<-which(dat[h,ik]<.1)
    dat[h,ik[jk]]<-0
    
  }
  
  if(nam=="F"){
    ## REmove spot A6_25 (mixture)
    h<-which(rownames(dat)=="A6_25")
    dat<-dat[-h,]
    
    #Remove problematic mutations in A5_28
    h<-which(rownames(dat)=="A5_28")
    ij<-which(colnames(dat)=="GRIK5")
    dat[h,ij]<-0
    
    #Remove problematic mutations in A5_12
    h<-which(rownames(dat)=="A5_12")
    ij<-which(colnames(dat)=="IGSF10")
    dat[h,ij]<-0
  }
  
  if(nam=="H"){
    
    ## Remove A11_17
    h<-which(rownames(dat)=="A11_17")
    dat<-dat[-h, ]
    
  }
  
  if(nam=="R"){
    
    h<-which(rownames(dat)=="A7_7")
    ik<-which(dat[h,]>.05)
    jk<-which(dat[h,ik]<.2)
    dat[h,ik[jk]]<-0
    
  }
  
  
  if(nam=="T"){
    
    # remove a spot that causes problems: A7_1
    h<-which(rownames(dat)=="A7_1")
    dat<-dat[-h,]
    
    # In A10_7 set two problem mutations to 0
    h<-which(rownames(dat)=="A10_7")
    ik<-c(which(colnames(dat)=="SLC8A1"), which(colnames(dat)=="DOPEY1"))
    dat[h,ik]<-0
    
  }
  
  
  if(nam=="W"){
    
    ## Problems with homoplasia in A6_24A
    h<-which(rownames(dat)=="A6_24A")
    ik<-which(colnames(dat)=="ZNF23")
    dat[h,ik]<-0
    
  }  
  
  
  
  out=dat
  
  
}

