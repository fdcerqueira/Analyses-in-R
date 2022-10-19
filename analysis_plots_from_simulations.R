###load packages, if they are not present, it will install them
if(!require(ggplot2)){install.packages("ggplot2", repos = "http://cran.us.r-project.org")}
if(!require(meanShiftR)){install.packages("meanShiftR",repos = "http://cran.us.r-project.org")}
if(!require(scatterplot3d)){install.packages("scatterplot3d",repos = "http://cran.us.r-project.org")}
if(!require(dplyr)){install.packages("dplyr",repos = "http://cran.us.r-project.org")}
if(!require(stringr)){install.packages("stringr",repos = "http://cran.us.r-project.org")}
if(!require(rawr)){install.packages("rawr",repos = "http://cran.us.r-project.org")}

##load the packages
require(ggplot2)
require(meanShiftR)
require(scatterplot3d)
require(dplyr)
require(stringr)
require(rawr)

###function to be able to run rscript from terminal and set up the working directory the source file location
sort_final <- function() {
  cmdArgs <- commandArgs(trailingOnly = FALSE)
  needle <- "--file="
  match <- grep(needle, cmdArgs)
  if (length(match) > 0) {
    # Rscript
    return(normalizePath(sub(needle, "", cmdArgs[match])))
  } else {
    # 'source'd via R console
    return(normalizePath(sys.frames()[[1]]$ofile))
  }
}

script.dir <- dirname(sort_final())
setwd(script.dir)
print(getwd())
wd<-getwd() 

###use all fixDiv files, change here the prefix
fich <- sort(list.files(wd, pattern="HGTNRCfixDiv*")) 
#remove extension
names_interest <- sapply(strsplit(fich, ".dat"), `[`, 1) 
# get names for the output files
names_final<-sapply(names_interest,  sub, pattern = "HGTNRCfixDiv", replacement = "") 

#sequencia com nr ficheiros, and empty lists
fileNumbers <- seq(fich)
interesse<-list() #list with  files
#lists for the different graphs
div<-list() 
seg<-list()
terc<-list()

##loop ggplots fromfixDiv files
dir.create("ggplot_graphs")
for (i in fileNumbers) {
  interesse[[i]] <- read.delim(fich[i], sep="", header =F)
  
   div[[i]]<-ggplot(interesse[[i]],aes(x=V2,y=V9,color=as.factor(V1),group=as.factor(V1)))+
   geom_line()+
   ggtitle("Div U=0.1 HGT=1.0")+
     theme_bw()+xlab("generations") + #ylab("aquilo que quiser")
     guides(color=guide_legend(title="replicates"))+
     theme(text = element_text(size = 15))  
   seg[[i]]<-ggplot(interesse[[i]],aes(x=V2,y=V4,color=as.factor(V1),group=as.factor(V1)))+
     geom_line()+
     ggtitle("Fix U=0.1 HGT=1.0_")+
     theme_bw()+xlab("generations") +#ylab("aquilo que quiser")
     guides(color=guide_legend(title="replicates"))+
     theme(text = element_text(size = 15))  
  terc[[i]]<-ggplot(interesse[[i]],aes(x=V2,y=V6,color=as.factor(V1),group=as.factor(V1)))+
  geom_line()+ggtitle("a1+a2 U=0.1 HGT=1.0")+
    theme_bw()+ xlab("generations") + #ylab("aquilo que quiser")
    guides(color=guide_legend(title="replicates"))+
    theme(text = element_text(size = 15))  
  
  ggsave(path="ggplot_graphs",div[[i]], filename=paste0("div",names_final[[i]],".png"), width = 8, height = 8)
  ggsave(path="ggplot_graphs",seg[[i]], filename=paste0("fix",names_final[[i]],".png"), width = 8, height = 8)
  ggsave(path="ggplot_graphs",terc[[i]], filename=paste0("U",names_final[[i]],".png"), width = 8, height = 8)
}

##HGT get files
fich_hgt <- sort(list.files(wd, pattern="HGTNRCSites*")) ## all Sites files 
names_interest_hgt <- sapply(strsplit(fich_hgt, ".txt"), `[`, 1) #para remover a extensao
names_final<-sapply(names_interest_hgt,  sub, pattern = "HGTNRCSites", replacement = "") #get names for output

## empty lists for loop
datsit<-list()
alf1<-list()
alf2<-list()
dat1<-list()
scatter<-list()

####loop to get list with alfa1, alf2, as well as list with replicates, from all files
for (i in seq(fich_hgt)){
  datsit[[i]]<- read.delim(fich_hgt[i], header=F)
  alf1[[i]]<-split(datsit[[i]]$V3, datsit[[i]]$V1)
  alf2[[i]]<-split(datsit[[i]]$V4, datsit[[i]]$V1)
  dat1[[i]]<-split(datsit[[i]]$V1, datsit[[i]]$V1)
}


## to print each scatterplot3d to one pdf
dir.create("scatter3D")
##suplementary text.tile created on the same folder, with number (x) andname of the file corresponding to the File(x) title of each graph
write.table(names_interest, "scatter3d/file_order.txt")
pdf(file="scatter3d/scatter3d.pdf")
lapply(X = seq(dat1), FUN = function(i)
  lapply(X = seq(dat1[[i]]), 
         FUN = function(j)  
           scatterplot3d(alf1[[i]][[j]], 
                         alf2[[i]][[j]],
                         dat1[[i]][[j]],
                         grid=T,box=T,
                         pch=21,
                         bg = "blue",
                         xlab = "Alfa1",
                         ylab = "Alfa2",
                         zlab = "", 
                         main=paste("File",i,"Replicate",j))))

dev.off()

##epsilon variable: to use the DELTA value acordding the the file, later on the meanShift() function 
epsilon<- sapply(str_extract(fich_hgt,"DELTA[:digit:].(\\d+)"),`[`, 1)
epsilon<- sapply(str_extract(epsilon,"[:digit:].(\\d+)"),`[`, 1)
epsilon<-as.character(as.numeric(epsilon))

#sequencia com nr ficheiros
fileNumbers_hgt <- seq(fich_hgt)

##empty lists for loop with alfas
datsites<-list()
alfa1<-list()
alfa2<-list()

#loop to get alf1 and alfa2 from each file per replicate
for (i in fileNumbers_hgt){
  datsites[[i]]<- read.delim(fich_hgt[i], header=F,sep = "")
  alfa1[[i]]<-split(datsites[[i]]$V3, datsites[[i]]$V1)
   alfa2[[i]]<-split(datsites[[i]]$V4, datsites[[i]]$V1)
}


###merge into data.frame (with alfa1 and alfa2) nested lists, per file, per replicate
aux<-lapply(X = seq(alfa1), FUN = function(i)
  lapply(X = seq(alfa1[[i]]), 
         FUN = function(j) data.frame(alfa1=alfa1[[i]][[j]], 
                                      alfa2 = alfa2[[i]][[j]])))

###apply the meanshift fucntion with the DELTA being accordoing epslion[i], defined on lines 93, 94 
cluster<-lapply(X = seq(aux), FUN = function(i)
  lapply(X = seq(aux[[i]]), FUN = function(j) 
    meanShift(as.matrix(aux[[i]][[j]]), bandwidth = c(0.05,0.05),epsilonCluster = epsilon[i])))

###get max of cluesters per replicate per file
assignments<-lapply(X = seq(cluster), FUN = function(i)
  lapply(X = seq(cluster[[i]]), FUN = function(j) 
     max(data.frame(assignment=cluster[[i]][[j]]$assignment))))

##unlist this mess
final<-lapply(X = seq(assignments), FUN = function(i)
  lapply(X = seq(assignments[[i]]), FUN = function(j) 
    unlist(assignments[[i]][[j]])))

###bind rows and convert to dataframe
final_csv<-do.call(rbind, final)%>%
  as.data.frame()

#convert the numbers to character, as lists cant be written 
final_csv <- apply(final_csv,2,as.character)

##write final csv, the rows are the file numbers, and each column a replicate.
dir.create("clusters")
write.csv(final_csv, "clusters/max_clusters.csv" )



