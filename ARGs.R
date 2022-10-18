library(ggplot2)
library(cowplot2)
library(dplyr)
library(data.table)
library(randomcolor)
library(cowplot)
library(vegan)
library(labdsv)
library(ggfortify)
library(factoextra)

setwd("~/Desktop/dorde/results")
#import table relative abundance

abundances<-read.csv("~/Desktop/dorde/Relative abundances.csv")

#remove NS samples
abundances<-abundances[-c(4:6)]
abundances<-abundances[-c(7:9)]
abundances<-abundances[-c(10:12)]

#remove 1st row
abundances<-abundances[-1,]
colnames(abundances)<-abundances[1,]
abundances<-abundances[-1,]

#subsitute empty cells with zeros

abundances_values <- abundances[-c(1:3)] 
abundances_values <-abundances_values%>%mutate(across(-1,  ~ as.numeric(replace(., . == '', 0))))
abundnce4A <- ifelse(abundances_values$`4A`=="", 0, abundances_values$`4A`)

#colbind
abundances1<-cbind(abundances[c(1:3)],abundnce4A,abundances_values[-1])

#change A4 column to numeric
abundances1$abundnce4A<-as.numeric(as.character(abundances1$abundnce4A))

#change name of colum 4A
names(abundances1)[names(abundances1) == "abundnce4A"] <-"4A"

#remove taxonomy
abundance_arg<-abundances1[!(abundances1$group=="Taxanomic" | abundances1$group=="16S rRNA"),]

#get colors for each group
no_of_colors <- length(unique(abundance_arg$group))
palette <- distinctColorPalette(no_of_colors) 

#set colours
abundance_arg$colour <-
  case_when(
    abundance_arg$group =="Other" ~ "black",
    abundance_arg$group =="Aminoglycoside" ~ "blue",
    abundance_arg$group =="Beta Lactam"~"gold2",
    abundance_arg$group =="Integrons"~"olivedrab1",
    abundance_arg$group =="MDR"~"skyblue",
    abundance_arg$group =="MGE"~"#D4DE52",
    abundance_arg$group =="MLSB"~"darkblue",
    abundance_arg$group =="Other"~"black",
    abundance_arg$group =="Quinolone"~"#E06667",
    abundance_arg$group =="Sulfonamide"~"khaki1",
    abundance_arg$group =="Tetracycline"~"#E1B096",
    abundance_arg$group =="Trimethoprim"~"darkviolet",
    TRUE ~ abundance_arg$group
  )

#change to long format table
long <- melt(setDT(abundance_arg), id.vars = c("group","gene","assay","colour"), variable.name = "sample")

#ggplot stacked barplot ARGs

args <- long%>%group_by(sample)%>%
  mutate(frequency = value / sum(value))

args$group<-factor(args$group)


args1 <- args %>% 
  group_by(group,gene, colour)

#create a vector with colors, and give each element a name corresponding to genus 
col <-as.character(args$colour)
names(col) <- as.character(args$group)

pdf("relative_ARGs_group.pdf",width = 10, height = 7)
ggplot(args, aes(x = sample, y = frequency, fill = group))+
  scale_y_continuous(labels = c("0","25","50","75","100"))+
  geom_bar(stat="identity", position = "stack")+
  scale_fill_manual(values=col)+
  guides(fill = guide_legend(reverse = TRUE, keywidth = 1, keyheight = 1)) +
  ylab("Relative Abundance (%)  \n")+
  theme_minimal(base_size =14)+
  #facet_wrap(~Treatment, scales="free_x", ncol = 5, labeller = label_both)+
  theme(axis.text.x = element_text(angle = 0,vjust = 1.2), 
        panel.grid.major = element_line(colour = "white"), 
  )+
  labs(title = NULL)+
  guides(fill=guide_legend(ncol=1))
dev.off()

#dotplot
pdf("dot_plot_ARGs.pdf",width = 10, height = 7)
ggplot(data = args, aes(x = sample, y = value, fill=group)) +
  geom_point(colour="black",pch=21, size=3)+
  labs(x = "", y = "Gene abundance/16S rDNA gene abundance")+ 
  geom_text_repel(data = .%>% 
                    mutate(label=gene),aes(label=label))+
  theme_cowplot(font_size = 14,)+
  theme(strip.background = element_blank())
dev.off()

#alpha and beta diversity of ARGs
#alpha
alpha1 <- subset(abundance_arg, select = -c(group,assay, colour))

gene<-alpha1$gene

alpha1<-alpha1[,-1]
alpha2<-t(alpha1)
colnames(alpha2)<-gene

#convert dataframe to numeric matrix
alpha3<-apply(as.matrix(alpha2), 2, as.numeric)
rownames(alpha3)<-rownames(alpha2)

group<-ifelse(grepl("[4-6]A",rownames(alpha2)),"S-plant",
              ifelse(grepl("10|11|12",rownames(alpha2)),"S-plant-T34",
                     ifelse(grepl("16|17|18",rownames(alpha2)),"S-control","inlet")))

div<-c("shannon","simpson","invsimpson")
alpha4<-data.frame(matrix(nrow=length(group)))

for (i in div){
  alpha4[i]<-diversity(alpha3, index =i)
} 

alpha4<-alpha4[,-1]
colnames(alpha4)<- div

data1 <- cbind(rownames(alpha3), group,alpha4)

rownames(data1)<-data1$`rownames(alpha3)`
data1<-as.data.frame(data1)


measures<-c("shannon","simpson","invsimpson")
groupy<-c("S-plant","S-plant-T34","S-control","inlet")


for (i in measures){ 
print(pairwise.wilcox.test(data1[,i], data1$group, p.adjust.method="hochberg"))
}

#beta diversity
bet_bray<-log10(alpha3+1)

#distances:bray and jacccard
bray<-vegdist(bet_bray, method="bray")
jaccard<-vegdist(alpha3, method="jaccard", binary = TRUE)

#dispersion

disper_bray<-betadisper(bray,group)
disper_jaccard <- betadisper(jaccard,group)
  
#permutests

permu_bray<-permutest(disper_bray)
permu_jaccard<-permutest(disper_jaccard)
  
#PERMANOVA
permanova_bray<-adonis(bray~ group)
permanova_jaccard<-adonis(jaccard ~ group)

#NMDS
NMDS_bray<-metaMDS(alpha3,distance = "bray")
NMDS_jaccard<-metaMDS(alpha3,distance = "jaccard", binary="true")

#stressplot
stress_bray<-stressplot(NMDS_bray) 
stress_jaccard<-stressplot(NMDS_jaccard)
 
#colours

gen1<-ifelse(grepl("^S-plant$",group),"red",
             ifelse(grepl("*T34",group),"gray48",
                    ifelse(grepl("*control",group),"darkgoldenrod2","skyblue2")))
#shapes

pchvec<-ifelse(grepl("^S-plant$",group),"21",
               ifelse(grepl("*T34",group),"22",
                      ifelse(grepl("*control",group),"23","24")))


groupp<-as.data.frame(group)

#plot
pdf("NMDS_bray.pdf",width = 9,height = 8)
plot(NMDS_bray$points, cex=1.7,pch=21, bg=gen1,
     col="black", xlab="NMDS 1",ylab="NMDS 2")
ordihull(NMDS_bray,groups=groupp$group,draw="polygon",col="grey90",
         label=FALSE)
abline(h = 0, lty = "dotted")
abline(v = 0, lty = "dotted")
legend("topright", "Stress = 0.03",cex=0.8, bty="n")
legend("bottomright", 
       legend = c("S-plant","S-plant-T34","S-control","inlet"), 
       pt.bg= c("red","gray48","darkgoldenrod2","skyblue2"),
       col = "black",
       pch=21,
       bty = "n", 
       pt.cex = 1.6, 
       cex = 0.65, 
       text.col = "black",
       border = "black",
       horiz = F)
dev.off()

#cca
#matrix with taxa
taxa<-abundances[c(1:8),]
rownames(taxa)<-taxa$gene
taxa<-taxa[,-c(1:3)]
taxa1<-apply(taxa, 2, function(x) gsub("^$|^ $", 0, x))
taxa2<-t(taxa1)
taxa2<-apply(as.matrix(taxa2), 2, as.numeric)
rownames(taxa2)<-rownames(alpha3)
taxa2<-taxa2[,colSums(taxa2)!=0]
scripsc
#convert to dataframe args table
alpha5<-as.data.frame(alpha3)
alpha5<-alpha5[,colSums(alpha5)!=0]
#cca
vare.cca3<-cca(alpha5~taxa2)

###check for colineatrity between constrained variables(ABs)
vif.cca(vare.cca3)

#plot
pdf("CCA.pdf", width =9 , height = 7)
plot(vare.cca3)
dev.off()

#change to long format table

#adjust pvalue
taxa3<-as.data.frame(cbind(group,taxa2))

#to long format
taxa_long <- melt(setDT(taxa3), id.vars = c("group"), variable.name = "taxa")
taxa_long$value<-as.numeric(taxa_long$value)

stat.test <- taxa_long %>%
  group_by(taxa) %>%
  wilcox_test(value ~ group) %>%
  adjust_pvalue(method = "fdr") %>%
  add_significance()


write.csv(stat.test,"~/Desktop/dorde/results/taxa_abundance_signficance.csv")

#PCA args

XB<- log10(alpha5+1)
let<-cbind(group, XB)

pca<-prcomp(let[,-which(names(let) %in% c("group"))], scale=TRUE)


###eigenvalues, cumulative variance with the PCs...
eig.val <- get_eigenvalue(pca)
eig.val

write.table(eig.val,"eigenvalues_pca_args.txt")
##scree plots
##cutoff for 70% variance presented by the PCA
cumpro <- cumsum(pca$sdev^2 / sum(pca$sdev^2))

pdf("variance_dimentions.pdf", width = 7, height = 6)
plot(cumpro[0:8], xlab = "PC #", ylab = "Amount of explained variance", main = "Cumulative variance plot", type="b")
abline(h = 0.7, col="blue", lty=5)
legend("topleft", legend=c("Cut-off 70 %"),
       col=c("blue"), lty=5, cex=0.6)
dev.off()


#pca plot

###PCA plot

# factor for colors 
PCAcolors <- gen1

###scores and loadings
PCAscores <- pca$x
PCAloadings <- pca$rotation

####multiply loadings  coordinates by x since is hard to distinguish them apart     
l.x <- PCAloadings[,1]*13
l.y <- PCAloadings[,2]*13

# Label position
l.pos <- l.y # Create a vector of y axis coordinates
lo <- which(l.y < 0) # Get the variables on the bottom half of the plot
hi <- which(l.y > 0) # Get variables on the top half
# Replace values in the vector
l.pos <- replace(l.pos, lo, "1")
l.pos <- replace(l.pos, hi, "3")

###to add variance atomaticly to the xlab and ylab
s <- summary(pca)

let$group<-as.factor(let$group)

pdf("PCA_args.pdf", width = 10, height = 9)
plot(PCAscores[,1:2],  
     pch = 21,           
     col="black",                                  
     bg=PCAcolors,                               
     cex=1.7,                                     
     main="PCA",                                    
     xlab=paste("PCA 1 (", round(s$importance[2]*100, 1), "%)", sep = ""),
     ylab=paste("PCA 2 (", round(s$importance[5]*100, 1), "%)", sep = "")
)
legend("bottomright",                                
       legend=levels(let$group),             
       pch=21,                                    
       pt.bg= c("skyblue2", "darkgoldenrod2","red", "gray48"),                           
       pt.cex=1.5,                                
       col = "black",
       bty = "n"
)
abline(h = 0, lty = "dotted", col="black")
abline(v = 0, lty = "dotted", col="black")
###plot loadings and Draw arrows
arrows(x0=0, x1=l.x, y0=0, y1=l.y, col="red", length=0.09, lwd=2)
text(l.x,l.y, rownames(PCAloadings), col="blue",  pos=l.pos )
dev.off()
