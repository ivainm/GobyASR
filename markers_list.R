setwd("C://Users//ivainm//AppData//Local//Packages//CanonicalGroupLimited.Ubuntu18.04LTS_79rhkp1fndgsc//LocalState//rootfs//home//ivainm//realdir//Radsex_files")

install.packages("devtools")
devtools::install_github("SexGenomicsToolkit/sgtr")

library(sgtr)
library(ggplot2)
library(ggExtra)


# Get the list of Radsex markers associated with sex 
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#Haplotype 1 of the reference genome
radmaps<-read.table("mapsex_hap1.tsv", header=T)
radhits<-radmaps[radmaps$Signif=="True",]
write_tsv(radhits,"radsex_hap1.tsv")

#Haplotype 2 of the reference genome
radmaps2<-read.table("mapsex_hap2.tsv", header=T)
radhits2<-radmaps2[radmaps2$Signif=="True",]
write_tsv(radhits2,"radsex_hap2.tsv")


# Get the list of SNPs from Stacks where allele frequency is associated with sex 
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#Haplotype 1-------------------------------------------------
assochap1<-read.table("assocsex_hap1.assoc.fisher", header =T)
str(assochap1)
subass<-assochap1[assochap1$P<0.001,] #avoid plotting everything, subset by p value
str(subass)

ggplot(data=subass[subass$P<0.0000000001,], aes(x=F_A, y=F_U))+ #find a satisfying p value threshold using the plot
  geom_point( aes(fill=P),shape=21, size=3)+
  labs(x="Minor allele frequency in males",y="Minor allele frequency in females",title="Association of allele frequency with sex")+
  scale_fill_gradientn(colours = c("darkblue", "blue", "lightblue", "red"),
                       values = c(0,0.001, 0.05, 0.1,  1))

#Haplotype 2-------------------------------------------------
assochap2<-read.table("assocsex_hap2.assoc.fisher", header =T)
str(assochap2)
subass2<-assochap2[assochap2$P<0.001,]
str(subass2)

ggplot(data=subass2[subass2$P<0.0000000001,], aes(x=F_A, y=F_U))+ #use same p value threshold as for hap1
  geom_point( aes(fill=P),shape=21, size=3)+
  labs(x="Minor allele frequency in males",y="Minor allele frequency in females",title="Association of allele frequency with sex")+
  scale_fill_gradientn(colours = c("darkblue", "blue", "lightblue", "red"),
                       values = c(0,0.001, 0.05, 0.1,  1))


#Using the p value threshold defined graphically above, export list of markers
candfreq_hap1_m<-subass[subass$P<0.0000000001 &  #markers polymorphic in males, haplotype 1
                          subass$F_A>0.3,]
candfreq_hap1_f<-subass[subass$P<0.0000000001 &  #markers polymorphic in females, hap1
                          subass$F_U>0.5,]
candfreq_hap2_m<-subass2[subass2$P<0.0000000001 & #markers polymorphic in males, haplotype 2
                           subass2$F_A>0.3,]
candfreq_hap2_f<-subass2[subass2$P<0.0000000001 & #markers polymorphic in females, hap2
                           subass2$F_U>0.5,]

write_tsv(candfreq_hap1_m,"stacks_freq_hap1_m.txt")
write_tsv(candfreq_hap2_m,"stacks_freq_hap2_m.txt")


# Get the list of SNPs from Stacks where SNP missingness is associated with sex 
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#Haplotype 1 ------------------------------------------------
missinghap1<-read.table("missingsex_hap1.missing", header =T)
submiss<-missinghap1[missinghap1$P<0.05,]
str(submiss)

ggplot(data=submiss[submiss$P<0.000001,], aes(x=F_MISS_A, y=F_MISS_U))+ #define p value threshold graphically
  geom_point( aes(fill=P),shape=21, size=3)+
  labs(x="Missing call frequency in males",y="Missing call frequency in females",title="Association of locus missingness with sex")+
  scale_fill_gradientn(colours = c("darkblue", "blue", "lightblue", "red"),
                       values = c(0,0.001, 0.05, 0.1,  1))
#Haplotype 2 ------------------------------------------------
missinghap2<-read.table("missingsex_hap2.missing", header =T)
str(missinghap2)
submiss2<-missinghap2[missinghap2$P<0.05,]
str(submiss2)


ggplot(data=submiss2[submiss2$P<0.000001,], aes(x=F_MISS_A, y=F_MISS_U))+ #use same p value threshold 
  geom_point( aes(fill=P),shape=21, size=3)+
  labs(x="Missing call frequency in males",y="Missing call frequency in females",title="Association of locus missingness with sex")+
  scale_fill_gradientn(colours = c("darkblue", "blue", "lightblue", "red"),
                       values = c(0,0.001, 0.05, 0.1,  1))



#Subset data to extract position of SNPs of interest and write file
candmiss_hap1_m<-submiss[submiss$P<0.000001 &
                           submiss$F_MISS_U>0.25,]
candmiss_hap1_f<-submiss[submiss$P<0.000001 &
                           submiss$F_MISS_A>0.4,]
candmiss_hap2_m<-submiss2[submiss2$P<0.000001 &
                            submiss2$F_MISS_U>0.25,]
candmiss_hap2_f<-submiss2[submiss2$P<0.000001 &
                            submiss2$F_MISS_A>0.4,]


write_tsv(candmiss_hap1_m,"stacks_miss_hap1_m.txt")
write_tsv(candmiss_hap1_f,"stacks_miss_hap1_f.txt")
write_tsv(candmiss_hap2_m,"stacks_miss_hap2_m.txt")
write_tsv(candmiss_hap2_f,"stacks_miss_hap2_f.txt")
