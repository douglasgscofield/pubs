#Ioana's script for calculating Similarity Indexes, NMDS

#reset R's brain (clear the data)
rm (list=ls())

#get the path
getwd ()

#set working directory
setwd ("/Users/ioana/Desktop/LICHEN_GUILD_PROJECT2016/OTU_clustering/clustering")

library (vegan)
library(ggplot2)


## step 1 rarefing the data
#make my tabel an object in R, I have a header and reads the first rows
OTUs_table<-read.csv("otutab.edgar_20160602_first_200_1.5.csv", header = T, row.names = 1,  sep= ",") 

OTUs_table_transp<-t(OTUs_table) ##transpose data set by swaping columns to raws

dim(OTUs_table_transp) #tells me how may rows and how many columns I have in the table

OTUs_table_transp<-OTUs_table_transp[1:6,] ##Removes negative controll column

otuRarefied<-rrarefy(OTUs_table_transp, min(rowSums(OTUs_table_transp))) ##rarefy to smallest integer

#rd4 = rd4[ ,colSums(rd4)!=0] ## it will eliminate the columns that have only 0s
#na.omit() ##removes the columns with NA 

write.table(otuRarefied, "oturrafey1.txt", sep = "\t") #writes the rarefied table
#write.table(OTUs_table_transp, "otus_Table_transp.txt", sep = "\t")


##step 2: calculating the Bray-Curtis distance index
#attention: the metaMDS function already calculates the B-C index, so I think I will the "otuRarefied" matrix and not "D_bray" matrix
#D_bray_Rarefied<-vegdist(otuRarefied, distance = "bray") #Bray-Curtis
#D_bray_nonRarefied<-vegdist(OTUs_table_transp, distance = "bray") #Bray-Curtis

############ chenck if i have enough dimenssion
#nmds_otuRarefied.mds1<-metaMDS(otuRarefied,k=1,autotransform=F,trymax=50)
#nmds_otuRarefied.mds2<-metaMDS(otuRarefied,k=2,autotransform=F,trymax=50)
#nmds_otuRarefied.mds3<-metaMDS(otuRarefied,k=3,autotransform=F,trymax=50)
#nmds_otuRarefied.mds4<-metaMDS(otuRarefied,k=4,autotransform=F,trymax=50)
#nmds_otuRarefied.mds5<-metaMDS(otuRarefied,k=5,autotransform=F,trymax=50)
#and now we can plot the stress value against the number of dimensions
#stress.val<-c(nmds_otuRarefied.mds1$stress,nmds_otuRarefied.mds2$stress,nmds_otuRarefied.mds3$stress,nmds_otuRarefied.mds4$stress,nmds_otuRarefied.mds5$stress)
#plot(1:5,stress.val,type='b')


##step 3: make the NMDS plots with the vegan ploting
#function used: metaMDS()
nmds_otuRarefied<-metaMDS(otuRarefied, k=2, trymax=2000)
#nmds_otuNonRarefied<-metaMDS(OTUs_table_transp, k=3, trymax=2000)


nmds_otuRarefied$stress
nmds_otuRarefied_scores<-scores(nmds_otuRarefied)
write.table(nmds_otuRarefied_scores, "nmds_otuRarefied_scores.txt", sep = "\t") #writes the NMDS scores table


#look at the Shepard plot: shows scatter around the regression between the interpoint distances in the final configuration 
#(i.e., the distances between each pair of communities) against their original dissimilarities. larger the scatters around the stress line (red), worse is the fit
stressplot(nmds_otuRarefied)


#plot the nmds with the vegan package
#plot(nmds_otuRarefied) #this is just very simple
#ordiplot(nmds_otuRarefied,type="n")
#orditorp(nmds_otuRarefied,display="species",col="red", air=0.01, cex=1)
#orditorp(nmds_otuRarefied,display="sites", col=c(rep("royalblue1",2),rep("gold2",2), rep("red4",2)),cex=1.25,air=0.01)
#legend("bottomright", pch = 19, cex=1.5, legend = c("Gardabaer", "Oland 1", "Oland 2"), col=c("royalblue1", "gold2","red4",0.01))



#or Use ggplot for the NMDS plot
#make my dataset
data.scores <- as.data.frame(scores(nmds_otuRarefied))  #Using the scores function from vegan to extract the site scores and convert to a data.frame
data.scores$site <- rownames(data.scores)  # create a column of site names, from the rownames of data.scores
data.scores[,"locality"]  <- c(rep("Gardabaer",2), rep("Oland1",2), rep("Oland2",2)) #added the column locality; rep = repeat Oland1 twice)
data.scores[,"lichen"]  <- c("Thamnolia Gardabaer", "Cetraria Gardabaer","Thamnolia Oland1", "Cetraria Oland1", "Thamnolia Oland2", "Cetraria Oland2")
data.scores #visualize data.score

species.scores <- as.data.frame(scores(nmds_otuRarefied, "species"))  #Using the scores function from vegan to extract the species scores and convert to a data.frame
species.scores$species <- rownames(species.scores)  # create a column of species, from the rownames of species.scores
Photobiont<- c("Trebouxia simplex C3", "Trebouxia simplex C3", "Trebouxia simplex C3", "Trebouxia simplex C3", "Trebouxia simplex C3", "Trebouxia simplex B", "Trebouxia simplex C3", "Trebouxia simplex B","Trebouxia simplex C3", "Trebouxia simplex B", "Trebouxia simplex C1", "Trebouxia simplex C3", "Trebouxia simplex B", "Trebouxia simplex B", "Trebouxia simplex B", "Trebouxia simplex B", "Trebouxia simplex B", "Trebouxia simplex B", "OTU25", "Trebouxia simplex B", "Trebouxia simplex B", "Trebouxia simplex B", "Trebouxia jamesii-vulpinae", "Trebouxia simplex C1", "Trebouxia simplex C1", "Trebouxia simplex C1", "Trebouxia simplex C1", "Trebouxia simplex C1", "Trebouxia gigantea-vagua", "Trebouxia simplex C1", "Trebouxia gigantea-vagua", "Trebouxia simplex C1")
species.scores$Photobiont <-Photobiont
species.scores  #look at the data

#make # hull values for groups
grp.TsB <- species.scores[species.scores$Photobiont == "Trebouxia simplex B", ][chull(species.scores[species.scores$Photobiont == 
                                                                                            "Trebouxia simplex B", c("NMDS1", "NMDS2")]), ]  
grp.TsC1 <- species.scores[species.scores$Photobiont == "Trebouxia simplex C1", ][chull(species.scores[species.scores$Photobiont == 
                                                                                              "Trebouxia simplex C1", c("NMDS1", "NMDS2")]), ]
grp.TsC3 <- species.scores[species.scores$Photobiont == "Trebouxia simplex C3", ][chull(species.scores[species.scores$Photobiont == 
                                                                                              "Trebouxia simplex C3", c("NMDS1", "NMDS2")]), ]
grp.Tgv <- species.scores[species.scores$Photobiont == "Trebouxia gigantea-vagua", ][chull(species.scores[species.scores$Photobiont == 
                                                                                          "Trebouxia gigantea-vagua", c("NMDS1", "NMDS2")]), ]
grp.Tjv <- species.scores[species.scores$Photobiont == "Trebouxia jamesii-vulpinae", ][chull(species.scores[species.scores$Photobiont == 
                                                                                          "Trebouxia jamesii-vulpinae", c("NMDS1", "NMDS2")]), ]
grp.OTU25 <- species.scores[species.scores$Photobiont == "OTU25", ][chull(species.scores[species.scores$Photobiont == 
                                                                                          "OTU25", c("NMDS1", "NMDS2")]), ]
hull.data <- rbind(grp.TsB, grp.TsC1, grp.TsC3, grp.Tgv, grp.Tjv, grp.OTU25)  #combine grp.a b and c
hull.data

ggplot() + 
#geom_text(data=data.scores,aes(x=NMDS1,y=NMDS2,label=lichen),size=5) +  # add the site labels
geom_polygon(data=hull.data,aes(x=NMDS1,y=NMDS2, group=Photobiont),alpha=0.2) + # add the convex hulls
geom_text(data=species.scores,aes(x=NMDS1,y=NMDS2,label=NA),alpha=0.5) +  # add the species labels  
geom_point(data=species.scores,aes(x=NMDS1,y=NMDS2,shape=Photobiont,colour=Photobiont),size=6) + # add the point markers
scale_shape_manual(values=c(20, 17, 17, 17, 17, 17))+
#scale_size_manual(values=c(2,4,6,8,10,12))+
scale_color_manual(values=c("Trebouxia simplex B" = "#ccff00", "Trebouxia simplex C1" = "#7137c8", "Trebouxia simplex C3" = "#ccaaff", "Trebouxia gigantea-vagua" = "#ff0066", "Trebouxia jamesii-vulpinae" = "#d4aa00", "OTU25" = "#00ffcc")) +
scale_x_continuous (limits = c(-0.6, 0.8))+
scale_y_continuous (limits = c(-0.6, 0.5))+
coord_equal()

ggplot() + 
geom_point(data=data.scores,aes(x=NMDS1,y=NMDS2,shape=lichen,colour=lichen), size=9)+ # add the point markers
#geom_text(data=data.scores,aes(x=NMDS1,y=NMDS2,label=lichen),size=3) +  # add the site labels
scale_shape_manual(values=c(18, 18, 18, 20, 20, 20))+
  #scale_size_manual(values=c(2,4,6,8,10,12))+
scale_color_manual(values=c("Thamnolia Gardabaer" = "#0066ff", "Cetraria Gardabaer" = "#0066ff", "Thamnolia Oland1" = "#ffcc00", "Cetraria Oland1" = "#ffcc00","Thamnolia Oland2" = "#AA0000", "Cetraria Oland2" = "#AA0000" )) +
scale_x_continuous (limits = c(-0.6, 0.8))+
scale_y_continuous (limits = c(-0.6, 0.5))+
coord_equal()

