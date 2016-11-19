
#Extracción de los nombres de cada línea de la conversación:
WhatsApp <- readLines("WhatsApp.txt", encoding = "UTF-8") #encoding = "UCS-2LE"

#########
#several files
WhatsApp <- c(readLines("WhatsApp.txt", encoding = "UTF-8"),
              readLines("WhatsApp2.txt", encoding = "UTF-8"),
              readLines("WhatsApp3.txt", encoding = "UTF-8"))
#########

wlines <- WhatsApp #extract lines

#extract names from lines
for (i in 1:length(WhatsApp)){
  #wlines[i] <- gsub("[0-9]+ (.*) [0-9]{2}:[0-9]{2} - ","", WhatsApp[i])
  wlines[i] <- gsub(": (.*)","", wlines[i])
  wlines[i] <- gsub("(.*) - ","", wlines[i])
  #wlines[i] <- gsub(" [A-z]*(.*)","",wlines[i])
  wlines[i] <- substr(wlines[i], 1, 10) #first characters
}

#list of names

table <- as.data.frame(table(wlines))
names <- table$wlines[table$Freq>length(WhatsApp)/300 & table$wlines!=""]

sellines <- wlines %in% names

nlines <- wlines[sellines]
as.data.frame(table(nlines))

#############

#run this if there are people with different names in each file
index1 <- 6
index2 <- 7
nlines[nlines==names[index1]] <- as.character(names[index2])
table <- as.data.frame(table(nlines))
names <- table$nlines
as.data.frame(table(nlines))
#############

write.csv(nlines, file = "nlines.csv")


#erase consecutive repeated names
collapselines = nlines
i=1
while (i < length(collapselines)){
  
  if (collapselines[i]==collapselines[i+1]){
    
    collapselines <- collapselines[-c(i)]
    
  }  else i=i+1
  
}
#nlines <-collapselines

#if repeated names are not deleted
#collapselines <- nlines

library(RColorBrewer)
colorpalette <- colorRampPalette(c("white","black"))(256)

#buckets

## bucket range - can be modified
r <- round(n/3)


B <- matrix(F, nrow=length(collapselines)-r+1, ncol=length(names))

for (i in 1:(length(collapselines)-r+1)){
  for (j in 1:length(names)){
    B[i,j] <- names[j] %in% collapselines[i:(i+r-1)]
  }
}


### affinity (correlation) clusters

CorrB <- cor(B)

#rownames(CorrB) <- names
#colnames(CorrB) <- names


diag(CorrB) <- NA
heatmap(1-CorrB, col=colorpalette, labRow=names, labCol=names)
hc <- hclust(as.dist(1-CorrB))
plot(hc, labels=names, cex=0.7)


ngroups <- 4
Mgroups <- cutree(hc, k=ngroups)


distCorr <- dist(cor(B))
hc <- hclust(distCorr)
plot(hc, labels=names, cex=0.7)
Mgroups <- cutree(hc, k=ngroups)

distCorrm <- as.matrix(distCorr)
diag(distCorrm)<-NA

heatmap(distCorrm, col=colorpalette, labRow=names, labCol=names)

cmdC <- cmdscale(distCorr, k = 18, eig=T, add=T)
km.out <- kmeans(cmdC$points,centers=ngroups,nstart=15)
km.out$cluster
indexcoords <- c(1,2)
plot(cmdC$points[,indexcoords])
#color tree clusters
text(cmdC$points[,indexcoords],labels=names, pos=4, cex=0.6, col=Mgroups)
#color kmeans clusters
text(cmdC$points[,indexcoords],labels=names, pos=4, cex=0.6, col=km.out$cluster)


### simultaneity clusters


distB <- dist(t(B), method="binary") #proportion of buckets with shared presence
hc <- hclust(distB, method="complete")
plot(hc,cex=.7, labels=names) #dendogram for simultaneity in conversation

Mgroups <- cutree(hc, k=ngroups)

cmdB <- cmdscale(distB, k = 2, eig=T, add=T)

km.out <- kmeans(cmdB$points,centers=ngroups,nstart=15)
km.out$cluster

indexcoords <- c(1,2)
plot(cmdB$points[,indexcoords])
#color tree clusters
text(cmdB$points[,indexcoords],labels=names, pos=4, cex=0.6, col=Mgroups)
#color kmeans clusters
text(cmdB$points[,indexcoords],labels=names, pos=4, cex=0.6, col=km.out$cluster)

plot(cmdB$eig)

#same result when mapping the points
cmdB <- cmdscale(distB, k = length(names)-1, eig=T)
hc <- hclust(dist(cmdB$points), method="complete")
plot(hc,cex=.7, labels=names)
#

distBm <- as.matrix(distB)
diag(distBm)<-NA
heatmap(distBm, col=colorpalette, labRow=names, labCol=names)

