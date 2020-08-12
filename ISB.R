#code adapted from 
#https://bioconductor.org/packages/release/bioc/vignettes/PCAtools/inst/doc/PCAtools.html#load-the-package-into-r-session

#install packages (only need to do this once)
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
#BiocManager::install("AnnotationDbi")
#BiocManager::install("illuminaHumanv4.db")

#call packages
library(PCAtools)
library(Biobase)
library(ggplot2)

#read files
setwd("/Users/anyas21/Desktop/RFiles")
pre_temporal_cortex <- read.csv("Temporal_Cortex.csv")
pre_cerebellum <- read.csv("Cerebellum.csv")
pre_frontal_cortex <- read.csv("Frontal_Cortex.csv")
pre_entorhinal_cortex <- read.csv("Entorhinal_Cortex.csv")

#TEMPORAL CORTEX
#modify the data
temporal_cortex <- pre_temporal_cortex[c(-1,-2,-3,-4), ] #create new data frame without first 4 rows
rownames(temporal_cortex) <- temporal_cortex[, 1] # set row names
colnames(temporal_cortex) <- temporal_cortex[1, ] # set column names
temporal_cortex <- temporal_cortex[-1, -1] #take out redundant first row and column

#create metadata table
tc_metadata_information <- pre_temporal_cortex[c(1,2,3,4,5),] #create table of metadata information
tc_metadata <- as.data.frame(t(tc_metadata_information)) #swap rows and columns
colnames(tc_metadata) <- c("Gender", "Age", "Individual", "Disease State") #add column names
tc_metadata <- tc_metadata[-1, ] #delete first row
rownames(tc_metadata) <- tc_metadata[, 5] #move ID to row name
tc_metadata <- tc_metadata[, -5] #delete fifth column
all(colnames(temporal_cortex) == rownames(tc_metadata)) #check to make sure no extra data has been deleted
temporal_cortex <- data.frame(apply(temporal_cortex, 2, function(x) as.numeric(as.character(x))), row.names = rownames(temporal_cortex)) #changing data frame type from character to numeric
tc_p <- pca(temporal_cortex, metadata = tc_metadata)

#plot data!
#screeplots
screeplot(tc_p) #create basic screeplot
#create a screeplot with Elbow Method and Horn's Analysis (for number of PCA) clearly labeled 
tc_horn <- parallelPCA(temporal_cortex)
tc_elbow <- findElbowPoint(tc_p$variance)
screeplot(tc_p,
          components = getComponents(tc_p, 1:20),
          vline = c(tc_horn$n, tc_elbow)) +
  geom_text(aes(tc_horn$n + 1, 50, label = "Horn's", vjust = -1)) +
  geom_text(aes(tc_elbow + 1, 50, label = "Elbow", vjust = -1))
#biplots
biplot(tc_p) #create basic biplot
#create a biplot with those with Alzheimer's in a different color
biplot(tc_p, 
       colby = 'Disease State', 
       hline = 0, vline = 0, 
       legendPosition = 'right',
       shape = 'Gender', 
       shapekey = c(' MALE'=15, ' FEMALE'=17), 
       labSize = 2.0,
       title = 'PCA bi-plot',
       subtitle = 'PC1 versus PC2') 
#loadings plots
dev.new(width=10, height=4, unit="in") #make the plot viewing window larger to avoid error in RStudio
plotloadings(tc_p) #create basic loadings plot
#create a loadings plot with only the top 1% of variables
dev.new(width=10, height=4, unit="in") #make the plot viewing window larger to avoid error in RStudio
plotloadings(tc_p,
             rangeRetain = 0.01,
             labSize = 2.0,
             title = 'PCA loadings plot',
             subtitle = 'Top 1% of variables',
             shape = 24,
             drawConnectors = TRUE)

#CEREBELLUM
#modify the data
cerebellum <- pre_cerebellum[c(-1,-2,-3,-4), ] #create new data frame without first 4 rows
rownames(cerebellum) <- cerebellum[, 1] # set row names
colnames(cerebellum) <- cerebellum[1, ] # set column names
cerebellum <- cerebellum[-1, -1] #take out redundant first row and column

#create metadata table
c_metadata_information <- pre_cerebellum[c(1,2,3,4,5),] #create table of metadata information
c_metadata <- as.data.frame(t(c_metadata_information)) #swap rows and columns
colnames(c_metadata) <- c("Gender", "Age", "Individual", "Disease State") #add column names
c_metadata <- c_metadata[-1, ] #delete first row
rownames(c_metadata) <- c_metadata[, 5] #move ID to row name
c_metadata <- c_metadata[, -5] #delete fifth row
all(colnames(cerebellum) == rownames(c_metadata)) #check to make sure no extra data has been deleted
cerebellum <- data.frame(apply(cerebellum, 2, function(x) as.numeric(as.character(x))), row.names = rownames(temporal_cortex)) #changing data frame type from character to numeric
c_p <- pca(cerebellum, metadata = c_metadata)

#plot data!
#screeplots
screeplot(c_p) #create basic screeplot
#create a screeplot with Elbow Method and Horn's Analysis (for number of PCA) clearly labeled 
c_horn <- parallelPCA(cerebellum)
c_elbow <- findElbowPoint(c_p$variance)
screeplot(c_p,
          components = getComponents(c_p, 1:25),
          vline = c(c_horn$n, c_elbow)) +
  geom_text(aes(c_horn$n + 1, 50, label = "Horn's", vjust = -1)) +
  geom_text(aes(c_elbow + 1, 50, label = "Elbow", vjust = -1))
#biplots
biplot(c_p) #create basic biplot
#create a biplot with those with Alzheimer's in a different color
biplot(c_p, 
       colby = 'Disease State', 
       hline = 0, vline = 0, 
       legendPosition = 'right',
       shape = 'Gender', 
       shapekey = c(' MALE'=15, ' FEMALE'=17), 
       labSize = 2.0,
       title = 'PCA bi-plot',
       subtitle = 'PC1 versus PC2') 
#loadings plots
dev.new(width=10, height=4, unit="in") #make the plot viewing window larger to avoid error in RStudio
plotloadings(c_p) #create basic loadings plot
#create a loadings plot with only the top 1% of variables
dev.new(width=10, height=4, unit="in") #make the plot viewing window larger to avoid error in RStudio
plotloadings(c_p,
             rangeRetain = 0.01,
             labSize = 2.0,
             title = 'PCA loadings plot',
             subtitle = 'Top 1% of variables',
             shape = 24,
             drawConnectors = TRUE)


#FRONTAL CORTEX
#modify the data
frontal_cortex <- pre_frontal_cortex[c(-1,-2,-3,-4), ] #create new data frame without first 4 rows
rownames(frontal_cortex) <- frontal_cortex[, 1] # set row names
colnames(frontal_cortex) <- frontal_cortex[1, ] # set column names
frontal_cortex <- frontal_cortex[-1, -1] #take out redundant first row and column

#create metadata table
fc_metadata_information <- pre_frontal_cortex[c(1,2,3,4,5),] #create table of metadata information
fc_metadata <- as.data.frame(t(fc_metadata_information)) #swap rows and columns
colnames(fc_metadata) <- c("Gender", "Age", "Individual", "Disease State") #add column names
fc_metadata <- fc_metadata[-1, ] #delete first row
rownames(fc_metadata) <- fc_metadata[, 5] #move ID to row name
fc_metadata <- fc_metadata[, -5] #delete fifth row
all(colnames(frontal_cortex) == rownames(fc_metadata)) #check to make sure no extra data has been deleted
frontal_cortex <- data.frame(apply(frontal_cortex, 2, function(x) as.numeric(as.character(x))), row.names = rownames(temporal_cortex)) #changing data frame type from character to numeric
fc_p <- pca(frontal_cortex, metadata = fc_metadata)

#plot data!
#screeplots
screeplot(fc_p) #create basic screeplot
#create a screeplot with Elbow Method and Horn's Analysis (for number of PCA) clearly labeled 
fc_horn <- parallelPCA(frontal_cortex)
fc_elbow <- findElbowPoint(fc_p$variance)
screeplot(fc_p,
          components = getComponents(fc_p, 1:25),
          vline = c(fc_horn$n, fc_elbow)) +
  geom_text(aes(fc_horn$n + 1, 50, label = "Horn's", vjust = -1)) +
  geom_text(aes(fc_elbow + 1, 50, label = "Elbow", vjust = -1))
#biplots
biplot(fc_p) #create basic biplot
#create a biplot with those with Alzheimer's in a different color
biplot(fc_p, 
       colby = 'Disease State', 
       hline = 0, vline = 0, 
       legendPosition = 'right',
       shape = 'Gender', 
       shapekey = c(' MALE'=15, ' FEMALE'=17), 
       labSize = 2.0,
       title = 'PCA bi-plot',
       subtitle = 'PC1 versus PC2') 
#loading plots
dev.new(width=10, height=4, unit="in") #make the plot viewing window larger to avoid error in RStudio
plotloadings(fc_p) #create basic loadings plot
#create a loadings plot with only the top 1% of variables
dev.new(width=10, height=4, unit="in") #make the plot viewing window larger to avoid error in RStudio
plotloadings(fc_p,
             rangeRetain = 0.01,
             labSize = 2.0,
             title = 'PCA loadings plot',
             subtitle = 'Top 1% of variables',
             shape = 24,
             drawConnectors = TRUE)

#ENTORHINAL CORTEX
#modify the data
entorhinal_cortex <- pre_entorhinal_cortex[c(-1,-2,-3,-4), ] #create new data frame without first 4 rows
rownames(entorhinal_cortex) <- entorhinal_cortex[, 1] # set row names
colnames(entorhinal_cortex) <- entorhinal_cortex[1, ] # set column names
entorhinal_cortex <- entorhinal_cortex[-1, -1] #take out redundant first row and column

#create metadata table
ec_metadata_information <- pre_entorhinal_cortex[c(1,2,3,4,5),] #create table of metadata information
ec_metadata <- as.data.frame(t(ec_metadata_information)) #swap rows and columns
colnames(ec_metadata) <- c("Gender", "Age", "Individual", "Disease State") #add column names
ec_metadata <- ec_metadata[-1, ] #delete first row
rownames(ec_metadata) <- ec_metadata[, 5] #move ID to row name
ec_metadata <- ec_metadata[, -5] #delete fifth row
all(colnames(entorhinal_cortex) == rownames(ec_metadata)) #check to make sure no extra data has been deleted
entorhinal_cortex <- data.frame(apply(entorhinal_cortex, 2, function(x) as.numeric(as.character(x))), row.names = rownames(temporal_cortex)) #changing data frame type from character to numeric
ec_p <- pca(entorhinal_cortex, metadata = ec_metadata)

#plot data!
#screeplots
screeplot(ec_p) #create basic screeplot
#create a screeplot with Elbow Method and Horn's Analysis (for number of PCA) clearly labeled 
ec_horn <- parallelPCA(entorhinal_cortex)
ec_elbow <- findElbowPoint(ec_p$variance)
screeplot(ec_p,
          components = getComponents(fc_p, 1:25),
          vline = c(ec_horn$n, ec_elbow)) +
  geom_text(aes(ec_horn$n + 1, 50, label = "Horn's", vjust = -1)) +
  geom_text(aes(ec_elbow + 1, 50, label = "Elbow", vjust = -1))
#biplots
biplot(ec_p) #create basic biplot
#create a biplot with those with Alzheimer's in a different color
biplot(ec_p, 
       colby = 'Disease State', 
       hline = 0, vline = 0, 
       legendPosition = 'right', 
       shape = 'Gender', 
       shapekey = c(' MALE'=15, ' FEMALE'=17), 
       labSize = 2.0,
       title = 'PCA bi-plot',
       subtitle = 'PC1 versus PC2') 
#loadings plots
dev.new(width=10, height=4, unit="in") #make the plot viewing window larger to avoid error in RStudio
plotloadings(ec_p) #create basic loadings plot
#create a loadings plot with only the top 1% of variables
dev.new(width=10, height=4, unit="in") #make the plot viewing window larger to avoid error in RStudio
plotloadings(ec_p,
             rangeRetain = 0.01,
             labSize = 2.0,
             title = 'PCA loadings plot',
             subtitle = 'Top 1% of variables',
             shape = 24,
             drawConnectors = TRUE)

#read in data
setwd("/Users/anyas21/Desktop/RFiles")
data<-read.csv("ISB_data_cerebellum_train.csv")
#gene_names<-read.csv("Gene_Names.csv")

#create a dataframe with illumina ids as gene names
#avector <- gene_names[,1]
#new_names <- data.frame(Gene=unlist(mget(x = avector,envir = illuminaHumanv4SYMBOL)))
