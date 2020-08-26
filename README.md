# Analyzing RNA-seq Alzheimer's Data
Tools used: Excel, R, Bioconductor. Make sure you have access to RStudio as well.
# About this project:
I first looked at this tutorial (https://bioconductor.org/packages/release/bioc/vignettes/PCAtools/inst/doc/PCAtools.html#load-the-package-into-r-session) from Bioconductor to read and analyze normalized RNA-sequencing data from a breast cancer study. I adapted the tutorial to fit this study (https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE144254), which analyzes the APOE2 differences in brain transcriptomics. For context, the APOE2 gene has been linked to Alzheimer's disease, which is the leading form of dementia. According to the Alzheimer's Association, Alzheimer's disease currently affects approximately 5.8 million individuals in the United States, two thirds of them women, and could affect as many as 13.8 million by 2050. To that end, finding a cure for Alzheimer's is incredibly important.
In contrast to the Bioconductor tutorial, I used Excel to modify my data before feeding it into RStudio; I then changed the tutorial code to better fit my environment. Then, I turned my process into my own tutorial for analyzing RNA-seq data from GEO using Excel, Bioconductor, and R. I've broken down this process into several steps, as seen below. 
# Downloading Data
It's possible - and sometimes more useful - to download the file directly into RStudio by using GEOquery and the GEO accession code. However, I found this approach doesn't always work for all kinds of data. If you're running into errors or not seeing the desired output, manually downloading and modifying the file may be best.
Steps:
1. Download the series matrix file from GEO. You can access that here. 
2. Unzip data and move it to an accessible folder. I used WinZip for this and created a folder, called "RFiles," specifically for this project.
# Modifying Data Using Excel
I chose to use Excel to quickly modify my data rather than utilizing RStudio for this part of the process. It's possible to run these modifications in R, but I found my computer slowed down quite a bit while doing so.  If you would like to run this in RStudio, simply hit "import dataset" in the environment tab and go from there. Note that if you're using the .csvs connected with this tutorial, you will have to combine temporal_cortex_1 and temporal_cortex_2. 
Steps:
1. Open Excel and create a new blank workbook.
2. Click on the "data" tab, then on "from text/csv."
3. Navigate to your unzipped file and open it.
4. After the data loads, click on "transform data."
5. Click on "split column," "by delimiter," and then select "tab" as your delimiter.
6. Delete unneeded rows from the top of the file. I deleted through row 36.
7. Click "close and load."
8. Delete extraneous rows from in the middle of your dataset (such as information about experimental procedure). I deleted rows 1-10 and 16-36. DO NOT delete the metadata (information about the age, gender, disease state, etc. of subjects) which you will analyze later. 
9. Change necessary names. For example, I changed a row name from "sample_characteristics_ch1" to "Gender." Other new row names included: "Age," "Disease State," and "Tissue." 
10. Some of the data says "tissue: temporal_cortex" instead of just "temporal_cortex." To fix this, use cmd+h on your keyboard and replace "tissue: " with "". Do this for the age, disease state, and gender rows as well. 
11. Search through the data for any NA values. Make sure not to delete any "entorhinal_cortex" data (since the name has "na" in it)! Then delete rows with NA values (the dataset is large enough that these individuals will not have a big impact).
12. Since there are four sections of the brain in this data, it's helpful to create a new .csv for each section of the brain. To do so, go to "data," "sort," "options," "sort left to right," and then choose the row that has your tissue data. Hit okay. Your data should now be sorted in alphabetical order by tissue type.
13. Save the original Excel file to folder "RFiles" as a .csv. 
14. Copy each tissue type to a seperate Excel spreadsheet and save to "RFiles" as a .csv. I named my files according to tissue type; i.e. "temporal_cortex.csv" or "cerebellum.csv". 
# Getting Set Up In RStudio
Before we import the data, we have a few packages to install and load beforehand. This is where the coding starts!
1. Install BiocManager, ggplot2, PCAtools, and Biobase. You only have to do this once.

'''
if (!requireNamespace("BiocManager", quietly = TRUE))
 	install.packages("BiocManager")
install.packages("ggplot2")
BiocManager::install("PCAtools")
BiocManager::install("Biobase")
'''

2. Call the necessary packages. You'll need to do this every time you restart RStudio.

'''
library(PCAtools)
library(Biobase)
library(ggplot2)
'''

# Reading and Modifying Data in RStudio
You should now have five .csv files: one big file with all the data and four smaller files divided up by part of the brain. I'll just demonstrate the next few lines of code on my temporal cortex files, named "temporal_cortex.csv". The code can be applicable to any section of the brain. In fact, I ran the code with all the .csv files (just remembering to change the variable names) in order to see differences between sections.

1. Get and set your working directory. Modify setwd() to correspond to the path to your folder.

'''
getwd()
setwd("/Users/anyas21/Desktop/RFiles")
'''

2. Read the .csv into RStudio and save it as a variable. This will create a data frame.

'''
pre_temporal_cortex <- read.csv("Temporal_Cortex.csv")
'''

Now modify your data! We want a data frame that just has the data collected from the experiment, not the metadata (which is in the first 4 rows of the data frame).

3. Create a new data frame without the first 4 rows of the original data frame.

'''
temporal_cortex <- pre_temporal_cortex[c(-1,-2,-3,-4), ]
'''

4. You'll see that the column and row titles aren't quite what we want them to be. Sub in gene names (which starts with "ILUMIN_") and subject ids (which start with "GSM") as the row titles and column titles, respectively.

'''
rownames(temporal_cortex) <- temporal_cortex[, 1]
colnames(temporal_cortex) <- temporal_cortex[1, ] 
'''

5. Delete the now redundant first row and column.

'''
temporal_cortex <- temporal_cortex[-1, -1]
'''

Now that we've finished modifying the data itself, we need to set up a metadata table.

6. Create a new variable from the metadata information (rows 1-4 on the original data frame, plus row 5 so we can see which subject ID matches with the metadata values).

'''
tc_metadata_information <- pre_temporal_cortex[c(1,2,3,4,5),]
'''

7. Swap the rows and columns of the metadata data frame.

'''
tc_metadata <- as.data.frame(t(tc_metadata_information))
'''

8. Change column/row titles to provide more information and delete some now redundant columns/rows. 

'''
colnames(tc_metadata) <- c("Gender", "Age", "Individual", "Disease State") #add column names
tc_metadata <- tc_metadata[-1, ] #delete first row
rownames(tc_metadata) <- tc_metadata[, 5] #move ID to row name
tc_metadata <- tc_metadata[, -5] #delete fifth column
'''

9. Accuracy check! Make sure the number of columns in your experimental data is equal to the number of rows in your metadata. If it's not, check your code to see where you went wrong.

'''
all(colnames(temporal_cortex) == rownames(tc_metadata))
'''

10. Convert the data frame with experimental data from a character object to a numeric object.

'''
temporal_cortex <- data.frame(apply(temporal_cortex, 2, function(x) as.numeric(as.character(x))), row.names = rownames(temporal_cortex))
'''

11. Create a new pca object with two arguments: the dataframe with experimental data and the metadata data frame.

'''
tc_p <- pca(temporal_cortex, metadata = tc_metadata)
'''

# Graphing Your Data!

This is the fun part! We'll run PCA analysis on this data. 
PCA is best used when dealing with large numbers of variables, which would have many "dimensions" when graphed. As we can't visualize much beyond 3D,  analysis reduces the "dimensionality" of the information without eliminating much of the data itself. In short, PCA analysis shows us what variables are important and which don't have much effect. PC1 describes the variables with the most variance, and PC2 describes the variables with the second-highest variance. 

Plot Type #1: Scree plots. Scree plots are used to demonstrate the variance associated with each PC. Generally, you don't want to use more PCs than absolutely necessary, because evaluating more components means more computing power. Scree plots help you understand how many PC you need to use.

Code for a basic scree plot:

'''
screeplot(tc_p)
'''

At some point, there is a group of PCs that account for the majority of the data, and incorporating new PCs won't tell us anything new about the data. There are several methods you can use to develop a better understanding of how many PCs to use. One is the Horn Analysis, another is the Elbow Method. You can run both on your scree plot. (See below.)

'''
tc_horn <- parallelPCA(temporal_cortex)
tc_elbow <- findElbowPoint(tc_p$variance)
'''

We've specified the components here, so the function will only graph PC 1 through 20. Note that the Elbow (which says to stop at PC8) and Horn's (which says to stop at PC17) analyses do not produce the same output! 

'''
screeplot(tc_p,
          components = getComponents(tc_p, 1:20),
          vline = c(tc_horn$n, tc_elbow)) +
  geom_text(aes(tc_horn$n + 1, 50, label = "Horn's", vjust = -1)) +
  geom_text(aes(tc_elbow + 1, 50, label = "Elbow", vjust = -1))
'''
  
Plot Type #2: Biplot. A biplot graphs data points on axes of PC1 and PC2, so you can see which one has more of an impact on the data. You can also add in metadata markers to learn more about metadata's impacts.

Here's a basic biplot from the data.

'''
biplot(tc_p) #create basic biplot
'''

The above biplot is a little hard to read. Let's modify the biplot with different colors for each disease state (top graph), as well as smaller fonts and  different shapes for male and females (bottom graph). Because males and females have different rates of Alzheimer's, I thought they would be a valuable group to look into. We can see there are clear groups of AD, AsymAD, and control subjects in the Temporal Cortex graph.

'''
biplot(tc_p, 
       colby = 'Disease State', 
       hline = 0, vline = 0, 
       legendPosition = 'right',
       shape = 'Gender', 
       shapekey = c(' MALE'=15, ' FEMALE'=17), 
       labSize = 2.0,
       title = 'PCA bi-plot',
       subtitle = 'PC1 versus PC2')
'''

Plot Type #3: Loadings Plot. A loadings plot gives us more insight into what exactly affects a PC, as in these loadings plots, we can see the genes that are associated with each. This will help us gain further understanding of which genes are responsible for variance.

We'll create a new window to avoid an RStudio error, as well as make a loadings plot.  Run both of these lines at the same time.

'''
dev.new(width=10, height=4, unit="in")
plotloadings(tc_p) #create basic biplot
'''

Now we'll edit the loadings plot to only show the top 1% of genes responsible for each PC. You'll still need to make a new window to run this function.

'''
dev.new(width=10, height=4, unit="in") 
plotloadings(tc_p,
             rangeRetain = 0.01,
             labSize = 2.0,
             title = 'PCA loadings plot',
             subtitle = 'Top 1% of variables',
             shape = 24,
             drawConnectors = TRUE)
'''
             
# Extra Step: Identifying Gene Names

In the loadings plot, each of the triangles corresponds to an Illumina ID... but not a gene name. If we want to find the gene names, we have to do some extra coding. 

1. Install and call illuminaHumanv4.db

'''
BiocManager::install("illuminaHumanv4.db")
library(illuminaHumanv4.db)
'''

2. Create a new variable with a group of all the variables retained from your loadings plot (you should see this output in the console). 

'''
tc_id<-c("ILMN_1676283", "ILMN_1732921", "ILMN_2218758", "ILMN_1670652", "ILMN_1839647", "ILMN_2207988", "ILMN_2255133", "ILMN_3245103", "ILMN_1736178", "ILMN_1653856", "ILMN_3248171", "ILMN_1679984", "ILMN_1691097", "ILMN_2172318", "ILMN_1788874")
'''

3. This line of code maps your inputted IDs to a gene name, then prints the resulting data frame. 

'''
data.frame(Gene=unlist(mget(x = tc_id,envir = illuminaHumanv4SYMBOL)))
'''

4. Since the entorhinal cortex shows the clearest clustering, let's find out more about the genes involved using genecards.org. 

KNCT1: Related to the creation of potassium channels. 
SERPINA3: Variations in this protein's sequence have been known to cause Alzheimer's.
CD44: Encodes a protein involved in cell-cell interactions.
CD24: Encodes a protein involved in cell growth and differentiation.
VCAM1: Encodes a protein involved in cell-cell interactions.
ZIC3: Encodes a protein involved in left-right body axis formation.
ENC1: Encodes a protein involved in the oxidative stress response. 
CHGB: Encodes a protein commonly found in neurons.
DDN: Encodes a protein involved in the transcription process.
SPP1: Involved in the regulation of interleukin-12, which is found in dendrites.
FOLH1: Encodes a protein involved in glutamate excitotoxicity, which negatively affects neurons.
UBASH3B: Encodes a protein that "attracts" receptors to  a cell.
IL1RL1: Encodes a protein that is involved in the interleukin receptor pathway.
CARNS1: Encodes a protein involved in the creation of carnosine, which is found in the central nervous system.
