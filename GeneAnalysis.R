#install packages (only need to do this once)
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("illuminaHumanv4.db")

library(illuminaHumanv4.db)

ec_id<-c("ILMN_1702255", "ILMN_1788874", "ILMN_1803429", "ILMN_2060413", "ILMN_2307903", "ILMN_1681805", "ILMN_1779147", "ILMN_1765966", "ILMN_1673450", "ILMN_1651354", "ILMN_1735115", "ILMN_1653856", "ILMN_2313672", "ILMN_3237981")
fc_id<-c("ILMN_1671777", "ILMN_1681805", "ILMN_1651496", "ILMN_2120624", "ILMN_1673450", "ILMN_1788874", "ILMN_1765966", "ILMN_1802633", "ILMN_1689088", "ILMN_2071809", "ILMN_2405756", "ILMN_1704376", "ILMN_3237981")
c_id<-c("ILMN_1712066", "ILMN_1910660", "ILMN_1673450", "ILMN_1803036", "ILMN_1670652", "ILMN_1712888", "ILMN_2100437", "ILMN_1742881", "ILMN_1765966", "ILMN_1813206", "ILMN_1764573", "ILMN_1685690")
tc_id<-c("ILMN_1676283", "ILMN_1732921", "ILMN_2218758", "ILMN_1670652", "ILMN_1839647", "ILMN_2207988", "ILMN_2255133", "ILMN_3245103", "ILMN_1736178", "ILMN_1653856", "ILMN_3248171", "ILMN_1679984", "ILMN_1691097", "ILMN_2172318", "ILMN_1788874")

#create a dataframe with illumina ids as gene names
data.frame(Gene=unlist(mget(x = ec_id,envir = illuminaHumanv4SYMBOL)))
data.frame(Gene=unlist(mget(x = fc_id,envir = illuminaHumanv4SYMBOL)))
data.frame(Gene=unlist(mget(x = c_id,envir = illuminaHumanv4SYMBOL)))
data.frame(Gene=unlist(mget(x = tc_id,envir = illuminaHumanv4SYMBOL)))
