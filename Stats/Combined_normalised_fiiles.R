getwd()
dir_g <- "./Glucose/HTSeq_results/"
# List normalized files
g_files <- list.files(path=dir_g,pattern = "_normalized.csv")
test<-read.csv("Glucose/HTSeq_results/SRR1166442_normalized.csv", header=FALSE)
test<-test[-1,]
sum <- grepl("^__",test[,1])
test_c <-test[!sum, ]
gene<-test_c[,1]
test<-test[,2]

g_files <- list.files(path="./Glucose/HTSeq_results/", pattern="*_normalized.csv", full.names = TRUE, recursive=FALSE)
the_files_g <- data.frame(gene=gene, stringsAsFactors = TRUE)

for (i in g_files){
  the_file_g <- as.data.frame(read.csv(i, header = F, stringsAsFactors = F))
  the_file_g <- the_file_g[-1,] #remove the column name
  sum_rows<-grepl("^__",the_file_g[,1])#Remove the HTSeq summary rows
  the_file_g <- the_file_g[!sum_rows, ]
  rpkm <- the_file_g[,3] #extra the rpkm value 
  sample_name <- sub("_normalized.csv", "", basename(i))
  the_files_g[[sample_name]] <- rpkm
}
dir_c <- "./Cellobiose/HTSeq_results/"
c_files <-list.files(path = dir_c, pattern = "*_normalized.csv", full.names = TRUE, recursive = FALSE)
the_files_c <- data.frame(gene=gene, stringsAsFactors = T)
for (m in c_files){
  the_file_c <- as.data.frame(read.csv(m, header = F, stringsAsFactors = F))
  the_file_c <- the_file_c[-1,]
  sum_rows_c <- grepl("^__", the_file_c[,1])
  the_file_c <- the_file_c[!sum_rows_c, ]
  rpkm_c <- the_file_c[,3]
  sample_name_c <- sub("_normalized.csv", "",basename(m))
  the_files_c[[sample_name_c]] <- rpkm_c
}

com_file_c <- the_files_c[,-1]
all_file <- cbind(the_files_g,com_file_c)

colnames(all_file) <- c("gene", "JCYL001B", "JCYL002B", "JCYL003B", "JCYL001D", "JCYL002D", "JCYL003D")

write.csv(all_file, "Combined_RPKM.csv", row.names = FALSE)

#Calculate means
data <- read.csv("Combined_RPKM.csv", header = T)
data$Mean_B <- rowMeans(data[,c("JCYL001B", "JCYL002B", "JCYL003B")], na.rm = T)
data$Mean_D <- rowMeans(data[, c("JCYL001D", "JCYL002D", "JCYL003D")], na.rm = T)
colnames(data)
data <- data[, c("gene", "JCYL001B", "JCYL002B", "JCYL003B", "Mean_B", 
                 "JCYL001D", "JCYL002D", "JCYL003D", "Mean_D")]
print(data)
write.csv(data, "combined_with_means.csv", row.names = FALSE)
