library(cooccur)
setwd("/Users/wyang/Desktop/phage_CT_project/dry_lab/052620_co_occurrence_network/")
df <- read.csv("top_20_enriched_pfam_PA.csv", row.names=1, header=TRUE)
df2 <- df[colSums(df)>0] # remove contigs containing no pfam
data <- cooccur(mat=df2, type="spp_site", thresh=FALSE, spp_names = TRUE)
table <- prob.table(data)
write.table(table, "top_20_pfam_prob_table.tsv", sep="\t", row.names = FALSE, quote=FALSE)

table$score <- -log10(table$p_gt)           # add a score column = -log10(p_gt)
table$score <- as.numeric(table$score)
table[sapply(table, is.infinite)] <- 100    # replace inf with 100
positive <- table[table$score >=1.3,]       # select significant positive related rows with p_gt < 0.05 (score >=1.3)
positive <- positive[,c("sp1_name","sp2_name","score")]   # edit data
write.table(positive, "significant_positive_correlation.tsv", sep="\t", row.names=FALSE, quote=FALSE)

table$score2 <- -log10(table$p_lt)          # add p_lt score column
table$score2 <- as.numeric(table$score2)
table[sapply(table, is.infinite)] <- 100    # replace inf with 100
negative <- table[table$score2 >=1.3,]      # select significant negative related rows p_lt < 0.05
negative <- negative[,c("sp1_name","sp2_name","score2")]
write.table(negative, "significant_negative_correlation.tsv", sep="\t", row.names=FALSE, quote=FALSE)


