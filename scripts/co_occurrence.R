library(cooccur)
df <- read.csv("enriched_pfam_a_matrix.csv", row.names=1, header=TRUE)
df2 <- df[colSums(df)>0] # remove contigs containing no pfam
data <- cooccur(mat=df2, type="spp_site", thresh=FALSE, spp_names = TRUE)
table <- prob.table(data)
write.table(table, "enriched_pfam_prob_table.csv", sep=",", row.names = FALSE, quote=FALSE)

table$score <- -log10(table$p_gt)           # add a score column = -log10(p_gt)
table$score <- as.numeric(table$score)
table[sapply(table, is.infinite)] <- 100    # replace inf with 100
positive <- table[table$score >=1.3,]       # select significant positive related rows with p_gt < 0.05 (score >=1.3)
positive <- positive[,c("sp1_name","sp2_name","score")]   # edit data
write.table(positive, "significant_positive_correlation.csv", sep=",", row.names=FALSE, quote=FALSE)

table$score2 <- -log10(table$p_lt)          # add p_lt score column
table$score2 <- as.numeric(table$score2)
table[sapply(table, is.infinite)] <- 100    # replace inf with 100
negative <- table[table$score2 >=1.3,]      # select significant negative related rows p_lt < 0.05
negative <- negative[,c("sp1_name","sp2_name","score2")]
write.table(negative, "significant_negative_correlation.csv", sep=",", row.names=FALSE, quote=FALSE)

# visualization can be done by cytospace
