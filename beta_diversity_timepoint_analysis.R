#Script takes a beta diversity file and a meta data file splitting samples into 2 groups, Pre and Post Industrial.
#Calculates beta diversity spreads within and between time points and statistical analysis of these separations
#PCOA is also performed

library(reshape2)
library(dplyr)
library(tidyr)
library(ggplot2)
library(ggrepel)
library(car)
library(rstatix)
library(vegan)
library(ape)

#set our working directory correctly
setwd("~/Dropbox/microArc/Gwyn/PrePostIndustrial/allVsAllBeta")

#Read in the mapping file
metadata_pre_post_orig <- read.csv("~/Dropbox/microArc/Gwyn/PrePostIndustrial/allVsAllBeta/metadata_pre_post.csv")



#Melt the beta diversity matrix
#read in the file, remove the top triangle and then melt to get 3 columns
beta <- as.matrix(read.csv("~/Dropbox/microArc/Gwyn/PrePostIndustrial/allVsAllBeta/beta_diversity_matrix_all.csv", row.names=1))
beta[upper.tri(beta)] = NA
beta_melted <- reshape2::melt(beta, na.rm=T)

#remove the outlier from the metadata dataframe and the beta diversity dataframe
metadata_pre_post <- metadata_pre_post_orig %>%
  filter (sample_name != "ERR1681523")

beta_melted <- beta_melted %>%
  filter (Var1 != "ERR1681523" & Var2 != "ERR1681523")



#Add a column that is 0 for  pre or 1 for post based on Var1 and Var2
merged_beta <- beta_melted %>%
  left_join(metadata_pre_post, by = c("Var1" = "sample_name")) %>%
  left_join(metadata_pre_post, by = c("Var2" = "sample_name")) %>%
  mutate(grouping.x = ifelse(grouping.x == "Pre-Industrial", 0, ifelse(grouping.x == "Post-Industrial", 1, grouping.x))) %>%
  mutate(grouping.y = ifelse(grouping.y == "Pre-Industrial", 0, ifelse(grouping.y == "Post-Industrial", 1, grouping.y)))

#Sum the columns: 0 is pre vs pre, 1 is between, 2 is post vs post
merged_beta <- merged_beta %>%
  mutate(combined = as.numeric(grouping.x) + as.numeric(grouping.y))

#Change the 0, 1,2 in the combined column into Pre-Industrial, Between Timepoints and Post-Industrial
merged_beta <- merged_beta %>%
  mutate(combined = ifelse(combined == 0, "Pre-Industrial", ifelse(combined == 2, "Post-Industrial", ifelse(combined == 1, "Between Timepoints", combined))))


#Box and whisker based on columns (need to clean up axis labels and colours)
ggplot(merged_beta, aes(x = as.character(combined), y = value)) +
  geom_boxplot(aes(x = factor(combined, level = c("Pre-Industrial", "Post-Industrial","Between Timepoints"))))

#Shapiro-Wilk test for normality on each of the 3 groups
shapiro_results <- merged_beta %>%
  group_by(combined) %>%
  summarise(
    shapiro_p_value = shapiro.test(value)$p.value,
    shapiro_statistic = shapiro.test(value)$statistic
  )
print(shapiro_results)

#Levene's test for homogeny of variance
levene = leveneTest(value ~ combined, data = merged_beta)
print(levene)

#Data is not normal nor has homogeneity of variance so Kruskal-Wallis test is used

KW_result <- kruskal.test(value ~ combined, data = merged_beta)
print(KW_result)

#P-value is significant so do pairwise Wilcoxon tests
posthoc_results <- merged_beta %>%
  pairwise_wilcox_test(value ~ combined, p.adjust.method = "bonferroni")
print(posthoc_results)

#All significantly different

#PCOA of data
#Read back in the beta diversity data and replace the NA with zero on the diagonals
beta <- as.data.frame(read.csv("~/Dropbox/microArc/Gwyn/PrePostIndustrial/allVsAllBeta/beta_diversity_matrix_all.csv", row.names=1))
diag(beta) <- 0

bray_curtis_pcoa <- pcoa(beta)

# Extract coordinates
beta_pcoa_coords <- as.data.frame(bray_curtis_pcoa$vectors[, 1:2])
colnames(beta_pcoa_coords) <- c("PCoA1", "PCoA2")
beta_pcoa_coords$sample_name <- rownames(beta_pcoa_coords)

#Merge the PCOA with the timepoint mapping file
beta_pcoa_coords <- left_join(beta_pcoa_coords, metadata_pre_post_orig, by = "sample_name")

ggplot(beta_pcoa_coords, aes(x = PCoA1, y = PCoA2)) +
  geom_point(aes(color = grouping), size = 3) +
  geom_text_repel(aes(label = sample_name), max.overlaps = 1) +
  labs(title = "PCoA Plot", x = "PCoA1", y = "PCoA2") +
  theme_minimal()


#remove outlier
beta <- subset(beta, select = -c(ERR1681523))  %>%
  filter(! rownames(.) == "ERR1681523")

#Get the order of Pre/Post industrial labels for the rows
sample_order<- as.data.frame(rownames(beta))
group_order <- sample_order %>%
  left_join(metadata_pre_post, by = c("rownames(beta)" = "sample_name")) 
group_counts <- factor(group_order$grouping)

#Perform Beta Dispersion, Adonis and ANOSIM
matrix.bray_curtis <-
  vegdist(beta, method = "bray", na.rm = TRUE)
bray_curtis.disp <-
  betadisper(matrix.bray_curtis, group_counts )

bray_curtis.disp.anova <- anova(bray_curtis.disp) 
bray_curtis.disp.anova

bray_curtis.disp.TukeyHSD <- TukeyHSD(bray_curtis.disp)
bray_curtis.disp.TukeyHSD


bray_curtis.adonis <- adonis2(matrix.bray_curtis ~ group_counts)
bray_curtis.adonis

bray_curtis.anosim <- anosim(matrix.bray_curtis, group_counts)
bray_curtis.anosim

#Dispersion is not significant so the avriances are the same.
#Adnois (PERMANOVA) showed significant difference so mean beta diversity values are different between groups
#ANOSIM also showed significant difference so there is more similarity within groups than between

