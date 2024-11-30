# Author: Dr. Elina Takola
# Date: 30 November 2024
# Publilcation reference: Manuscript in prep.

#load("SOMA-Aug24.RData")

sessionInfo()
options(scipen = 999, digits = 2) # disable scientific notation and ask for 2 decimals
library(metafor)
library(plyr) # always load before dplyr
library(dplyr)
library(tidyverse)
library(tidyr)
library(ggplot2)
library(ggthemes)
library(paletteer)
library(ggpubr)

###################################
# Import data
###################################

# Correlation matrix of overlap between studies
A <- read.csv(file.choose(), header = TRUE, row.names = 1, sep = ",") # Load the file "overlap_matrix_primary_percentage.csv". It contains the count of common primary studies for each pair of meta-analyses. Diagonal is based on the sample size of two meta-analyses that have almost 100% overlap.
A1 <- data.matrix(A, rownames.force = TRUE) # matrix shoud be square, symmetric, numeric

# Import dataset with pairs
df <- read.csv(file.choose(), header = TRUE, sep = ",") # soma_data_final.csv
df3 = df

# We also need the dataset with the effect sizes for Biodiversity
biodiv <- read.csv(file.choose(), header = TRUE, sep = ",")  # soma_biodiv_final.csv
biodiv3 = biodiv

# We will also need the dataset with yield effect sizes
yield <- read.csv(file.choose(), header = TRUE, sep = ",") # soma_yield_final.csv
yield3 = yield

###################################
# 1)	What is the impact of sustainable farming practices on biodiversity and yield?
###################################
df_lrr <- subset(df3, Biodiv_ESMetric %in% c("LRR","Percentage change", "RR"))

line <- lm(df_lrr$Prod_ES_homog ~ df_lrr$Biodiv_ES_homog)
summary(line)
# Plot the linear regression
ggplot(data = df_lrr, aes(x = Biodiv_ES_homog, y = Prod_ES_homog, ymax = 1.5)) +
  labs(y = "LRR(Yield)", x = "LRR(Biodiversity)") +
  geom_smooth(method = "lm", 
              formula = y ~ x, fill = "#8babf1", color = "#f57600") +
  geom_point() +
  geom_hline(yintercept=0, linetype="solid", color = "gray", linewidth=0.7) + # Draw zero line
  geom_vline(xintercept=0, linetype="solid", color = "grey", linewidth=0.7) + # Draw zero line
  theme_few()

# Plot the linear regression with Management colour codes
ggplot(data = df_lrr, aes(x = Biodiv_ES_homog, y = Prod_ES_homog, ymax = 1.5, color = Management_grouped)) +
  labs(y = "LRR(Yield)", x = "LRR(Biodiversity)") +
  geom_smooth(method = "lm", 
              formula = y ~ x, 
              se = FALSE) +
  scale_color_manual(values = c("diversification" = "#c44601",                     
                                "lower intensity" = "#f57600", 
                                "reduced resource addition" = "#8babf1",           
                                "grazing pause" = "#0073e6", 
                                "organic" = "#054fb9")) + 
  geom_point() +
  geom_hline(yintercept=0, linetype="solid", color = "gray", linewidth=0.7) + # Draw zero line
  geom_vline(xintercept=0, linetype="solid", color = "grey", linewidth=0.7) + # Draw zero line
  theme_few()

# Overall meta-analyses for biodiversity and yield
yield <- subset(yield, Prod_ESMetric %in% c("LRR","Percentage change", "RR"))
biodiv <- subset(biodiv, Biodiv_ESMetric %in% c("LRR","Percentage change", "RR"))
A1_yield_lrr <- A1[unique(yield$StudyID), unique(yield$StudyID)]
A1_biodiv_lrr <- A1[unique(biodiv$StudyID), unique(biodiv$StudyID)]

meta_yield <- rma.mv(yi = Prod_ES_homog, V = Prod_SE3, 
                       random = list(~1|StudyID, ~1|ES_ID),
                       method = "REML", R = list(StudyID = A1_yield_lrr),
                       data = yield)
summary(meta_yield)

meta_bio <- rma.mv(yi = Biodiv_ES_homog, V = Biodiv_SE3, 
                                random = list(~1|StudyID, ~1|ES_ID),
                                method = "REML", R = list(StudyID = A1_biodiv_lrr),
                                data = biodiv)
summary(meta_bio)

##############################
# 2)	Which management practices lead to win-win outcomes for biodiversity and yield?
##############################

# Create a subset for each taxon
invertebrates <- biodiv[biodiv$Taxonomic_group == "invertebrates",]
plants <- biodiv[biodiv$Taxonomic_group == "plants",]
fungi <- biodiv[biodiv$Taxonomic_group == "fungi",]
#microbes <- biodiv[biodiv$Taxonomic_group == "microbes",]
# Can't do a meta-regression for microbes because Management is only organic (1 level)
vertebrates <- biodiv[biodiv$Taxonomic_group == "vertebrates",]

# Subset the matrix of primary studies for each subset
A1_invertebrates <- A1[unique(invertebrates$StudyID), unique(invertebrates$StudyID)]
A1_plants <- A1[unique(plants$StudyID), unique(plants$StudyID)]
A1_fungi <- A1[unique(fungi$StudyID), unique(fungi$StudyID)]
A1_vertebrates <- A1[unique(vertebrates$StudyID), unique(vertebrates$StudyID)]
# And then we will fit random effects-only meta-analytic models and meta-regressions separately for each taxon using management as a moderator

# Overall meta-regression for all biodiversity taxa and management practices
metareg_bio <- rma.mv(yi = Biodiv_ES_homog, V = Biodiv_SE3, mods = ~ 0 + Management_grouped,
                                random = list(~1|StudyID, ~1|ES_ID),
                                method = "REML", R = list(StudyID = A1_biodiv_lrr),
                                data = biodiv)
summary(metareg_bio)

#Invertebrates
meta_invertebrates <- rma.mv(yi = Biodiv_ES_homog, V = Biodiv_SE3,
                   random = list(~1|StudyID, ~1|ES_ID),
                   method = "REML", R = list(StudyID = A1_invertebrates),
                   data = invertebrates)
summary(meta_invertebrates)
predict(meta_invertebrates)

metareg_invertebrates <- rma.mv(yi = Biodiv_ES_homog, V = Biodiv_SE3, mods = ~ 0 + Management_grouped,
                     random = list(~1|StudyID, ~1|ES_ID),
                     method = "REML", R = list(StudyID = A1_invertebrates),
                     data = invertebrates)
summary(metareg_invertebrates)

#Plants
meta_plants <- rma.mv(yi = Biodiv_ES_homog, V = Biodiv_SE3,
                   random = list(~1|StudyID, ~1|ES_ID),
                   method = "REML", R = list(StudyID = A1_plants),
                   data = plants)
summary(meta_plants)
predict(meta_plants)

metareg_plants <- rma.mv(yi = Biodiv_ES_homog, V = Biodiv_SE3, mods = ~ 0 + Management_grouped,
                      random = list(~1|StudyID, ~1|ES_ID),
                      method = "REML", R = list(StudyID = A1_plants),
                      data = plants)
summary(metareg_plants)

#Fungi
meta_fungi <- rma.mv(yi = Biodiv_ES_homog, V = Biodiv_SE3,
                      random = list(~1|StudyID, ~1|ES_ID),
                      method = "REML", R = list(StudyID = A1_fungi),
                      data = fungi)
summary(meta_fungi)
predict(meta_fungi)

metareg_fungi <- rma.mv(yi = Biodiv_ES_homog, V = Biodiv_SE3, mods = ~ 0 + Management_grouped,
                         random = list(~1|StudyID, ~1|ES_ID),
                         method = "REML", R = list(StudyID = A1_fungi),
                         data = fungi)
summary(metareg_fungi)

#Vertebrates
meta_vertebrates <- rma.mv(yi = Biodiv_ES_homog, V = Biodiv_SE3,
                             random = list(~1|StudyID, ~1|ES_ID),
                             method = "REML", R = list(StudyID = A1_vertebrates),
                             data = vertebrates)
summary(meta_vertebrates)
predict(meta_vertebrates)

metareg_vertebrates <- rma.mv(yi = Biodiv_ES_homog, V = Biodiv_SE3, mods = ~ 0 + Management_grouped,
                                random = list(~1|StudyID, ~1|ES_ID),
                                method = "REML", R = list(StudyID = A1_vertebrates),
                                data = vertebrates)
summary(metareg_vertebrates)

# Gather all moderator estimates in a data frame
meta_bio <- as.data.frame(cbind(rownames = rownames(metareg_bio), metareg_bio$b, metareg_bio$ci.lb, metareg_bio$ci.ub))
colnames(meta_bio) <- c("Estimate", "CI.LB", "CI.UB")
meta_bio <- tibble::rownames_to_column(meta_bio, "Moderator")
meta_bio$Moderator <- c("Diversification", "Grazing pause", "Lower intensity", "Organic", "Reduced resource addition")
meta_bio$Taxon <- "All taxa"

meta_bio_inv <- as.data.frame(cbind(rownames = rownames(metareg_invertebrates), metareg_invertebrates$b, metareg_invertebrates$ci.lb, metareg_invertebrates$ci.ub))
colnames(meta_bio_inv) <- c("Estimate", "CI.LB", "CI.UB")
meta_bio_inv <- tibble::rownames_to_column(meta_bio_inv, "Moderator")
meta_bio_inv$Moderator <- c("Diversification", "Lower intensity", "Organic")
meta_bio_inv$Taxon <- "Invertebrates"

meta_bio_pl <- as.data.frame(cbind(rownames = rownames(metareg_plants), metareg_plants$b, metareg_plants$ci.lb, metareg_plants$ci.ub))
colnames(meta_bio_pl) <- c("Estimate", "CI.LB", "CI.UB")
meta_bio_pl <- tibble::rownames_to_column(meta_bio_pl, "Moderator")
meta_bio_pl$Moderator <- c("Diversification", "Grazing pause", "Lower intensity", "Organic", "Reduced resource addition")
meta_bio_pl$Taxon <- "Plants"

meta_bio_fungi <- as.data.frame(cbind(rownames = rownames(metareg_fungi), metareg_fungi$b, metareg_fungi$ci.lb, metareg_fungi$ci.ub))
colnames(meta_bio_fungi) <- c("Estimate", "CI.LB", "CI.UB")
meta_bio_fungi <- tibble::rownames_to_column(meta_bio_fungi, "Moderator")
meta_bio_fungi$Moderator <- c("Diversification", "Organic", "Reduced resource addition")
meta_bio_fungi$Taxon <- "Fungi"

meta_bio_vert <- as.data.frame(cbind(rownames = rownames(metareg_vertebrates), metareg_vertebrates$b, metareg_vertebrates$ci.lb, metareg_vertebrates$ci.ub))
colnames(meta_bio_vert) <- c("Estimate", "CI.LB", "CI.UB")
meta_bio_vert <- tibble::rownames_to_column(meta_bio_vert, "Moderator")
meta_bio_vert$Moderator <- c("Lower intensity", "Organic")
meta_bio_vert$Taxon <- "Vertebrates"

mod_bio <- rbind(meta_bio, meta_bio_inv, meta_bio_pl, meta_bio_fungi, meta_bio_vert)

# Plot moderators for Biodiversity
library(rcartocolor)
mod_bio_plot <- ggplot(mod_bio, aes(x = Estimate, y = Moderator)) +
  geom_errorbar(aes(xmin = CI.LB, xmax = CI.UB, width = 0.1)) +
  geom_point(aes(colour = Moderator), size = 5) +
  facet_grid(. ~ Taxon, switch = "y", scales = "free") +
  geom_vline(xintercept = 0, linetype = "solid", color = "black", linewidth = 0.7) + # Draw zero line vertical
  theme_minimal() +  
  scale_color_manual(values = c("Diversification" = "#c44601",                     
                                "Lower intensity" = "#f57600", 
                                "Reduced resource addition" = "#8babf1",           
                                "Grazing pause" = "#0073e6", 
                                "Organic" = "#054fb9")) 
mod_bio_plot


######################################
#  3) Which taxa and crop systems are benefited by specific agricultural management practices?
######################################
crop <- yield[yield$CropType_grouped == "cropland",]
agrof <- yield[yield$CropType_grouped == "(agro)forestry",]
grass <- yield[yield$CropType_grouped == "grassland",]

A1_crop <- A1[unique(crop$StudyID), unique(crop$StudyID)]
A1_agrof <- A1[unique(agrof$StudyID), unique(agrof$StudyID)]
A1_grass <- A1[unique(grass$StudyID), unique(grass$StudyID)]

metareg_yield <- rma.mv(yi = Prod_ES_homog, V = Prod_SE3, mods = ~ 0 + Management_grouped,
                       random = list(~1|StudyID, ~1|ES_ID),
                       method = "REML", R = list(StudyID = A1_yield_lrr),
                       data = yield)
summary(metareg_yield)

# Cropland
metareg_crop <- rma.mv(yi = Prod_ES_homog, V = Prod_SE3, mods = ~ 0 + Management_grouped,
                              random = list(~1|StudyID, ~1|ES_ID),
                              method = "REML", R = list(StudyID = A1_crop),
                              data = crop)
summary(metareg_crop)

meta_crop <- rma.mv(yi = Prod_ES_homog, V = Prod_SE3, 
                       random = list(~1|StudyID, ~1|ES_ID),
                       method = "REML", R = list(StudyID = A1_crop),
                       data = crop)
summary(meta_crop)
predict(meta_crop)


# Agroforestry
meta_agrof <- rma.mv(yi = Prod_ES_homog, V = Prod_SE3,
                    random = list(~1|StudyID, ~1|ES_ID),
                    method = "REML", R = list(StudyID = A1_agrof),
                    data = agrof)
summary(meta_agrof)
predict(meta_agrof)

metareg_agrof <- rma.mv(yi = Prod_ES_homog, V = Prod_SE3, mods = ~ 0 + Management_grouped,
                       random = list(~1|StudyID, ~1|ES_ID),
                       method = "REML", R = list(StudyID = A1_agrof),
                       data = agrof)
summary(metareg_agrof)

# Grassland
meta_grass <- rma.mv(yi = Prod_ES_homog, V = Prod_SE3,
                     random = list(~1|StudyID, ~1|ES_ID),
                     method = "REML", R = list(StudyID = A1_grass),
                     data = grass)
summary(meta_grass)
predict(meta_grass)

metareg_grass <- rma.mv(yi = Prod_ES_homog, V = Prod_SE3, mods = ~ 0 + Management_grouped,
                        random = list(~1|StudyID, ~1|ES_ID),
                        method = "REML", R = list(StudyID = A1_grass),
                        data = grass)
summary(metareg_grass)

# Gather all moderator estimates in a data frame
meta_yield <- as.data.frame(cbind(rownames = rownames(metareg_yield), metareg_yield$b, metareg_yield$ci.lb, metareg_yield$ci.ub))
colnames(meta_yield) <- c("Estimate", "CI.LB", "CI.UB")
meta_yield <- tibble::rownames_to_column(meta_yield, "Moderator")
meta_yield$Moderator <- c("Diversification", "Grazing pause", "Lower intensity", "Organic", "Reduced resource addition")
meta_yield$Crop_system <- "All crop systems"


meta_yield_crop <- as.data.frame(cbind(rownames = rownames(metareg_crop), metareg_crop$b, metareg_crop$ci.lb, metareg_crop$ci.ub))
colnames(meta_yield_crop) <- c("Estimate", "CI.LB", "CI.UB")
meta_yield_crop <- tibble::rownames_to_column(meta_yield_crop, "Moderator")
meta_yield_crop$Moderator <- c("Diversification", "Lower intensity", "Organic", "Reduced resource addition")
meta_yield_crop$Crop_system <- "Cropland"

meta_yield_agrof <- as.data.frame(cbind(rownames = rownames(metareg_agrof), metareg_agrof$b, metareg_agrof$ci.lb, metareg_agrof$ci.ub))
colnames(meta_yield_agrof) <- c("Estimate", "CI.LB", "CI.UB")
meta_yield_agrof <- tibble::rownames_to_column(meta_yield_agrof, "Moderator")
meta_yield_agrof$Moderator <- c("Diversification", "Organic")
meta_yield_agrof$Crop_system <- "Agroforestry"

meta_yield_grass <- as.data.frame(cbind(rownames = rownames(metareg_grass), metareg_grass$b, metareg_grass$ci.lb, metareg_grass$ci.ub))
colnames(meta_yield_grass) <- c("Estimate", "CI.LB", "CI.UB")
meta_yield_grass <- tibble::rownames_to_column(meta_yield_grass, "Moderator")
meta_yield_grass$Moderator <- c("Grazing pause", "Lower intensity", "Organic", "Reduced resource addition")
meta_yield_grass$Crop_system <- "Grassland"

mod_yield <- rbind(meta_yield, meta_yield_crop, meta_yield_agrof, meta_yield_grass)
mod_yield$Crop_system <- factor(mod_yield$Crop_system, levels = unique(mod_yield$Crop_system)) # Convert this column to a factor because I want them to appear in the facet grid as they are in the dataframe, not alphabetically

library(rcartocolor)
mod_yie_plot <- ggplot(mod_yield,aes(x = Estimate, y = Moderator)) +
  geom_errorbar(aes(xmin = CI.LB, xmax = CI.UB, width = 0.1)) +
  geom_point(aes(colour = Moderator), size = 5) +
  facet_grid(. ~ Crop_system, switch = "y", scales = "free") +
  geom_vline(xintercept = 0, linetype = "solid", color = "black", linewidth = 0.7) + # Draw zero line vertical
  theme_minimal() +  
  scale_color_manual(values = c("Diversification" = "#c44601",                     
                                "Lower intensity" = "#f57600", 
                                "Reduced resource addition" = "#8babf1",           
                                "Grazing pause" = "#0073e6", 
                                "Organic" = "#054fb9")) 

mod_bio_plot2 <- mod_bio_plot + theme(legend.position = "none", 
                                      axis.text.x = element_text(size = 14),
                                      axis.text.y = element_text(size = 18),
                                      strip.text = element_text(size = 16),
                                      panel.spacing = unit(0.8, "cm"))

mod_yie_plot2 <- mod_yie_plot + theme(legend.position = "none", 
                                      axis.ticks.x=element_blank(),
                                      axis.text.y=element_blank(),
                                      axis.text.x = element_text(size = 14),
                                      axis.ticks.y=element_blank(),
                                      axis.title.y=element_blank(),
                                      strip.text = element_text(size = 16),
                                      panel.spacing = unit(0.8, "cm"))
mods <- ggarrange(mod_bio_plot2, mod_yie_plot2, ncol = 2, nrow = 1) #labels = c("A)", "B)"),
mods

#######################################
# Knowledge gaps
#######################################
library(pheatmap)
heat <- table(biodiv$Taxonomic_group, biodiv$Management_grouped)
pheatmap(heat, col = colorRampPalette(c("#FCEEEA", "#c44601"))(100), 
         main = "Evidence map
         Sustainable farming impact on biodiversity", cellwidth = 100, cellheight = 40, fontsize = 12, 
         number_format = "%.0f", display_numbers = TRUE, cluster_cols = FALSE, cluster_rows = FALSE,
         angle_col = 45)

heat <- table(yield$CropType_grouped, yield$Management_grouped)
pheatmap(heat, col = colorRampPalette(c("aliceblue", "#0073e6"))(100), 
         main = "Evidence map
         Sustainable farming impact on yield", cellwidth = 100, cellheight = 40, fontsize = 12, 
         number_format = "%.0f", display_numbers = TRUE, cluster_cols = FALSE, cluster_rows = FALSE,
         angle_col = 45)

heat <- table(biodiv$Taxonomic_group, biodiv$CropType)
pheatmap(heat, col = colorRampPalette(c("#EFEBF5", "#6A4C93"))(100), 
         main = "Evidence map
         Sustainable farming impact on specific taxa and crop types", cellwidth = 100, cellheight = 40, fontsize = 12, 
         number_format = "%.0f", display_numbers = TRUE, cluster_cols = FALSE, cluster_rows = FALSE,
         angle_col = 45)

#####################################
# Above- vs Below-ground
#####################################
abg <- ggplot(data = df_lrr,aes(x = Prod_ES_homog,y = Biodiv_ES_homog, color = AboveBelowGround, label=StudyID)) + 
  geom_point(aes(size = Prod_sample_size), alpha = 0.95) + 
  geom_errorbar(aes(ymin = Biodiv_ES_homog-Biodiv_SE3,ymax = Biodiv_ES_homog+Biodiv_SE3), width = 0.05, linewidth = 0.7) + 
  geom_errorbarh(aes(xmin = Prod_ES_homog-Prod_SE3, xmax = Prod_ES_homog+Prod_SE3), height = 0.05, linewidth = 0.7) +
  scale_color_manual(values = c("Below" = "#054fb9", "Both" = "#f57600",
                                "Above" = "#8babf1")) + 
  geom_hline(yintercept=0, linetype="solid", color = "black", linewidth=0.7) + # Draw zero line
  geom_vline(xintercept=0, linetype="solid", color = "black", linewidth=0.7) + # Draw zero line
  theme_few()
abg



sub <- subset(df_lrr, AboveBelowGround %in% "Above")
above <- ggplot(data = sub,aes(x = Prod_ES_homog,y = Biodiv_ES_homog, color = AboveBelowGround, label=StudyID)) + 
  ggtitle("A. Above-ground") +
  geom_point(aes(size = Prod_sample_size), alpha = 0.95) + 
  geom_errorbar(aes(ymin = Biodiv_ES_homog-Biodiv_SE3,ymax = Biodiv_ES_homog+Biodiv_SE3), width = 0.05, linewidth = 0.7) + 
  geom_errorbarh(aes(xmin = Prod_ES_homog-Prod_SE3, xmax = Prod_ES_homog+Prod_SE3), height = 0.05, linewidth = 0.7) +
  scale_color_manual(values = c("Below" = "#054fb9", "Both" = "#054fb9",
                                "Above" = "#054fb9")) + 
  geom_hline(yintercept=0, linetype="solid", color = "black", linewidth=0.7) + # Draw zero line
  geom_vline(xintercept=0, linetype="solid", color = "black", linewidth=0.7) + # Draw zero line
  theme_few() + theme(legend.position = "none")
above

sub <- subset(df_lrr, AboveBelowGround %in% "Below")
below <- ggplot(data = sub,aes(x = Prod_ES_homog,y = Biodiv_ES_homog, color = AboveBelowGround, label=StudyID)) + 
  ggtitle("B. Below-ground") +
  geom_point(aes(size = Prod_sample_size), alpha = 0.95) + 
  geom_errorbar(aes(ymin = Biodiv_ES_homog-Biodiv_SE3,ymax = Biodiv_ES_homog+Biodiv_SE3), width = 0.05, linewidth = 0.7) + 
  geom_errorbarh(aes(xmin = Prod_ES_homog-Prod_SE3, xmax = Prod_ES_homog+Prod_SE3), height = 0.05, linewidth = 0.7) +
  scale_color_manual(values = c("Below" = "#054fb9", "Both" = "#054fb9",
                                "Above" = "#054fb9")) + 
  geom_hline(yintercept=0, linetype="solid", color = "black", linewidth=0.7) + # Draw zero line
  geom_vline(xintercept=0, linetype="solid", color = "black", linewidth=0.7) + # Draw zero line
  theme_few() + theme(legend.position = "none")
below

ggarrange(above, below)
