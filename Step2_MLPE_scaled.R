# load libraries
library(ecodist)
library(tidyverse)
library(lme4)
library(Matrix)
library(MuMIn)
library(usdm)

setwd("~/Box Sync/FLL/Analyses/Models")

###IMPORT FIXED EFFECTS### this should be repeated for each HS transformation calculated in Step 1a

#Least Cost Path Distance
LCP_csv <- read.csv("~/Box Sync/FLL/Analyses/Models/Fixed/LCPnlin8.dist.csv", header = TRUE, row.names = 1)
lcprnames<- rownames(LCP_csv) #overwrite row names to correct format (also use later if needed)
colnames(LCP_csv) <- lcprnames
LCP_mtrx <- as.matrix(LCP_csv) #change to matrix
LCP_ord <- arrange(LCP_mtrx, order.rows = TRUE, order.cols = TRUE, comparator = NULL)
df_LCP <- matrixConvert(LCP_ord, colname = c("site_1", "site_2", "LCP"))
df_data <- df_LCP ##dataframe to which all other data will be added

#Circuitscape resistance distance
CS_csv <- read.csv("~/Box Sync/FLL/Analyses/Models/Fixed/rdist_nlin8CS.csv", header = TRUE, row.names = 1)
csrnames<- rownames(CS_csv) #overwrite row names to correct format (also use later if needed)
colnames(CS_csv) <- csrnames
CS_mtrx <- as.matrix(CS_csv) #change to matrix
CS_ord <- arrange(CS_mtrx, order.rows = TRUE, order.cols = TRUE, comparator = NULL)
df_CS <- matrixConvert(CS_ord, colname = c("site_1", "site_2", "CS"))

#Pairwise geographic euclidean distance (easier calculation below):
#EUC_csv <- read.csv("~/Box Sync/FLL/Analyses/Models/Fixed/euc_dist142.csv", header = TRUE, row.names = 1)
#eucrnames<- rownames(EUC_csv) #overwrite row names to correct format (also use later if needed)
#colnames(EUC_csv) <- eucrnames
#EUC_mtrx <- as.matrix(EUC_csv) #change to matrix
#EUC_ord <- arrange(EUC_mtrx, order.rows = TRUE, order.cols = TRUE, comparator = NULL)
#df_EUC <- matrixConvert(EUC_ord, colname = c("site_1", "site_2", "EUC"))

#View(df_EUC)

sample_points <- 
  paste0("~/Box Sync/FLL/Analyses/Circuitscape/142nDNA.shp") %>%
  raster::shapefile() %>%
  tibble::as_tibble()

#View(sample_points)
EUC_matrix <- as.matrix(dist(cbind(sample_points$Lat,sample_points$Long)))
rownames(EUC_matrix) <- sample_points$Ele_ID
colnames(EUC_matrix) <- sample_points$Ele_ID
EUC_ord2 <- arrange(EUC_matrix, order.rows = TRUE, order.cols = TRUE, comparator = NULL)
df_EUC <- matrixConvert(EUC_ord2, colname = c("site_1", "site_2", "EUC"))



##IMPORT GENETIC DISTANCE METRICS FROM STEP 1B###

#GD as 1 minus the Portion of Shared Alleles
PS_csv <- read.csv("~/Box Sync/FLL/Analyses/Models/Response/GD_1minPS.csv", header = TRUE, row.names = 1)
PSrnames<- rownames(PS_csv) #overwrite row names to correct format (also use later if needed)
colnames(PS_csv) <- PSrnames
PS_mtrx <- as.matrix(PS_csv) #change to matrix
PS_ord <- arrange(PS_mtrx, order.rows = TRUE, order.cols = TRUE, comparator = NULL)
df_PS <- matrixConvert(PS_ord, colname = c("site_1", "site_2", "PS"))

#GD as number of allelic differences between two individuals
aldiff_csv <- read.csv("~/Box Sync/FLL/Analyses/Models/Response/GD_aldiff.csv", header = TRUE, row.names = 1)
aldiffrnames<- rownames(aldiff_csv) #overwrite row names to correct format (also use later if needed)
colnames(aldiff_csv) <- aldiffrnames
aldiff_mtrx <- as.matrix(aldiff_csv) #change to matrix
aldiff_ord <- arrange(aldiff_mtrx, order.rows = TRUE, order.cols = TRUE, comparator = NULL)
df_aldiff <- matrixConvert(aldiff_ord, colname = c("site_1", "site_2", "aldiff"))

#GD as euclidean distance among vector of allele frequencies
gdeuc_csv <- read.csv("~/Box Sync/FLL/Analyses/Models/Response/GD_euc.csv", header = TRUE, row.names = 1)
gdeucrnames<- rownames(gdeuc_csv) #overwrite row names to correct format (also use later if needed)
colnames(gdeuc_csv) <- gdeucrnames
gdeuc_mtrx <- as.matrix(gdeuc_csv) #change to matrix
gdeuc_ord <- arrange(gdeuc_mtrx, order.rows = TRUE, order.cols = TRUE, comparator = NULL)
df_gdeuc <- matrixConvert(gdeuc_ord, colname = c("site_1", "site_2", "gdeuc"))

#GD as Reynolds's distance with the reynolds.dist function in poppr
gdtot_csv <- read.csv("~/Box Sync/FLL/Analyses/Models/Response/GD_reyn.csv", header = TRUE, row.names = 1)
gdtotrnames<- rownames(gdtot_csv) #overwrite row names to correct format (also use later if needed)
colnames(gdtot_csv) <- gdtotrnames
gdtot_mtrx <- as.matrix(gdtot_csv) #change to matrix
gdtot_ord <- arrange(gdtot_mtrx, order.rows = TRUE, order.cols = TRUE, comparator = NULL)
df_gdtot <- matrixConvert(gdtot_ord, colname = c("site_1", "site_2", "gdtot"))


###create master data table as df_data###
# add distance datasets
df_data <- dplyr::left_join(df_data, df_CS, 
                               by = c("site_1", "site_2"))
df_data <- dplyr::left_join(df_data, df_EUC, 
                            by = c("site_1", "site_2"))

# add genetic datasets
df_data <- dplyr::left_join(df_data, df_PS, 
                               by = c("site_1", "site_2"))
df_data <- dplyr::left_join(df_data, df_aldiff, 
                            by = c("site_1", "site_2"))
df_data <- dplyr::left_join(df_data, df_gdeuc, 
                            by = c("site_1", "site_2"))
df_data <- dplyr::left_join(df_data, df_gdtot, 
                            by = c("site_1", "site_2"))

df_data_scaled <- dplyr::mutate(df_data, EUC = as.numeric(scale(EUC)), LCP = as.numeric(scale(LCP)), CS = as.numeric(scale(CS)))
View(df_data)
View(df_data_scaled) #verify data table has all information and fixed effects have been scaled

# Create ZZ matrix
#ZZ_surface <- create_ZZ(data = df_surface, 
#                        col_names = c("site_1", "site_2"))


#fit MLPE without REML for AIC and BIC calculation
#with EUC dist & 4 gen datasets:
EUC_PS_noreml <- lmer(PS~EUC + (1|site_1), data = df_data_scaled, REML = FALSE )
summary(EUC_PS_noreml)
EUC_aldiff_noreml <- lmer(aldiff~EUC + (1|site_1), data = df_data_scaled, REML = FALSE )
summary(EUC_aldiff_noreml)
EUC_gdeuc_noreml <- lmer(gdeuc~EUC + (1|site_1), data = df_data_scaled, REML = FALSE )
summary(EUC_gdeuc_noreml)
EUC_gdtot_noreml <- lmer(gdtot~EUC + (1|site_1), data = df_data_scaled, REML = FALSE )
summary(EUC_gdtot_noreml)

#with LCP dist & 4 gen datasets:
LCP_PS_noreml <- lmer(PS~LCP +EUC + (1|site_1), data = df_data_scaled, REML = FALSE )
summary(LCP_PS_noreml)
LCP_aldiff_noreml <- lmer(aldiff~LCP +EUC + (1|site_1), data = df_data_scaled, REML = FALSE )
summary(LCP_aldiff_noreml)
LCP_gdeuc_noreml <- lmer(gdeuc~LCP +EUC + (1|site_1), data = df_data_scaled, REML = FALSE )
summary(LCP_gdeuc_noreml)
LCP_gdtot_noreml <- lmer(gdtot~LCP +EUC + (1|site_1), data = df_data_scaled, REML = FALSE )
summary(LCP_gdtot_noreml)

#with CS dist & 4 gen datasets:
#with LCP dist & 4 gen datasets:
CS_PS_noreml <- lmer(PS~ CS +EUC + (1|site_1), data = df_data_scaled, REML = FALSE )
summary(CS_PS_noreml)
CS_aldiff_noreml <- lmer(aldiff~CS +EUC + (1|site_1), data = df_data_scaled, REML = FALSE )
summary(CS_aldiff_noreml)
CS_gdeuc_noreml <- lmer(gdeuc~CS +EUC + (1|site_1), data = df_data_scaled, REML = FALSE )
summary(CS_gdeuc_noreml)
CS_gdtot_noreml <- lmer(gdtot~CS +EUC + (1|site_1), data = df_data_scaled, REML = FALSE )
summary(CS_gdtot_noreml)

#fit MLPE with REML for R2 calculations in MuMIn:
#with EUC dist & 4 gen datasets:
EUC_PS_reml <- lmer(PS~EUC + (1|site_1), data = df_data_scaled, REML = TRUE )
r2_EUC_PS <- r.squaredGLMM(EUC_PS_reml)
r2_EUC_PS
EUC_aldiff_reml <- lmer(aldiff~EUC + (1|site_1), data = df_data_scaled, REML = TRUE )
r2_EUC_aldiff <- r.squaredGLMM(EUC_aldiff_reml)
r2_EUC_aldiff
EUC_gdeuc_reml <- lmer(gdeuc~EUC + (1|site_1), data = df_data_scaled, REML = TRUE )
r2_EUC_gdeuc <- r.squaredGLMM(EUC_gdeuc_reml)
r2_EUC_gdeuc
EUC_gdtot_reml <- lmer(gdtot~EUC + (1|site_1), data = df_data_scaled, REML = TRUE )
r2_EUC_gdtot <- r.squaredGLMM(EUC_gdtot_reml)
r2_EUC_gdtot

#with LCP dist & 4 gen datasets:
LCP_PS_reml <- lmer(PS~LCP +EUC + (1|site_1), data = df_data_scaled, REML = TRUE )
r2_LCP_PS <- r.squaredGLMM(LCP_PS_reml)
r2_LCP_PS
LCP_aldiff_reml <- lmer(aldiff~LCP +EUC + (1|site_1), data = df_data_scaled, REML = TRUE )
r2_LCP_aldiff <- r.squaredGLMM(LCP_aldiff_reml)
r2_LCP_aldiff
LCP_gdeuc_reml <- lmer(gdeuc~LCP +EUC + (1|site_1), data = df_data_scaled, REML = TRUE )
r2_LCP_gdeuc <- r.squaredGLMM(LCP_gdeuc_reml)
r2_LCP_gdeuc
LCP_gdtot_reml <- lmer(gdtot~LCP +EUC + (1|site_1), data = df_data_scaled, REML = TRUE )
r2_LCP_gdtot <- r.squaredGLMM(LCP_gdtot_reml)
r2_LCP_gdtot

#with CS dist & 4 gen datasets:
#with LCP dist & 4 gen datasets:
CS_PS_reml <- lmer(PS~CS +EUC + (1|site_1), data = df_data_scaled, REML = TRUE )
r2_CS_PS <- r.squaredGLMM(CS_PS_reml)
r2_CS_PS
CS_aldiff_reml <- lmer(aldiff~CS +EUC + (1|site_1), data = df_data_scaled, REML = TRUE )
r2_CS_aldiff <- r.squaredGLMM(CS_aldiff_reml)
r2_CS_aldiff
CS_gdeuc_reml <- lmer(gdeuc~CS +EUC + (1|site_1), data = df_data_scaled, REML = TRUE )
r2_CS_gdeuc <- r.squaredGLMM(CS_gdeuc_reml)
r2_CS_gdeuc
CS_gdtot_reml <- lmer(gdtot~CS +EUC + (1|site_1), data = df_data_scaled, REML = TRUE )
r2_CS_gdtot <- r.squaredGLMM(CS_gdtot_reml)
r2_CS_gdtot


