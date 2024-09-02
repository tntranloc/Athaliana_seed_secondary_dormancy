### This workflow allows us to run GLM with kinship as covariate
### which is not possible with base R glm function 
## required packages: Matrix, lme4, lme4qtl


# 'dorm01' is  data frame, 'germ_rate' is the response variable,
# 'bio9' and 'treatment' are fixed effects, and 'ID' corresponds to individual identifiers
# 'kin01_pd' is kinship matrix
# example kinship matrix
library(Matrix)
kin01 = read.csv("/path/dorm01_kinship_for_glm.csv", header = T)
kin02 = read.csv("/path/dorm02_kinship_for_glm.csv", header = T)
kin03 = read.csv("/path/dorm03_kinship_for_glm.csv", header = T)


## some cleaning up
rownames(kin01) = kin01$X
kin01 = kin01[,-1]
colnames(kin01) = rownames(kin01)
kin01 = as.matrix(kin01)
kin01_pd = nearPD(kin01, keepDiag = T, maxit = 1000, ensureSymmetry = T, conv.tol = 1e-6)$mat

rownames(kin02) = kin02$X
kin02 = kin02[,-1]
colnames(kin02) = rownames(kin02)
kin02 = as.matrix(kin02)
kin02_pd = nearPD(kin02, keepDiag = T, maxit = 1000, ensureSymmetry = T, conv.tol = 1e-6)$mat

rownames(kin03) = kin03$X
kin03 = kin03[,-1]
colnames(kin03) = rownames(kin03)
kin03 = as.matrix(kin03)
kin03_pd = nearPD(kin03, keepDiag = T, maxit = 1000, ensureSymmetry = T, conv.tol = 1e-6)$mat

## create matrices for the next step 
dorm01$treatment = as.factor(dorm01$treatment)
dorm01 = na.omit(dorm01)
summary(dorm01)
dorm01 = subset(dorm01, dorm01$ID %in% rownames(kin01))
mat1 = data.frame(dorm01$germ, dorm01$non_germ)
head(mat1)
colnames(mat1) = c("germ", "nongerm")
mat1 = as.matrix(mat1)

dorm02$treatment = as.factor(dorm02$treatment)
dorm02 = na.omit(dorm02)
summary(dorm02)
dim(dorm02)
dorm02 = dorm02[!duplicated(dorm02$ID),]
dorm02 = subset(dorm02, dorm02$ID %in% rownames(kin02))
mat2 = data.frame(dorm02$germ, dorm02$non_germ)
head(mat2)
colnames(mat2) = c("germ", "nongerm")
mat2 = as.matrix(mat2)

dorm03$treatment = as.factor(dorm03$treatment)
dorm03 = na.omit(dorm03)
summary(dorm03)
dorm03 = subset(dorm03, dorm03$ID %in% rownames(kin03))
mat3 = data.frame(dorm03$germ, dorm03$non_germ)
head(mat3)
colnames(mat3) = c("germ", "nongerm")
mat3 = as.matrix(mat3)

#############################################
## Load required libraries
library(lme4)
library(lme4qtl)


# Include kinship as a random effect in the model
model_bio4.2 = relmatGlmer(
  mat2 ~ treatment * scale(bio4) + (1 | ID),
  data = dorm02,
  family = binomial(link = "logit"),
  relmat = list(ID = kin02_pd)
)

summary(model_bio4.2)

### If things work fine, run the loop ####
### Initialize a list to store models
models_list1 = list()
# Iterate over columns from 8 to 26 where predictor variables are
for (i in 8:26) {
  # Extract the current bio column name
  bio_column = names(dorm02)[i]
  
  # Construct the model formula dynamically
  formula = as.formula(paste("mat1 ~ treatment * (", bio_column, ") + (1 | ID)", sep = ""))
  
  # Fit the model using relmatGlmer
  model = relmatGlmer(
    formula, 
    data = dorm01,
    family = binomial(link = "logit"),
    relmat = list(ID = kin01_pd)
  )
  
  # Save the model into the list using the bio column name
  model_name = paste("model_", bio_column, ".1", sep = "")
  models_list1[[model_name]] = model
}
summary(models_list1[[1]]) ## access model list by this

#############################################
### Tidy and export model result ####
### try one model first
library(broom.mixed)
tidy_fixed = broom.mixed::tidy(models_list1[[1]]) #fixed effects and random effect
tidy_fixed = broom.mixed::tidy(models_list1[[1]], effects = "ran_pars") #only random effects

print(tidy_fixed)

### Export all statistics at once from model list ###
library(broom.mixed)
library(dplyr)

# List to hold each tidy summary
tidy_summaries_1 = list()
tidy_summaries_2 = list()
tidy_summaries_3 = list()

# Loop through each model in the list
for (i in 1:length(models_list1)) {
    # Tidy the i-th model
    tidy_model = tidy(models_list1[[i]])
    # Add a column to identify the model
    tidy_model$model_index = paste("Model", i)
    # Append to the list
    tidy_summaries_1[[i]] = tidy_model
}

# Combine all summaries into one dataframe
all_summaries_1 = bind_rows(tidy_summaries_1)

# Write the combined dataframe to a CSV file
write.csv(all_summaries_1, "all_model_summaries_dorm01.csv", row.names = FALSE)


#############################################
######### Plot the odd ratios  #########

tidy_summaries_1[[8]]

plot1.8 = ggplot(tidy_summaries_1[[8]], aes(x = term, y = exp(estimate), ymin = exp(estimate - 1.96 * std.error), ymax = exp(estimate + 1.96 * std.error))) +
  geom_pointrange() +
  geom_hline(yintercept = 1, linetype = "dashed") +
  scale_y_log10() +  # This sets the scale to log10 for easier interpretation of odds ratios
  labs(y = "Odds Ratios", x = "") +
  coord_flip()  # This flips the coordinates so terms are on the y-axis and easier to read

######################################################
########## Plot fitted curve #########

#### Predicted probability of germination # germination is the response variable here
# Create a data frame for predictions
new_data = expand.grid(
  treatment = levels(dorm02$treatment),
  bio18 = seq(min(dorm02$bio18), max(dorm02$bio18), length.out = 100), 
  ID = dorm02$ID
)

# Add predictions to the data frame
new_data$predicted = predict(models_list2[[18]], new_data, type = "response")

#### Plot the predicted probabilities

  ## chose the palette
library(RColorBrewer)
display.brewer.pal(n = 8, name = "Dark2") # personal favourite and colourblind friendly
brewer.pal(n = 8, name = "Dark2")

p1 = ggplot(new_data, aes(x = bio18, y = predicted, color = treatment)) + #geom_point(alpha = 0.5) +
  geom_smooth(method = "glm", method.args = list(family = "binomial"), se = FALSE) +
  #geom_line() +
  labs(
    title = "Predicted Germination Probability by Treatment and Bio18",
    x = "BIO18",
    y = "Predicted Probability of Germination"
  ) +
  scale_color_manual(values = c("#1B9E77","#7570B3","#D95F02")) +
  theme(plot.title=element_text(family = "Lucida Grande", face = "bold", size=13,hjust=0.5,vjust=1), 
        legend.position = "right", 
        legend.title = element_text(colour = "black", size = 13, family = "Lucida Grande"),
        legend.text = element_text(colour = "black", size = 11, family = "Lucida Grande"),
        axis.title.y = element_text(colour = 'black',  size = 13, family = "Lucida Grande", angle=90), 
        axis.title.x = element_text(colour = "black", size = 13, family = "Lucida Grande"),
        axis.text.x = element_text(colour = "black", size = 13, family = "Lucida Grande"),
        axis.text.y = element_text(colour = "black", size = 13, family = "Lucida Grande"),)

