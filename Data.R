# Install and load the package
install.packages("readxl")
install.packages("lme4")
install.packages("lmerTest")
install.packages("ggplot2")
install.packages("ggeffects")
install.packages("mvabund")
install.packages("vegan")
install.packages("stringr")
install.packages("patchwork")
library(readxl)
library(lme4)
library(lmerTest)
library(ggplot2)
library(ggeffects)
library(mvabund)
library(vegan)
library(stringr)
library(patchwork)

# read data
data <- read_excel("/Users/yehong/Desktop/Data.xlsx", sheet = "Sheet1")
# Converts the field name
names(data) <- make.names(names(data))
# Converting Site variables to factor types
data$Site <- as.factor(data$Site)
# Clean up and convert SOM and Bulk.Density variables to numeric values
data$SOM <- as.numeric(gsub("[^0-9.]", "", data$SOM))
data$Bulk.Density <- as.numeric(gsub("[^0-9.]", "", data$Bulk.Density))
# Check data
str(data)


# Analyse SOM and Bulk Density
# Separate SOM and Bulk Density models
model_SOM <- lmer(SOM ~ Farmlet * Site + (1 | Field), data = data)
summary(model_SOM)
predicted_SOM <- ggpredict(model_SOM, terms = c("Farmlet", "Site"))
plot(predicted_SOM, show_data = T, jitter = T) +   labs(y = "Soil Organic Matter (%)", x = "Farmlet") +   theme_minimal() +   theme(plot.title = element_blank())

model_BulkDensity <- lmer(Bulk.Density ~ Farmlet * Site + (1 | Field), data = data)
summary(model_BulkDensity)
predicted_BulkDensity <- ggpredict(model_BulkDensity, terms = c("Farmlet", "Site"))
plot(predicted_BulkDensity, show_data = T, jitter = T) +   labs(y = expression(Bulk~Density~(g/cm^3)), x = "Farmlet") +   theme_minimal() +   theme(plot.title = element_blank())

# Standardise Bulk Density and rebuild the SOM model
data$Bulk.Density.z <- scale(data$Bulk.Density)[,1]
model_SOM1 <- lmer(SOM ~ Farmlet * Site * Bulk.Density.z + (1 | Field), data = data)
summary(model_SOM1)
anova(model_SOM1)

# Choose the best SOM model
model_SOM2 <- update(model_SOM1, .~.-Farmlet:Site:Bulk.Density.z)
anova(model_SOM2)
model_SOM3 <- update(model_SOM2, .~.-Site:Bulk.Density.z)
anova(model_SOM3)
model_SOM4 <- update(model_SOM3, .~.-Farmlet:Bulk.Density.z)
anova(model_SOM4)
model_SOM5 <- update(model_SOM4, .~.-Farmlet:Site)
anova(model_SOM5)
summary(model_SOM5)

# Plotting SOM model predictions for standardised BDs
plot(ggpredict(model_SOM5, terms = c("Bulk.Density.z", "Site", "Farmlet")), add.data = TRUE) +
  labs(title = "", x = "Z-Standardised Bulk Density", y = "SOM") +
  theme_minimal()

# Analyse Earthworm data
# Prepare multivariate response data
abund <- mvabund(data[, c("Juvenile_Count", "Epigeic_Count", "Anecic_Count", "Endogeic_Count")])
# Fit a multivariable generalized linear model
model_abund <- manyglm(abund ~ Farmlet + Site, data = data, family = "negative.binomial")
# plot a residual graph
plot(model_abund, which = 1)
summary(model_abund)

# Perform multivariate analysis of variance
anova.manyglm(model_abund, p.uni = "adjusted")

# Modelling - population analysis of juvenile earthworms
model_juvenile <- glmer(Juvenile_Count ~ Farmlet * Site + (1|Field), data = data, family = "poisson")
summary(model_juvenile)

model_juvenile_nb <- glmer.nb(Juvenile_Count ~ Farmlet + Site + (1|Field), data = data)
summary(model_juvenile_nb)
anova(model_juvenile_nb)

predicted_model_juvenile_nb <- ggpredict(model_juvenile_nb, terms = c("Farmlet", "Site"))
plot1 <-plot(predicted_model_juvenile_nb, show_data = T, jitter = T) +   labs(y = "Juvenile Count", x = "Farmlet") +   theme_minimal() +   theme(plot.title = element_blank())

# Modelling - Adult Earthworm Population Analysis
model_adult <- glmer(Adult_Count ~ Farmlet * Site + (1|Field), data = data, family = "poisson")
summary(model_adult)

model_adult_nb <- glmer.nb(Adult_Count ~ Farmlet + Site + (1|Field), data = data)
summary(model_adult_nb)
anova(model_adult_nb)

predicted_model_adult_nb <- ggpredict(model_adult_nb, terms = c("Farmlet", "Site"))
plot2 <-plot(predicted_model_adult_nb, show_data = T, jitter = T) +   labs(y = "Adult Count", x = "Farmlet") +   theme_minimal() +   theme(plot.title = element_blank())

# combined plot
combined_plot <- plot1 + plot2
print(combined_plot)

# Community Structure Analysis (NMDS)
nmds_abund <- metaMDS(abund, distance = "manhattan")

# Plotting the NMDS of Farmlet
farmlet_plot <- ordiplot(nmds_abund, type = "n", main = "")
orditorp(nmds_abund, display = "sites", col = ifelse(data$Farmlet == "Blue", "blue", "darkgreen"), air = 0.01)
ordiellipse(nmds_abund, data$Farmlet, col = c("blue", "darkgreen"))
legend("topright", legend = c("Blue Farmlet", "Green Farmlet"), col = c("blue", "darkgreen"), pch = 19, bty = "n")

# Plotting the NMDS of Site
site_plot <- ordiplot(nmds_abund, type = "n", main = "")
orditorp(nmds_abund, display = "sites", col = ifelse(data$Site == "Edge", "red", "black"), air = 0.01)
ordiellipse(nmds_abund, data$Site, col = c("red", "black"))
legend("topright", legend = c("Edge Site", "Middle Site"), col = c("red", "black"), pch = 19, bty = "n")

# Calculate relative proportions
data$Epigeic_proportion <- (data$Epigeic_Count / data$Adult_Count) * (data$Epigeic_Weight / data$Epigeic_Count)
data$Epigeic_proportion[is.na(data$Epigeic_proportion)] <- 0

data$Endogeic_proportion <- (data$Endogeic_Count / data$Adult_Count) * (data$Endogeic_Weight / data$Endogeic_Count)
data$Endogeic_proportion[is.na(data$Endogeic_proportion)] <- 0

data$Anecic_proportion <- (data$Anecic_Count / data$Adult_Count) * (data$Anecic_Weight / data$Anecic_Count)
data$Anecic_proportion[is.na(data$Anecic_proportion)] <- 0

data$cwm <- rowSums(data[, c("Epigeic_proportion", "Endogeic_proportion", "Anecic_proportion")])

# Modelling and plotting weighted mean body weights of communities
model_cwm <- lmer(cwm ~ Farmlet * Site + (1|Field), data = data)
summary(model_cwm)
anova(model_cwm)

# Analyzing the relationship between SOM, Bulk Density, and Earthworms.
model_SOM_Earthworm <- lmer(SOM ~ Farmlet * Site * Bulk.Density.z + Juvenile_Count + Adult_Count + (1|Field), data = data)
summary(model_SOM_Earthworm)
anova(model_SOM_Earthworm)

# Plot the relationship between SOM and earthworm populations.
plot(ggpredict(model_SOM_Earthworm, terms = c("Bulk.Density.z", "Juvenile_Count", "Adult_Count")), show_data = TRUE) +
  labs(title = "", x = "Z-Standardised Bulk Density", y = "SOM") +
  theme_minimal()

# summary
summary(model_juvenile)
summary(model_juvenile_nb)
summary(model_adult)
summary(model_adult_nb)
summary(model_cwm)

# Read data
compaction_data <- read_excel("/Users/yehong/Desktop/Compaction.xlsx")

# Convert field names
names(compaction_data) <- make.names(names(compaction_data))

# Converts related variables to factor types
compaction_data$Site <- as.factor(compaction_data$Site)
compaction_data$Field <- as.factor(compaction_data$Field)
compaction_data$Farmlet <- as.factor(compaction_data$Farmlet)
compaction_data$Site2 <- ifelse(str_detect(compaction_data$Site, "E"), "Edge", "Middle")
compaction_data$Depth.z <- scale(compaction_data$Depth)[,1]

# check data
str(compaction_data)
head(compaction_data)

# grouped data
edge_data <- subset(compaction_data, Site2 == "Edge")
middle_data <- subset(compaction_data, Site2 == "Middle")

# The relationship between pressure and depth
model <- lmer(Pressure ~ Depth * Site2 + (1 | Field), data = compaction_data)
summary(model)
plot(ggpredict(model, terms = c("Depth", "Site2")), add.data = TRUE, jitter = TRUE) +
  labs(title = "", x = "Depth(cm)", y = "Pressure(MPa)") +
  theme_minimal()

mean(compaction_data$Depth)
plot(model)

# Plot scatter plots
ggplot(compaction_data, aes(x = Depth, y = Pressure)) +
  geom_point(color = 'blue') +
  labs(title = "Depth vs Pressure", x = "Depth (cm)", y = "Pressure(MPa)") +
  theme_minimal()

# Calculate the descriptive statistics of the SOM.
mean_SOM <- mean(data$SOM, na.rm = TRUE)  
sd_SOM <- sd(data$SOM, na.rm = TRUE) 
min_SOM <- min(data$SOM, na.rm = TRUE)  
max_SOM <- max(data$SOM, na.rm = TRUE)  
median_SOM <- median(data$SOM, na.rm = TRUE)  

# Calculate the descriptive statistics of the BD.
mean_BD <- mean(data$Bulk.Density, na.rm = TRUE)
sd_BD <- sd(data$Bulk.Density, na.rm = TRUE)
min_BD <- min(data$Bulk.Density, na.rm = TRUE)
max_BD <- max(data$Bulk.Density, na.rm = TRUE)
median_BD <- median(data$Bulk.Density, na.rm = TRUE)

# Calculate the descriptive statistics of the Earthworm.
mean_juvenile <- mean(data$Juvenile_Count, na.rm = TRUE)
sd_juvenile <- sd(data$Juvenile_Count, na.rm = TRUE)
min_juvenile <- min(data$Juvenile_Count, na.rm = TRUE)
max_juvenile <- max(data$Juvenile_Count, na.rm = TRUE)
median_juvenile <- median(data$Juvenile_Count, na.rm = TRUE)

mean_adult <- mean(data$Adult_Count, na.rm = TRUE)
sd_adult <- sd(data$Adult_Count, na.rm = TRUE)
min_adult <- min(data$Adult_Count, na.rm = TRUE)
max_adult <- max(data$Adult_Count, na.rm = TRUE)
median_adult <- median(data$Adult_Count, na.rm = TRUE)

# output result
mean_SOM
sd_SOM
min_SOM
max_SOM
median_SOM

mean_BD
sd_BD
min_BD
max_BD
median_BD

mean_juvenile
sd_juvenile
min_juvenile
max_juvenile
median_juvenile

mean_adult
sd_adult
min_adult
max_adult
median_adult
