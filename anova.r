library(tibble)
library(tidyverse)
library(ggpubr)
library(rstatix)

path_data<-"./"
name_data_file <-"SummaryResults.txt"
data_file<-paste(path_data, name_data_file,sep="")

# read the table it seems that header = TRUE does not set the colunm name but with tibble, it works
mydata <- tibble(read.table(data_file, header = TRUE))
data_anova<-mydata[,c("Obj","Slack", "Demand", "Capacity")]
data_anova <- data_anova %>% mutate(Demand=as.factor(Demand), Capacity=as.factor(Capacity), Slack=as.factor(Slack))

data_anova  %>%  group_by(Demand, Capacity,Slack) %>%
  get_summary_stats(Obj,type = "mean_sd")

bxp <- ggboxplot(
  data_anova, x = "Demand", y = "Obj",   color = "Capacity", palette = "jco", facet.by = "Slack"
  )
bxp

data_anova %>%
  group_by(Slack) %>%
  identify_outliers(Obj)

model  <- lm(Obj ~ Demand*Capacity*Slack, data = data_anova)
# Créer un QQ plot des résidus
ggqqplot(residuals(model))
# Calculer le test de normalité de Shapiro-Wilk
shapiro_test(residuals(model))

res.anova <-data_anova %>% anova_test(Obj ~ Slack*Demand*Capacity)
res.anova


