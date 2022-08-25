####### R Code for Summarising Accelerometer Measured Physical #######
#######    Activity Profiles for Modelling Health Outcomes     #######

#######         Dissertation Presented for the Degree of       #######
#######        MSc in Statistics and Operational Research      #######
#######          by  Minqing Li  (Maeve Li)   s2167017         #######

# Attach libraries and define used in this thesis
library(nhanesaccel)
library(SASxport)
library(foreign)
library(data.table)
library(dplyr)
library(VIM)
library(deltacomp)
library(compositions)
library(stringr)
library(ggplot2)
library(ggpubr)
library(gridExtra)
library(zCompositions)
library(psych)
library(boot)
library(ggfortify)
library(zoo)
library(mgcv)
library(Hmisc)
library(corrplot)
library(mgcViz)

se <- function(x) sqrt(var(x)/length(x))

############### Derive NHANES  data ############### 
# derive nhanes1 data using the package
nhanes1 <- process_nhanes(waves = 2,  # 05-06 data
                          brevity = 2, # include 5 intensity levels
                          valid_days = 5, # minimum valid days
                          nonwear_window = 60, 
                          # >60 is considered to be non-wear time
                          weartime_minimum = 600, 
                          # daily wear time minimum 10 hours
                          int_cuts = c(100,760,1952,5999))
# int_cuts is a numeric vector with four cutpoints from 
# which five intensity ranges are derived.
# intensity 1-5 are typically viewed as sedentary, light, 
# lifestyle, moderate, and vigorous.
# as we only study sendentary, LPA and MVPA, we can add up 
# light & lifestyle, and moderate & vigorous
# which is already done by the function (in column lightlife and mvpa)

mycolumns = c("seqn","include","valid_min","cpm",
              "sed_min","lightlife_min","mvpa_min")
# nhanes 1.1 eligible for analysis data & wanted columns
nhanes1.1 <- data.table(nhanes1) %>% 
  .[, ..mycolumns]
colnames(nhanes1.1)[6] <- "light_min"

# Preprocess raw data provided by the package
data("w2", envir = environment())
data("wave2_paxinten", envir = environment())

## remove values that are larger than 15000
length.wave2 <- length(wave2_paxinten)
for (i in 1:length.wave2){
  if (wave2_paxinten[i] >= 15000)
    wave2_paxinten[i] <- NA
} 

## Gives two examples of accelerometer PA data
start31199 <- w2[which(w2[,1]==31199, arr.ind = TRUE),2]
end31199 <- w2[which(w2[,1]==31199, arr.ind = TRUE),3]
acce31199 <- wave2_paxinten[start31199:end31199]

start31223 <- w2[which(w2[,1]==31223, arr.ind = TRUE),2]
end31223 <- w2[which(w2[,1]==31223, arr.ind = TRUE),3]
acce31223 <- wave2_paxinten[start31223:end31223]

d31199 <- data.frame(time = 1:length(acce31199), intensity = acce31199)
plot31199 <- ggplot(d31199, aes(time, intensity)) + geom_line() +
  ggtitle("PA profile for participant No.31199") +
  theme(plot.title = element_text(hjust = 0.5)) +
  scale_x_continuous("Day of the week", 
                     breaks=c(seq(1000,10080,1450)), 
                     labels=c("Sun", "Mon", "Tue", "Wed", "Thu", "Fri", "Sat")) +
  annotate("text", x=8000, y=12000, 
           label= "Female, age 25") + 
  annotate("text", x=8000, y=11000, 
           label= "CPM=681.33, BMI=24.20") + 
  ylab("intensity (count per minute)") +
  geom_hline(yintercept = 100, linetype="dotted", 
             color = "grey", size=1) +
  geom_hline(yintercept = 1951, linetype="dotted", 
             color = "grey", size=1) +
  geom_hline(yintercept = 5724, linetype="dotted", 
             color = "grey", size=1) +
  theme_bw()

d31223 <- data.frame(time = 1:length(acce31223), intensity = acce31223)
plot31223 <- ggplot(d31223, aes(time, intensity)) + geom_line() +
  ggtitle("PA profile for participant No.31223") +
  theme(plot.title = element_text(hjust = 0.5)) +
  scale_x_continuous("Day of the week", 
                     breaks=c(seq(1000,10080,1450)), 
                     labels=c("Mon", "Tue", "Wed", "Thu", "Fri", "Sat", "Sun")) +
  annotate("text", x=8000, y=12000, 
           label= "Male, age 61") + 
  annotate("text", x=8000, y=11000, 
           label= "CPM=457.78, BMI=30.22") +
  ylab("intensity (count per minute)") +
  geom_hline(yintercept = 100, linetype="dotted", 
             color = "grey", size=1) +
  geom_hline(yintercept = 1951, linetype="dotted", 
             color = "grey", size=1) +
  geom_hline(yintercept = 5724, linetype="dotted", 
             color = "grey", size=1) +
  theme_bw()

ggarrange(plot31199, plot31223, ncol=2, nrow=1)

############### Merge BMI, sleep and Demographic Data ###############
BMX <- read.xport("E:/Mae/EDIN/Dissertation/Data/NHANES/BMX_D.XPT")
DEMO <- read.xport("E:/Mae/EDIN/Dissertation/Data/NHANES/DEMO_D.XPT")
SLQ <- read.xport("E:/Mae/EDIN/Dissertation/Data/NHANES/SLQ_D.XPT")

Sleep <-SLQ[,c(1,2)] %>% data.table()
colnames(Sleep) <- c("seqn","sleep")
BMXBMI <-BMX[,c(1,11)] %>% data.table()
colnames(BMXBMI) <- c("seqn","BMI")
Demographics <- DEMO[,c("SEQN","RIDAGEYR","RIAGENDR","RIDRETH1",
                        "INDHHINC","DMDEDUC2")] %>% data.table()
colnames(Demographics) <- c("seqn","age","gender","race","income","education")

DEMO.common <- Demographics[Demographics$seqn %in% nhanes1.1$seqn,]
all.equal(DEMO.common$seqn,nhanes1.1$seqn)
BMI.common <- BMXBMI[BMXBMI$seqn %in% nhanes1.1$seqn,]
all.equal(BMI.common$seqn,nhanes1.1$seqn) 
Sleep.common <- Sleep[Sleep$seqn %in% nhanes1.1$seqn,]
all.equal(Sleep.common$seqn,nhanes1.1$seqn) # not equal (only >=16 have sleep data)

nhanes1.2 <- merge(BMI.common, nhanes1.1, by="seqn", all=TRUE) %>%
  merge(DEMO.common, ., by="seqn", all=TRUE) %>% 
  merge(Sleep.common, ., by="seqn", all=TRUE) %>%
  data.table()

nhanes1.21 <- nhanes1.2[include==1,] # 4339 wear-time eligible

# target group: >=20 yrs adults
nhanes1.3 <- nhanes1.21[age>=20,] # 2754 adults

# Change numbers that indicates missing values to NA
nhanes1.3 <- nhanes1.3[, income:=ifelse(c(income==77|income==99),NA,income)] %>%
  .[,education:=ifelse(c(education==7|education==9),NA,education)] %>%
  .[,sleep:=ifelse(c(sleep==77|sleep==99),NA,sleep)]

# Inspect Missingness
v <- aggr(nhanes1.3, plot = FALSE)
plot(v, col=c("mistyrose","indianred2"), prop = FALSE, numbers = TRUE)

# Keep complete rows
nhanes2 <- nhanes1.3[complete.cases(nhanes1.3),] # 2664 participants left

# Change numbers to factors
nhanes2 <- nhanes2[,income:=ifelse(c(income==1|income==2|income==3|income==4
                                     |income==13),
                                   "Low","High")] %>%
  .[,gender:=ifelse(gender==1,"Male","Female")] %>%
  .[,race:=ifelse(race==1,"Mexican American",
                  ifelse(race==2,"Other Hispanic",
                         ifelse(race==3,"Non-Hispanic White",
                                ifelse(race==4,"Non-Hispanic Black",
                                       "Other Races"))))] %>%
  .[,education:=ifelse(education==1|education==2,"No Highschool Diploma",
                       ifelse(education==3,"Highschool Grad",
                              "College Graduate or Above"))]
nhanes2$gender <- as.factor(nhanes2$gender)
nhanes2$race <- as.factor(nhanes2$race)
nhanes2$education <- as.factor(nhanes2$education)
nhanes2$income <- as.factor(nhanes2$income)

# Impute Zeros
nhanes2.subset <- lrEM(nhanes2[,c(10:14)],label=0,dl=rep(1,5),
                       ini.cov="multRepl") %>% data.table()
nhanes2[,mvpa_min:=nhanes2.subset$mvpa_min]

# Recalculate proportions of time use behaviour
nhanes3 <- nhanes2 %>% copy() %>% 
  .[,valid_min:=sed_min+light_min+mvpa_min] %>%
  .[,sed_percent:=sed_min/valid_min] %>%
  .[,light_percent:=light_min/valid_min] %>%
  .[,mvpa_percent:=mvpa_min/valid_min] %>%
  .[,sleep:=60*sleep] %>%  # transfer into minutes
  .[,sleep_percent:=sleep/1440] %>%
  .[,sed_percent:=sed_percent*(1-sleep_percent)] %>%
  .[,light_percent:=light_percent*(1-sleep_percent)] %>%
  .[,mvpa_percent:=mvpa_percent*(1-sleep_percent)]
all.equal(rep(1,2664),nhanes3$light_percent+nhanes3$sed_percent+
            nhanes3$mvpa_percent+nhanes3$sleep_percent)

### Correlation of edu with income
income_edu <- table(nhanes3$education, nhanes3$income)
chisqres <- chisq.test(income_edu)
chisqres 
# Crammer's V for nominal vars
sqrt(chisqres$statistic/sum(income_edu)*min(dim(income_edu)-1))

############### Merge for Descriptive Stats ###############

### 1. Inspecting the whole population ###

# did not merge education here because children and adults have 
# different vars for edu and also children's education level largely 
# increase with age (since they're still in school)
nhanes.des <- merge(BMXBMI, nhanes1.1, by="seqn", all=TRUE) %>%
  .[complete.cases(.),] %>% data.table() %>%
  .[include==1,]
DEMO.common.des <- DEMO[DEMO$SEQN %in% nhanes.des$seqn,]
all.equal(DEMO.common.des[,1],nhanes.des$seqn)
nhanes.des <- nhanes.des[,age:=DEMO.common.des$RIDAGEYR] %>%
  .[,gender:=DEMO.common.des$RIAGENDR] %>%
  .[,gender:=ifelse(gender==1,"Male","Female")] %>%
  .[,agegroup:=ifelse(age<12,"Children",ifelse(age<20,"Youth","Adult"))] %>%
  .[,race:=DEMO.common.des$RIDRETH1] %>%
  .[,race:=ifelse(race==1,"Mexican American",
                  ifelse(race==2,"Other Hispanic",
                         ifelse(race==3,"Non-Hispanic White",
                                ifelse(race==4,"Non-Hispanic Black",
                                       "Other Races"))))] %>%
  .[,income:=DEMO.common.des$INDHHINC] %>%
  .[,income:=ifelse(c(income==1|income==2|income==3|income==4|income==13),
                    "Low","High")] %>%
  .[complete.cases(.),]

# for better text presentation in the boxplot
nhanes.des$race <- str_wrap(nhanes.des$race, width = 13) 
nhanes.des$agegroup <- factor(nhanes.des$agegroup,levels=c("Children",
                                                           "Youth","Adult"))

### boxplots
## BMI
# gender
bmi1 <- ggplot(nhanes.des, aes(x=gender, y=BMI)) +
  geom_boxplot(color="#66C2A5", fill="#66C2A5", alpha=0.2, 
               size=1, width = 0.6, fatten = 1) +
  ggtitle("BMI stratified by gender") +
  theme(plot.title = element_text(hjust = 0.5)) +
  theme(legend.position="bottom") 
# age
bmi2 <- ggplot(nhanes.des, aes(x=agegroup, y=BMI)) +
  geom_boxplot(color="#FC8D62", fill="#FC8D62", alpha=0.2,  
               size=1, width = 0.6, fatten = 1) +
  ggtitle("BMI stratified by age group") +
  theme(plot.title = element_text(hjust = 0.5))+
  theme(legend.position="bottom") 
# race
bmi3 <- ggplot(nhanes.des, aes(x=race, y=BMI)) +
  geom_boxplot(color="#56B4E9", fill="#56B4E9", alpha=0.2, 
               size=1, width = 0.6, fatten = 1) +
  ggtitle("BMI stratified by race") +
  theme(plot.title = element_text(hjust = 0.5)) +
  theme(legend.position="right") 
# income
bmi4 <- ggplot(nhanes.des, aes(x=income, y=BMI)) +
  geom_boxplot(color="#E78AC3", fill="#E78AC3", alpha=0.2, 
               size=1, width = 0.6, fatten = 1) +
  ggtitle("BMI stratified by income") +
  theme(plot.title = element_text(hjust = 0.5)) +
  theme(legend.position="bottom") 
plt1 <- ggarrange(bmi1, bmi2, bmi3, bmi4, ncol=2, nrow=2)
plt1

# Mean count per minute
cpm1 <- ggplot(nhanes.des, aes(x=gender, y=cpm)) +
  geom_boxplot(color="#66C2A5", fill="#66C2A5", alpha=0.2, 
               size=1, width = 0.6, fatten = 1) +
  ggtitle("CPM stratified by gender") +
  theme(plot.title = element_text(hjust = 0.5)) +
  theme(legend.position="bottom") 

cpm2 <- ggplot(nhanes.des, aes(x=agegroup, y=cpm)) +
  geom_boxplot(color="#FC8D62", fill="#FC8D62", alpha=0.2,
               size=1, width = 0.6, fatten = 1) +
  ggtitle("CPM stratified by age group") +
  theme(plot.title = element_text(hjust = 0.5)) +
  theme(legend.position="bottom") 

cpm3 <- ggplot(nhanes.des, aes(x=race, y=cpm)) +
  geom_boxplot(color="#56B4E9", fill="#56B4E9", alpha=0.2,
               size=1, width = 0.6, fatten = 1) +
  ggtitle("CPM stratified by race") +
  theme(plot.title = element_text(hjust = 0.5)) +
  theme(legend.position="bottom") 

cpm4 <- ggplot(nhanes.des, aes(x=income, y=cpm)) +
  geom_boxplot(color="#E78AC3", fill="#E78AC3", alpha=0.2,
               size=1, width = 0.6, fatten = 1) +
  ggtitle("CPM stratified by income") +
  theme(plot.title = element_text(hjust = 0.5)) +
  theme(legend.position="bottom") 
plt2 <- ggarrange(cpm1, cpm2, cpm3, cpm4, ncol=2, nrow=2)
plt2 

## histogram
## cpm
rects <- data.frame(ymin = -Inf, 
                    ymax = Inf,
                    xmin = c(0,100,1951),  
                    xmax = c(100,1951,Inf))

hist1 <- ggplot(nhanes.des, aes(x=cpm)) +
  xlab("Intensity of PA (cpm)")  +  
  ylab("Counts") +
  geom_rect(data = rects, aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax), 
            fill =  c("green", "DodgerBlue", "purple"),
            # Control the shading opacity here
            inherit.aes = FALSE, alpha = 0.11) +
  geom_histogram(color="black", fill="white",binwidth = 20) +
  ggtitle("CPM Histogram") +
  theme(plot.title = element_text(hjust = 0.5)) +
  scale_x_continuous( expand = c(0,0) , limits = c(0,2300) ) +
  scale_y_continuous( expand = c(0,0), limits = c(0,230) ) +
  geom_vline(xintercept = 100, linetype="dotted", 
             color = "grey", size=1) +
  geom_vline(xintercept = 1951, linetype="dotted", 
             color = "grey", size=1)

hist1

### 2. Inspecting adults that were excluded from the weartime requirements ###

nhanes0.0 <- nhanes1.2[include==0,] %>% copy() %>%
  .[age>=20, ] %>%
  .[,income:=NULL] %>%
  .[,education:=ifelse(c(education==7|education==9),NA,education)] %>%
  .[,sleep:=ifelse(c(sleep==77|sleep==99),NA,sleep)] %>%
  .[complete.cases(.),] %>%
  .[,gender:=ifelse(gender==1,"Male","Female")] %>%
  .[,race:=ifelse(race==1,"Mexican American",
                  ifelse(race==2,"Other Hispanic",
                         ifelse(race==3,"Non-Hispanic White",
                                ifelse(race==4,"Non-Hispanic Black",
                                       "Other Races"))))] %>%
  .[,education:=ifelse(education==1|education==2,"No Highschool Diploma",
                       ifelse(education==3,"Highschool Grad",
                              "College Graduate or Above"))]


nhanes0.subset <- lrEM(nhanes0.0[,c(9:13)],label=0,dl=rep(1,5),
                       ini.cov="multRepl") %>% data.table()

nhanes0.0 <- nhanes0.0[,mvpa_min:=nhanes0.subset$mvpa_min] %>% 
  .[,valid_min:=sed_min+light_min+mvpa_min] %>%
  .[,sed_percent:=sed_min/valid_min] %>%
  .[,light_percent:=light_min/valid_min] %>%
  .[,mvpa_percent:=mvpa_min/valid_min] %>%
  .[,sleep:=60*sleep] %>%  
  .[,sleep_percent:=sleep/1440] %>%
  .[,sed_percent:=sed_percent*(1-sleep_percent)] %>%
  .[,light_percent:=light_percent*(1-sleep_percent)] %>%
  .[,mvpa_percent:=mvpa_percent*(1-sleep_percent)]

nhanes0.0$gender <- as.factor(nhanes0.0$gender)
nhanes0.0$race <- as.factor(nhanes0.0$race)
nhanes0.0$education <- as.factor(nhanes0.0$education)
nhanes0.0$income <- as.factor(nhanes0.0$income)

summary(nhanes0.0) 

se(nhanes0.0$age)
se(nhanes0.0$BMI)

# calculate geometric means of PA variables
sleep.gm0 <- geometric.mean(nhanes0.0$sleep_percent)
sed.gm0 <- geometric.mean(nhanes0.0$sed_percent)
lpa.gm0 <- geometric.mean(nhanes0.0$light_percent)
mvpa.gm0 <- geometric.mean(nhanes0.0$mvpa_percent)

sleep.gm0.1 <- sleep.gm0/(sleep.gm0 + sed.gm0 + lpa.gm0 + mvpa.gm0)
sed.gm0.1 <- sed.gm0/(sleep.gm0 + sed.gm0 + lpa.gm0 + mvpa.gm0)
lpa.gm0.1 <- lpa.gm0/(sleep.gm0 + sed.gm0 + lpa.gm0 + mvpa.gm0)
mvpa.gm0.1 <- mvpa.gm0/(sleep.gm0 + sed.gm0 + lpa.gm0 + mvpa.gm0)

############### CoDA with age, gender, edu, race and income as covariates ###############
# derive compositional values
# study sleep
comp1 <- nhanes3 %>% copy() %>%
  .[,c("sed_percent","light_percent","mvpa_percent","sleep_percent")]  
# study sed
comp2 <- nhanes3 %>% copy() %>%
  .[,c("light_percent","mvpa_percent","sleep_percent","sed_percent")]
# study lpa
comp3 <- nhanes3 %>% copy() %>%
  .[,c("mvpa_percent","sleep_percent","sed_percent","light_percent")] 
# study mvpa
comp4 <- nhanes3 %>% copy() %>%
  .[,c("sleep_percent","sed_percent","light_percent","mvpa_percent")] 

# Create isometric log ratio coordinates
ilr1 <-ilr(comp1) %>% data.table() %>%
  .[,BMI:=nhanes3$BMI] %>%
  .[,age:=nhanes3$age] %>%
  .[,gender:=nhanes3$gender] %>%
  .[,race:=nhanes3$race] %>%
  .[,income:=nhanes3$income] %>%
  .[,education:=nhanes3$education]

ilr2 <-ilr(comp2) %>% data.table() %>%
  .[,BMI:=nhanes3$BMI] %>%
  .[,age:=nhanes3$age] %>%
  .[,gender:=nhanes3$gender] %>%
  .[,race:=nhanes3$race] %>%
  .[,income:=nhanes3$income] %>%
  .[,education:=nhanes3$education]

ilr3 <-ilr(comp3) %>% data.table() %>%
  .[,BMI:=nhanes3$BMI] %>%
  .[,age:=nhanes3$age] %>%
  .[,gender:=nhanes3$gender] %>%
  .[,race:=nhanes3$race] %>%
  .[,income:=nhanes3$income] %>%
  .[,education:=nhanes3$education]

ilr4 <-ilr(comp4) %>% data.table() %>%
  .[,BMI:=nhanes3$BMI] %>%
  .[,age:=nhanes3$age] %>%
  .[,gender:=nhanes3$gender] %>%
  .[,race:=nhanes3$race] %>%
  .[,income:=nhanes3$income] %>%
  .[,education:=nhanes3$education]

# Build CoDA model
# the last regression coefficient (V3) is of interest according to
# https://stats.stackexchange.com/questions/244118/
# how-to-use-isometric-logratio-ilr-from-a-package-compositions

# this can be confirmed by looking at the value obtained by the deltacomp package)

# studies sleep
coda1 <- lm(BMI~V1+V2+V3+
              age+gender+race+income+education,data=ilr1)
summary(coda1)
drop1(coda1, test="F")
par(mfrow = c(2, 2))
plot(coda1)

# studies SB
coda2 <- lm(BMI~V1+V2+V3+
              age+gender+race+income+education,data=ilr2)
summary(coda2)

# studies LPA
coda3 <- lm(BMI~V1+V2+V3+
              age+gender+race+income+education,data=ilr3)
summary(coda3) 

# studies MVPA
coda4 <- lm(BMI~V1+V2+V3+
              age+gender+race+income+education,data=ilr4)
summary(coda4) 

############### CoDA with income removed ###############
## Create income-removed dataset
nhanes1.31 <- nhanes1.3 %>% copy()%>% .[,income:=NULL]
nhanes2.1 <- nhanes1.31[complete.cases(nhanes1.31),]

nhanes2.1 <- nhanes2.1[,gender:=ifelse(gender==1,"Male","Female")] %>%
  .[,race:=ifelse(race==1,"Mexican American",
                  ifelse(race==2,"Other Hispanic",
                         ifelse(race==3,"Non-Hispanic White",
                                ifelse(race==4,"Non-Hispanic Black",
                                       "Other Races"))))] %>%
  .[,education:=ifelse(education==1|education==2,"No Highschool Diploma",
                       ifelse(education==3,"Highschool Grad",
                              "College Graduate or Above"))]
nhanes2.1$gender <- as.factor(nhanes2.1$gender)
nhanes2.1$race <- as.factor(nhanes2.1$race)
nhanes2.1$education <- as.factor(nhanes2.1$education)

# Impute Zeros
nhanes2.1.subset <- lrEM(nhanes2.1[,c(9:13)],label=0,dl=rep(1,5),
                         ini.cov="multRepl") %>% data.table()
nhanes2.1[,mvpa_min:=nhanes2.1.subset$mvpa_min]

# Recalculate proportions of time use behaviour
nhanes3.1 <- nhanes2.1 %>% copy() %>% 
  .[,valid_min:=sed_min+light_min+mvpa_min] %>%
  .[,sed_percent:=sed_min/valid_min] %>%
  .[,light_percent:=light_min/valid_min] %>%
  .[,mvpa_percent:=mvpa_min/valid_min] %>%
  .[,sleep:=60*sleep] %>%  # transfer into minutes
  .[,sleep_percent:=sleep/1440] %>%
  .[,sed_percent:=sed_percent*(1-sleep_percent)] %>%
  .[,light_percent:=light_percent*(1-sleep_percent)] %>%
  .[,mvpa_percent:=mvpa_percent*(1-sleep_percent)]
all.equal(rep(1,2738),nhanes3.1$light_percent+nhanes3.1$sed_percent+
            nhanes3.1$mvpa_percent+nhanes3.1$sleep_percent)

summary(nhanes3.1)

se(nhanes3.1$age)
se(nhanes3.1$BMI)

# calculate geometric means of PA variables
sleep.gm1 <- geometric.mean(nhanes3.1$sleep_percent)
sed.gm1 <- geometric.mean(nhanes3.1$sed_percent)
lpa.gm1 <- geometric.mean(nhanes3.1$light_percent)
mvpa.gm1 <- geometric.mean(nhanes3.1$mvpa_percent)

sleep.gm1.1 <- sleep.gm1/(sleep.gm1 + sed.gm1 + lpa.gm1 + mvpa.gm1)
sed.gm1.1 <- sed.gm1/(sleep.gm1 + sed.gm1 + lpa.gm1 + mvpa.gm1)
lpa.gm1.1 <- lpa.gm1/(sleep.gm1 + sed.gm1 + lpa.gm1 + mvpa.gm1)
mvpa.gm1.1 <- mvpa.gm1/(sleep.gm1 + sed.gm1 + lpa.gm1 + mvpa.gm1)

## Refit CoDA model
# derive compositional values
# study sleep
comp1.1 <- nhanes3.1 %>% copy() %>%
  .[,c("sed_percent","light_percent","mvpa_percent","sleep_percent")] 
# study sed
comp2.1 <- nhanes3.1 %>% copy() %>%
  .[,c("light_percent","mvpa_percent","sleep_percent","sed_percent")] 
# study lpa
comp3.1 <- nhanes3.1 %>% copy() %>%
  .[,c("mvpa_percent","sleep_percent","sed_percent","light_percent")] 
# study mvpa
comp4.1 <- nhanes3.1 %>% copy() %>%
  .[,c("sleep_percent","sed_percent","light_percent","mvpa_percent")] 

# Create isometric log ratio coordinates
ilr1.1 <-ilr(comp1.1) %>% data.table() %>%
  .[,BMI:=nhanes3.1$BMI] %>%
  .[,age:=nhanes3.1$age] %>%
  .[,gender:=nhanes3.1$gender] %>%
  .[,race:=nhanes3.1$race] %>%
  .[,education:=nhanes3.1$education]

ilr2.1 <-ilr(comp2.1) %>% data.table() %>%
  .[,BMI:=nhanes3.1$BMI] %>%
  .[,age:=nhanes3.1$age] %>%
  .[,gender:=nhanes3.1$gender] %>%
  .[,race:=nhanes3.1$race] %>%
  .[,education:=nhanes3.1$education]

ilr3.1 <-ilr(comp3.1) %>% data.table() %>%
  .[,BMI:=nhanes3.1$BMI] %>%
  .[,age:=nhanes3.1$age] %>%
  .[,gender:=nhanes3.1$gender] %>%
  .[,race:=nhanes3.1$race] %>%
  .[,education:=nhanes3.1$education]

ilr4.1 <-ilr(comp4.1) %>% data.table() %>%
  .[,BMI:=nhanes3.1$BMI] %>%
  .[,age:=nhanes3.1$age] %>%
  .[,gender:=nhanes3.1$gender] %>%
  .[,race:=nhanes3.1$race] %>%
  .[,education:=nhanes3.1$education]

# Build CoDA model
# Studies Sleep
coda1.1 <- glm(BMI~V1+V2+V3+
                 age+gender+race+education,data=ilr1.1)
summary(coda1.1)
100*with(summary(coda1.1), 1 - deviance/null.deviance)

# Studies SB
coda2.1 <- lm(BMI~V1+V2+V3+
                age+gender+race+education,data=ilr2.1)
summary(coda2.1) 

# Studies LPA
coda3.1 <- lm(BMI~V1+V2+V3+
                age+gender+race+education,data=ilr3.1)
summary(coda3.1) 

# Studies MVPA
coda4.1 <- lm(BMI~V1+V2+V3+
                age+gender+race+education,data=ilr4.1)
summary(coda4.1) 

AIC(coda1.1)

## diagnosis & checking residuals
par(mfrow = c(2, 2))
plot(coda1.1)
autoplot(coda1.1, label.size = 3) + theme_bw()

par(mfrow = c(1, 1))
diag1 <- qplot(.fitted, .resid, data = coda1.1) +
  geom_hline(yintercept = 0) +
  geom_smooth(se = FALSE) +
  ggtitle("Residuals vs Fitted plot") +
  theme(plot.title = element_text(hjust = 0.5))

diag2 <- qplot(sample =.stdresid, data = coda1.1, stat = "qq") + 
  geom_abline() +
  xlab("theoretcial quantiles") +
  ylab("standardized residuals") +
  ggtitle("Q-Q plot") +
  theme(plot.title = element_text(hjust = 0.5))

ggarrange(diag1, diag2, ncol=2, nrow=1)

############### Using Deltacomp Package: Estimating differences and plotting changes ###############
# study sleep
pred_df1 <- predict_delta_comps(
  dataf = nhanes3.1,
  y = "BMI",
  comps = c("sleep_percent", "sed_percent", "light_percent", "mvpa_percent"),
  # only provide one regression model?
  covars = c("age", "gender","education","race"),
  deltas = seq(-12.5, 30, by = 2) / (24 * 60),
  comparisons = "prop-realloc",
  alpha = 0.05
)
# study sedentary
pred_df2 <- predict_delta_comps(
  dataf = nhanes3.1,
  y = "BMI",
  comps = c("sed_percent", "sleep_percent", "light_percent", "mvpa_percent"),
  # only provide one regression model?
  covars = c("age", "gender","education","race"),
  deltas = seq(-12.5, 30, by = 2) / (24 * 60),
  comparisons = "prop-realloc",
  alpha = 0.05
)
# study LPA
pred_df3 <- predict_delta_comps(
  dataf = nhanes3.1,
  y = "BMI",
  comps = c("light_percent", "sleep_percent", "sed_percent", "mvpa_percent"),
  # only provide one regression model?
  covars = c("age", "gender","education","race"),
  deltas = seq(-12.5, 30, by = 2) / (24 * 60),
  comparisons = "prop-realloc",
  alpha = 0.05
)
# study MVPA
pred_df4 <- predict_delta_comps(
  dataf = nhanes3.1,
  y = "BMI",
  comps = c("mvpa_percent", "sleep_percent", "sed_percent", "light_percent"),
  # only provide one regression model?
  covars = c("age", "gender","education","race"),
  deltas = seq(-12.5, 30, by = 2) / (24 * 60),
  comparisons = "prop-realloc",
  alpha = 0.05
)

# Plots
plot_delta_comp(
  pred_df1, # provide the returned object from predict_delta_comps()
  # x-axis can be converted from propotion of composition to meaningful units
  comp_total = 24 * 60, # minutes available in the composition
  units_lab = "min" # just a label for plotting
)

pred_df1_gg <- pred_df1 %>% copy() %>% data.table() %>%
  .[,Change_delta_in_composition_min:=delta*1440] 
names(pred_df1_gg)[1] <- "Physical_Activity"
pred_df1_gg[,Physical_Activity:=ifelse(Physical_Activity=="sleep_percent","Sleep",
                                       ifelse(Physical_Activity=="sed_percent",
                                              "SB",ifelse(Physical_Activity==
                                                       "light_percent",
                                                       "LPA","MVPA")))]
pred_df1_gg$Physical_Activity <- factor(pred_df1_gg$Physical_Activity,
                                        levels=c("Sleep","SB","LPA","MVPA"))

ggplot(pred_df1_gg, aes(x=Change_delta_in_composition_min,y=delta_pred,
                        color = Physical_Activity)) + 
  geom_point(aes(color = Physical_Activity), size = 1) +
  geom_line(aes(color = Physical_Activity, linetype = Physical_Activity), 
            size = 0.5) + 
  xlab("Change(delta) in minutes") + 
  ylab("Estimates change in BMI with 95% CI") +
  geom_vline(xintercept=0) + geom_hline(yintercept = 0) +
  geom_ribbon(aes(ymin = ci_lo, ymax = ci_up), size = 0.5,alpha = 0.1)+ 
  scale_x_continuous(breaks = seq(-12, 30, by = 2)) + 
  scale_y_continuous(breaks = seq(-2, 3, by=0.2)) +
  ggtitle("Predicted Changes in BMI When Components of PA Change Proportionally") +
  theme_light() + 
  theme(plot.title = element_text(hjust = 0.5))

# 1v1 change
pred_df5 <- predict_delta_comps(
  dataf = nhanes3.1,
  y = "BMI",
  comps = c("sleep_percent", "sed_percent", "light_percent", "mvpa_percent"),
  # only provide one regression model?
  covars = c("age", "gender","education","race"),
  deltas = seq(-12.5, 12.5, by = 2) / (24 * 60),
  comparisons = "one-v-one",
  alpha = 0.05
)

plot_delta_comp(
  pred_df5, # provide the returned object from predict_delta_comps()
  # x-axis can be converted from propotion of composition to meaningful units
  comp_total = 24 * 60, # minutes available in the composition
  units_lab = "min" # just a label for plotting
)


#### Baseline model ####
# derive time-use variables (min/day)
nhanes3.1[,SB:=1440*sed_percent] %>%
  .[,LPA:=1440*light_percent] %>%
  .[,MVPA:=1440*mvpa_percent]

all.equal(rep(1440,2738),nhanes3.1$LPA+nhanes3.1$SB+nhanes3.1$MVPA+
            nhanes3.1$sleep)

# Leave out mvpa
basemodel1 <- glm(BMI ~ sleep + SB + LPA + 
                   age + gender + race + education, data = nhanes3.1)
# Leave out lpa
basemodel2 <- lm(BMI ~ sleep + SB + MVPA + 
                   age + gender + race + education, data = nhanes3.1)
# Leave out sedentary
basemodel3 <- lm(BMI ~ sleep + MVPA + LPA + 
                   age + gender + race + education, data = nhanes3.1)
# Leave out Sleep
basemodel4 <- lm(BMI ~ SB + LPA + MVPA + 
                   age + gender + race + education, data = nhanes3.1)
summary(basemodel1)
summary(basemodel2)
summary(basemodel3)
summary(basemodel4)
AIC(basemodel2)

PAcorr <- rcorr(data.matrix(nhanes3.1[,c(2,18:20)]))
PAcorrp <- PAcorr$r
corrplot.mixed(PAcorrp, order = "original", 
               lower ="square", upper = "number",
               tl.cex = 1.2, tl.col = "black")

# calculate explained variance
100*with(summary(basemodel1), 1 - deviance/null.deviance)

############### CoDA searching for a new family distribution ###############

# gamma family with log link (studie sleep)
coda.gamma.log <- glm(BMI ~ V1 + V2 + V3 +
                        age + gender + race + education, 
                      family = "Gamma" (link = "log"), data=ilr1.1)
summary(coda.gamma.log)
100*with(summary(coda.gamma.log), 1 - deviance/null.deviance)
# diagnosis plots
coda.diag1 <- glm.diag(coda.gamma.log)
glm.diag.plots(coda.gamma.log, coda.diag1) # better fit, same as the histogram method
autoplot(coda.gamma.log, label.size = 3) + theme_bw()

# gaussian family with log link
coda.gaussian.log <- glm(BMI ~ V1 + V2 + V3 +
                           age + gender + race + education, 
                         family = "gaussian" (link = "log"), data=ilr1.1)
summary(coda.gaussian.log)
# diagnosis plots
coda.diag2 <- glm.diag(coda.gaussian.log)
glm.diag.plots(coda.gaussian.log, coda.diag2)

#### all four CoDA models using the gamma distribution and link function
## study SB
coda.gamma.log2 <- glm(BMI ~ V1 + V2 + V3 +
                         age + gender + race + education, 
                       family = "Gamma" (link = "log"), data=ilr2.1)
summary(coda.gamma.log2)

## study LPA
coda.gamma.log3 <- glm(BMI ~ V1 + V2 + V3 +
                         age + gender + race + education, 
                       family = "Gamma" (link = "log"), data=ilr3.1)
summary(coda.gamma.log3) 

## study MVPA 
coda.gamma.log4 <- glm(BMI ~ V1 + V2 + V3 +
                         age + gender + race + education, 
                       family = "Gamma" (link = "log"), data=ilr4.1)
summary(coda.gamma.log4)

#### Preprocessing Histogram Data ####

#### tests with cutpoints
test1 <- hist(wave2_paxinten,plot = FALSE)
test1[[1]] # see cutpoints
test1[[2]] # see counts

#### initial cut points
## For more conciseness in code, the bins of cut points can be directly changed
## here and rerun the code in Preprocessing Histogram Data Section 
test2 <- hist(wave2_paxinten, breaks=c(seq(0,10000,100),15000), 
              # 101 bins here
              plot=FALSE)
# best performed bins here
# test2 <- hist(wave2_paxinten, breaks=c(0,seq(100,10100,200),15000),plot=FALSE)
cuts.2 <- test2[[1]]

#### Build midpoints data frame (P) (same dimension as the frequency df)
midpoints <- data.frame(matrix(0, nrow=1, ncol=length(cuts.2)-1))
for (i in 1:(length(cuts.2)-1)){
  midpoints[1,i] <- (cuts.2[i+1] + cuts.2[i])/2
}
midpoints <- midpoints[rep(seq_len(nrow(midpoints)), each = 2738), ]

#### obtain common sequence numbers for the histogram method
w2.0 <- w2[w2[,1] %in% nhanes3.1$seqn,] %>% data.table()

## build histogram frequency matrix
histdf <- data.frame()
for (i in 1:nrow(w2.0)){ # for each participant
  acc_i <- wave2_paxinten[w2.0[i,start]:w2.0[i,stop]]
  
  # remove 60 consecutive zeros (non-wear time)
  acc_i <- c(na.omit(na.locf0(ifelse(acc_i == 0, NA, 0), maxgap = 60) + acc_i))
  
  hist_i <- hist(acc_i, breaks=cuts.2, plot=FALSE)
  num.bin <- length(cuts.2)-1
  for (j in 1:num.bin){
    histdf[i,j] <- hist_i[[2]][j]
  }
}
histmt <- data.matrix(histdf)
midpointsmt <- data.matrix(midpoints)

## build dataframe for covariates and response
datBMI <- nhanes3.1[,c("seqn","cpm","age","gender","race",
                       "education","BMI","sleep")]
nhanes1sub <- nhanes1[nhanes1$seqn %in% nhanes3.1$seqn,]
datBMI[,valid_days:=nhanes1sub$valid_days] %>% .[,sleepcounts:=sleep*valid_days]

## Add sleep to the frequency matrix
histmt[,1] <- histmt[,1] + datBMI$sleepcounts

############### Correlation Plot ###############
histcorr <- rcorr(histmt[,-101])
histcorr_1 <- histcorr$r
colnames(histcorr_1) <- c("0","","","","","500","","","","",
                          "1000","","","","","1500","","","","",
                          "2000","","","","","2500","","","","",
                          "3000","","","","","3500","","","","",
                          "4000","","","","","4500","","","","",
                          "5000","","","","","5500","","","","",
                          "6000","","","","","6500","","","","",
                          "7000","","","","","7500","","","","",
                          "8000","","","","","8500","","","","",
                          "9000","","","","","9500","","","","")
rownames(histcorr_1) <- c("0","","","","","500","","","","",
                          "1000","","","","","1500","","","","",
                          "2000","","","","","2500","","","","",
                          "3000","","","","","3500","","","","",
                          "4000","","","","","4500","","","","",
                          "5000","","","","","5500","","","","",
                          "6000","","","","","6500","","","","",
                          "7000","","","","","7500","","","","",
                          "8000","","","","","8500","","","","",
                          "9000","","","","","9500","","","","")
corrplot(histcorr_1, type = "lower", order = "original", 
         method = "square",
         tl.cex = 1.2, 
         tl.col = "black", tl.srt = 45,
         col=colorRampPalette(c("darkred", "white", "steelblue"))(200))

############### GAM Model Fitting ###############

### GAM, using default basis function "thin plate"
## normal with no link 
model1 <- gam(BMI ~ age + gender + race + education + 
                s(midpointsmt[,-1], k=10, by=histmt[,-1]), 
              data=datBMI, method="REML") 
summary.gam(model1)
AIC(model1)

# diagnosis plots
par(mfrow=c(2,2))
gam.check(model1)
check(modelviz.1,
      a.qq = list(method = "tnorm", 
                  a.cipoly = list(fill = "light blue")), 
      a.respoi = list(size = 0.5), 
      a.hist = list(bins = 10))

# plot smooth curve
par(mfrow=c(1,1))
modelviz.1 <- getViz(model1)
viz.1 <- plot(sm(modelviz.1, 1))
viz.1 +
  xlab("Intensity of PA (cpm)")  +  
  ylab("Estimates change in BMI with 95% CI") +
  scale_x_continuous( expand = c(0,0) , limits = c(0,10300) ) +
  geom_rect(data = rects, aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
            fill =  c("green", "DodgerBlue", "purple"),
            # Control the shading opacity here
            inherit.aes = FALSE, alpha = 0.1) +
  geom_vline(xintercept = 100, linetype="dotted", 
             color = "grey", size=0.8) +
  geom_vline(xintercept = 1951, linetype="dotted", 
             color = "grey", size=0.8) +
  ggtitle("Estimated Smooth Function for PA with 
          Normal Distribution and Identity Link Function ") +
  theme_light() + 
  theme(plot.title = element_text(hjust = 0.5))

## gaussian with log link
model2 <- gam(BMI ~ age + gender + race + education + 
                s(midpointsmt[,-1], k=10, by=histmt[,-1]),
              family = gaussian(link = log),data=datBMI, method="REML") 
par(mfrow=c(2,2))

# diagnosis plots
gam.check(model2)

# plot smooth curve
par(mfrow=c(1,1))
plot(model2)

## gamma with log link
model3 <- gam(BMI ~ age + gender + race + education + 
                s(midpointsmt[,-1], k=10, by=histmt[,-1]),
              family = Gamma(link = log),data=datBMI, method="REML") 
summary.gam(model3)
AIC(model3)

# diagnosis plots
par(mfrow=c(2,2))
gam.check(model3)
# plot smooth curve
modelviz.3 <- getViz(model3)
check(modelviz.3,
      a.qq = list(method = "tnorm", 
                  a.cipoly = list(fill = "light blue")), 
      a.respoi = list(size = 0.5), 
      a.hist = list(bins = 10))

# plot smooth curve
par(mfrow=c(1,1))
plot(model3)
viz.3 <- plot(sm(modelviz.3, 1))
viz.3 +
  xlab("Intensity of PA (cpm)")  +  
  ylab("Estimates change in BMI with 95% CI") +
  scale_x_continuous( expand = c(0,0) , limits = c(0,10300) ) +
  geom_rect(data = rects, aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
            fill =  c("green", "DodgerBlue", "purple"),
            # Control the shading opacity here
            inherit.aes = FALSE, alpha = 0.1) +
  geom_vline(xintercept = 100, linetype="dotted", 
             color = "grey", size=0.8) +
  geom_vline(xintercept = 1951, linetype="dotted", 
             color = "grey", size=0.8) +
  ggtitle("Estimated Smooth Function for PA 
          with Gamma Distribution and Log Link Function ") +
  theme_light() + 
  theme(plot.title = element_text(hjust = 0.5))

## removing linear part of functional predictor and add it as mean cpm for model 3
model3.cpm <- gam(BMI ~ age + cpm + gender + race + education + 
                    s(midpointsmt, k=10, by=histmt, bs="tp", m=c(2,0)),
                  family = Gamma(link = log),data=datBMI, method="REML") 
summary.gam(model3.cpm)
AIC(model3.cpm)

# diagnosis plots
par(mfrow=c(2,2))
gam.check(model3.cpm)
modelviz.3cpm <- getViz(model3.cpm)
check(modelviz.3cpm,
      a.qq = list(method = "tnorm", 
                  a.cipoly = list(fill = "light blue")), 
      a.respoi = list(size = 0.5), 
      a.hist = list(bins = 10))

# plot smooth curve
par(mfrow=c(1,1))
plot.gam(model3.cpm) 

viz.3cpm <- plot(sm(modelviz.3cpm, 1))
viz.3cpm +
  xlab("Intensity of PA (cpm)")  +  
  ylab("Estimates change in BMI with 95% CI") +
  scale_x_continuous( expand = c(0,0) , limits = c(0,10300) ) +
  geom_rect(data = rects, aes(xmin = xmin, xmax = xmax, 
                              ymin = ymin, ymax = ymax),
            fill =  c("green", "DodgerBlue", "purple"),
            # Control the shading opacity here
            inherit.aes = FALSE, alpha = 0.1) +
  geom_vline(xintercept = 100, linetype="dotted", 
             color = "grey", size=0.8) +
  geom_vline(xintercept = 1951, linetype="dotted", 
             color = "grey", size=0.8) +
  ggtitle("Estimated Smooth Function for PA (CPM removed) with 
          Gamma Distribution and Log Link Function ") +
  theme_light() + 
  theme(plot.title = element_text(hjust = 0.5))
