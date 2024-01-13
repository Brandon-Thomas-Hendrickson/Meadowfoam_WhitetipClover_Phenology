###################
#Download Packages#
###################
install.packages(c("nortest","readr","sjPlot","ggplot2","lmtest","AICcmodavg","DescTools","rcompanion"))
lapply(c("nortest","readr","sjPlot","ggplot2","lmtest","AICcmodavg","DescTools","rcompanion"), require, character.only = TRUE)
list.of.packages <- c("nortest","readr","sjPlot","ggplot2","lmtest","AICcmodavg","DescTools","rcompanion")
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)
##############
#Upload Files#
##############

merced_climate <- read.csv("~/Documents/Github_Projects/Meadowfoam_WhitetipClover_Phenology/Data/csv_files/MERCED_CLIMATE_2000_2022.csv")
climate_stations <- read.csv("~/Documents/Github_Projects/Meadowfoam_WhitetipClover_Phenology/Data/csv_files/CIMIS_SITES.csv")
Phenomerged <- read.csv("~/Documents/Github_Projects/Meadowfoam_WhitetipClover_Phenology/Data/csv_files/Phenomerged.csv")
ClayFlower <- read.csv("~/Documents/Github_Projects/Meadowfoam_WhitetipClover_Phenology/Data/csv_files/ClayFlower.csv")
merced_daily_climate <- read.csv("~/Documents/Github_Projects/Meadowfoam_WhitetipClover_Phenology/Data/csv_files/merced_daily_climate.csv")

#################################################
#Calculate Growing Degree Days and Precipitation#
#################################################

#Create Growing Degree Hours function
growing_season <- na.omit(merced_climate[(merced_climate$MONTH>=10 | merced_climate$MONTH<4),])
GDH<-function(x,t,GDH){
    if(any(class(x)=="data.frame")==TRUE){
        output <- character (nrow(x))
        for (i in 1:nrow(x)){
            if(x[i,"TEMP"] > t){
                output[i]<-1
            } else {
                output[i]<-0
            }
        }
    } else {
        print("Not a data.frame")
    }
    x[[GDH]]<-output
    return(x)
}

#Acquire Accumulated GDH and Precipitation in Merced For all Years (2000-2022)
growing_season<-GDH(growing_season,9.44,GDH="GDH_MERCED")
EARLY_climate <- subset(growing_season,growing_season$MONTH>=10)
EARLY_gdh<-split(as.numeric(EARLY_climate$GDH_MERCED),EARLY_climate$YEAR)
EARLY_precip <- split(as.numeric(EARLY_climate$PRECIP),EARLY_climate$YEAR)
LATE_climate <- subset(growing_season,growing_season$MONTH<=3)
LATE_gdh<-split(as.numeric(LATE_climate$GDH_MERCED),LATE_climate$YEAR)
LATE_precip <- split(as.numeric(LATE_climate$PRECIP),LATE_climate$YEAR)
tab_sum <- data.frame("Accumulated_GDH"=NULL,"Accumulated_PRECIP"=NULL,"SEASON"=NULL)
for(i in 1:length(EARLY_gdh)){
    sum_gdh <- sum(EARLY_gdh[[i]])
    sum_precip <- sum(EARLY_precip[[i]])
    tab_sum[nrow(tab_sum)+1,"Accumulated_GDH"]<-sum_gdh
    tab_sum[nrow(tab_sum),"Accumulated_PRECIP"]<-sum_precip
    tab_sum[nrow(tab_sum),"SEASON"]<-"EARLY"
    tab_sum[nrow(tab_sum),"YEAR"]<-as.numeric(names(EARLY_gdh[i]))
    sum_gdh <- sum(LATE_gdh[[i]])
    sum_precip <- sum(LATE_precip[[i]])
    tab_sum[nrow(tab_sum)+1,"Accumulated_GDH"]<-sum_gdh
    tab_sum[nrow(tab_sum),"Accumulated_PRECIP"]<-sum_precip
    tab_sum[nrow(tab_sum),"SEASON"]<-"LATE"
    tab_sum[nrow(tab_sum),"YEAR"]<-as.numeric(names(EARLY_gdh[i]))
}

#Calculate Diural Temperature Range
merced_daily_climate$TempRange<-merced_daily_climate$Max.Air.Temp..C. - merced_daily_climate$Min.Air.Temp..C.

#Calculate CMI 
merced_daily_climate$CMI <- (merced_daily_climate$Precip..mm. - merced_daily_climate$PM.ETo..mm.) + 100

###########################################
#Analyze Fidelity of CIMIS Climate Station#
###########################################
#Aquire total precipation for all stations each year
climate_stations<-climate_stations[!is.na(climate_stations$PRECIP),]
EARLY_climate <- subset(climate_stations,climate_stations$MONTH>=10)
EARLY_precip <- split(as.numeric(EARLY_climate$PRECIP), list(EARLY_climate$YEAR,EARLY_climate$GROUP))
LATE_climate <- subset(climate_stations,climate_stations$MONTH<=3)
LATE_precip <- split(as.numeric(LATE_climate$PRECIP), list(LATE_climate$YEAR,LATE_climate$GROUP))
tab_sum <- data.frame("Accumulated_PRECIP"=0,"SEASON"=0,"GROUP"=0,"YEAR"=0,"SITE"=0)
for(i in 1:length(EARLY_precip)){
    sum_precip <- sum(EARLY_precip[[i]])
    tab_sum[nrow(tab_sum)+1,"Accumulated_PRECIP"]<-sum_precip
    tab_sum[nrow(tab_sum),"SEASON"]<-"EARLY"
    tab_sum[nrow(tab_sum),"GROUP"]<-names(EARLY_precip[i])
    sum_precip <- sum(LATE_precip[[i]])
    tab_sum[nrow(tab_sum)+1,"Accumulated_PRECIP"]<-sum_precip
    tab_sum[nrow(tab_sum),"SEASON"]<-"LATE"
    tab_sum[nrow(tab_sum),"GROUP"]<-names(EARLY_precip[i])
}
for(i in 1:nrow(tab_sum)){
    spli <- strsplit(tab_sum[i,"GROUP"],split="[.]")
    tab_sum[i,"YEAR"]<-spli[[1]][1]
    tab_sum[i,"SITE"]<-spli[[1]][2]
}
tab_sum<-tab_sum[-1,]
tab_sum$SITE<-as.factor(tab_sum$SITE)
tab_early<-tab_sum[tab_sum$SEASON=="EARLY",]
tab_late<-tab_sum[tab_sum$SEASON=="LATE",]
model<-aov(Accumulated_PRECIP~SITE,data=tab_early)
DunnettTest(x=tab_early$Accumulated_PRECIP,g=tab_early$SITE)
model<-aov(Accumulated_PRECIP~SITE,data=tab_late)
DunnettTest(x=tab_late$Accumulated_PRECIP,g=tab_late$SITE)

################################
#Summarize Climate Variables####
################################
climate_vars<-list("GDDEARLY","GDDLATE","MEANEARLY","MEANLATE","MinEARLY","MinLate","MAXEarly","MAXLate","PRECIP_EARLY","PRECIP_LATE","CMI_EARLY","CMI_LATE","TempRange_EARLY","TempRange_LATE")
climate_summary <- list()
for(i in 1:length(climate_vars)){
    name = climate_vars[[i]][1]
    climate_summary[[name]]<-summary(Phenomerged[,name])
}

#Run Linear Regressions of Climate Variables and Year
pd<-Phenomerged[Phenomerged$Pool==1,]
lm_list_climate<-list()
for(i in 1:length(climate_vars)){
    name = climate_vars[[i]]
    lm_list_climate[[name]]<-summary(lm(as.matrix(pd[,name])~Year,data=pd))
}
#Run Linear Regression of GDH with Diural Temperature Range
summary(lm(GDDEARLY~TempRange_EARLY+MAXLate,data=pd))
summary(lm(GDDLATE~TempRange_LATE+MAXLate,data=pd))
#Run Linear Regression of Diural Temperature Range with CMI
summary(lm(TempRange_EARLY~CMI_EARLY,data=pd))
summary(lm(TempRange_LATE~CMI_LATE,data=pd))

################################################################
#Determine Normality of Phenological Variables#
################################################################
var_pheno <- as.matrix(Phenomerged[,3:14],header=T)
AD_list <- list()
var_pheno1<-as.data.frame(var_pheno)
for (i in colnames(var_pheno)){
    name = paste(i,"_normality_test",sep="")
    ad_test <- ad.test(var_pheno[,i])
    if(ad_test$p.value < 0.05){
        BOXCOX<-transformTukey(var_pheno[,i])
        var_pheno1[,ncol(var_pheno1)+1]<-BOXCOX
        colnames(var_pheno1)[ncol(var_pheno1)]<- paste("BOX_COX",i,sep=".")
    }
    AD_list[[name]]<-ad_test
}
#Of the variables that are non-normally distributed, do transformed variables produce a lower AIC score for polynomial regressions?
log_group<-grepl("BOX_COX.",colnames(var_pheno1))
log_data <- var_pheno1[log_group]
len <- length(log_data)
lm_list <- list()
lm_list_Untransformed <- list()
x <- sub("BOX_COX.","",colnames(var_pheno1[log_group]))
for (i in seq(1,len)){
    for (q in seq(1,4)){
    name <- paste(colnames(log_data[i]),q,sep="_")
    lm <- lm(log_data[,i]~poly(Phenomerged$Year,q,raw=TRUE))
    lm_list[[name]] <- lm
    }
}
AIC_Scores <- data.frame(BOX_COX=NULL)
for (i in lm_list){
    aic <- AIC(i)
    AIC_Scores[nrow(AIC_Scores)+1,"BOX_COX"] <- aic
}
for (i in seq(1:len)){
    for (q in seq(1,4)){
        name <- paste(x[i],q,sep=".Poly_")
        lm <- lm(as.matrix(var_pheno1[x[i]])~poly(Phenomerged$Year,q,raw=TRUE))
        lm_list_Untransformed[[name]]<-lm
    }
}
AIC_Scores_UnTransformed<-data.frame(UnTransformed=NULL)
for (i in lm_list_Untransformed){
    aic <- AIC(i)
    AIC_Scores_UnTransformed[nrow(AIC_Scores_UnTransformed)+1,"UnTransformed"]<-aic
}
AIC_Scores <- cbind("UnTransformed"=AIC_Scores_UnTransformed$UnTransformed,AIC_Scores)

for(i in 1:nrow(AIC_Scores)){
    if(AIC_Scores$UnTransformed[i]<AIC_Scores$BOX_COX[i]){
        AIC_Scores$Winner <- print("UnTransformed")
    } else {
        AIC_Scores$Winner <- print("Transformed")
    }
}
row.names(AIC_Scores)<-names(lm_list_Untransformed)

##################################################################
#Run Multiple Regression Models: Phenology and Population Density#
##################################################################
lm_list <- list()
BP_table_PHENO<-data.frame("Heteroscedastic"=NULL)
Breusch_Pagan_list_Pheno<-list()
for(W in names(Phenomerged[,3:14])){
    Climate_List<-list(Phenomerged$GDDEARLY,Phenomerged$GDDLATE,Phenomerged$PRECIP_EARLY,Phenomerged$PRECIP_LATE)
    names(Climate_List)<-c("GDDEARLY","GDDLATE","PRECIP_EARLY","PRECIP_LATE")
    for(pp in seq(1,5)){
        a <- length(Climate_List)
        b <- a - 1
        c <- b - 1
        d <- c - 1
        ax <- names(Climate_List[a])
        bx <- names(Climate_List[b])
        cx <- names(Climate_List[c])
        dx <- names(Climate_List[d])
        if(a>=1 && b>=1 && c>=1 && d>= 1 ){
            lm <- lm(as.matrix(Phenomerged[,W]) ~ Climate_List[[a]]+Climate_List[[b]]+Climate_List[[c]]+Climate_List[[d]],data=Phenomerged)
            names(lm$coefficients)<-c("Intercept",ax,bx,cx,dx)
            fail <- "NO"
        }else if(a>=1 && b>=1 && c>=1){
            lm <- lm(as.matrix(Phenomerged[,W]) ~ Climate_List[[a]]+Climate_List[[b]]+Climate_List[[c]])
            names(lm$coefficients)<-c("Intercept",ax,bx,cx)
            fail <- "NO"
        }else if(a>=1 && b>=1){
            lm <- lm(as.matrix(Phenomerged[,W]) ~ Climate_List[[a]]+Climate_List[[b]])
            names(lm$coefficients)<-c("Intercept",ax,bx)
            fail <- "NO"
        }else if(a>=1){
            lm <- lm(as.matrix(Phenomerged[,W]) ~ Climate_List[[a]])
            names(lm$coefficients)<-c("Intercept",ax)
            fail <- "NO"
        }else {
            fail <- "YES"
        }
        temp <- data.frame(Variable=NULL,Value=NULL)
        for(i in 1:(length(Climate_List))+1){
            if(summary(lm)$coefficients[i,4]>0.05){
            Var<-(names(summary(lm)$coefficients[,4][i]))
            Val<-summary(lm)$coefficients[Var,4]
            temp[nrow(temp)+1,"Variable"]<-Var
            temp[nrow(temp),"Value"]<-Val
            } 
    }
        if(any(temp$Value>0.05)==TRUE){
        x <- temp[which.max(temp$Value),]
        x <- x[,1]
        Climate_List<-Climate_List[names(Climate_List) != x]
        }   
    }
    if(fail == "NO"){
    lm_list[[W]] <- lm
    bpmodel <- bptest(lm)
    Breusch_Pagan_list_Pheno[[W]] <- bpmodel
    }
    if(bpmodel$p.value<0.05){
    BP_table_PHENO[nrow(BP_table_PHENO)+1,"Heteroscedastic"]<-print("True")
        }else   {
    BP_table_PHENO[nrow(BP_table_PHENO)+1,"Heteroscedastic"]<-print("False")
    }
}
#CMI 
lm_list_CMI <- list()
BP_table_PHENO<-data.frame("Heteroscedastic"=NULL)
Breusch_Pagan_list_Pheno<-list()
for(W in names(Phenomerged[,3:14])){
    Climate_List<-list(Phenomerged$GDDEARLY,Phenomerged$GDDLATE,Phenomerged$CMI_EARLY,Phenomerged$CMI_LATE)
    names(Climate_List)<-c("GDDEARLY","GDDLATE","CMI_EARLY","CMI_LATE")
    for(pp in seq(1,5)){
        a <- length(Climate_List)
        b <- a - 1
        c <- b - 1
        d <- c - 1
        ax <- names(Climate_List[a])
        bx <- names(Climate_List[b])
        cx <- names(Climate_List[c])
        dx <- names(Climate_List[d])
        if(a>=1 && b>=1 && c>=1 && d>= 1 ){
            lm <- lm(as.matrix(Phenomerged[,W]) ~ Climate_List[[a]]+Climate_List[[b]]+Climate_List[[c]]+Climate_List[[d]],data=Phenomerged)
            names(lm$coefficients)<-c("Intercept",ax,bx,cx,dx)
            fail <- "NO"
        }else if(a>=1 && b>=1 && c>=1){
            lm <- lm(as.matrix(Phenomerged[,W]) ~ Climate_List[[a]]+Climate_List[[b]]+Climate_List[[c]])
            names(lm$coefficients)<-c("Intercept",ax,bx,cx)
            fail <- "NO"
        }else if(a>=1 && b>=1){
            lm <- lm(as.matrix(Phenomerged[,W]) ~ Climate_List[[a]]+Climate_List[[b]])
            names(lm$coefficients)<-c("Intercept",ax,bx)
            fail <- "NO"
        }else if(a>=1){
            lm <- lm(as.matrix(Phenomerged[,W]) ~ Climate_List[[a]])
            names(lm$coefficients)<-c("Intercept",ax)
            fail <- "NO"
        }else {
            fail <- "YES"
        }
        temp <- data.frame(Variable=NULL,Value=NULL)
        for(i in 1:(length(Climate_List))+1){
            if(summary(lm)$coefficients[i,4]>0.05){
            Var<-(names(summary(lm)$coefficients[,4][i]))
            Val<-summary(lm)$coefficients[Var,4]
            temp[nrow(temp)+1,"Variable"]<-Var
            temp[nrow(temp),"Value"]<-Val
            } 
    }
        if(any(temp$Value>0.05)==TRUE){
        x <- temp[which.max(temp$Value),]
        x <- x[,1]
        Climate_List<-Climate_List[names(Climate_List) != x]
        }   
    }
    if(fail == "NO"){
    lm_list_CMI[[W]] <- lm
    bpmodel <- bptest(lm)
    Breusch_Pagan_list_Pheno[[W]] <- bpmodel
    }
    if(bpmodel$p.value<0.05){
    BP_table_PHENO[nrow(BP_table_PHENO)+1,"Heteroscedastic"]<-print("True")
        }else   {
    BP_table_PHENO[nrow(BP_table_PHENO)+1,"Heteroscedastic"]<-print("False")
    }
}
#Mean
lm_list_MEAN <- list()
BP_table_PHENO<-data.frame("Heteroscedastic"=NULL)
Breusch_Pagan_list_Pheno<-list()
for(W in names(Phenomerged[,3:14])){
    Climate_List<-list(Phenomerged$MEANEARLY,Phenomerged$MEANLATE,Phenomerged$PRECIP_EARLY,Phenomerged$PRECIP_LATE)
    names(Climate_List)<-c("MEANEarly","MEANLate","PRECIP_Early","PRECIP_Late")
    for(pp in seq(1,5)){
        a <- length(Climate_List)
        b <- a - 1
        c <- b - 1
        d <- c - 1
        ax <- names(Climate_List[a])
        bx <- names(Climate_List[b])
        cx <- names(Climate_List[c])
        dx <- names(Climate_List[d])
        if(a>=1 && b>=1 && c>=1 && d>= 1 ){
            lm <- lm(as.matrix(Phenomerged[,W]) ~ Climate_List[[a]]+Climate_List[[b]]+Climate_List[[c]]+Climate_List[[d]],data=Phenomerged)
            names(lm$coefficients)<-c("Intercept",ax,bx,cx,dx)
            fail <- "NO"
        }else if(a>=1 && b>=1 && c>=1){
            lm <- lm(as.matrix(Phenomerged[,W]) ~ Climate_List[[a]]+Climate_List[[b]]+Climate_List[[c]])
            names(lm$coefficients)<-c("Intercept",ax,bx,cx)
            fail <- "NO"
        }else if(a>=1 && b>=1){
            lm <- lm(as.matrix(Phenomerged[,W]) ~ Climate_List[[a]]+Climate_List[[b]])
            names(lm$coefficients)<-c("Intercept",ax,bx)
            fail <- "NO"
        }else if(a>=1){
            lm <- lm(as.matrix(Phenomerged[,W]) ~ Climate_List[[a]])
            names(lm$coefficients)<-c("Intercept",ax)
            fail <- "NO"
        }else {
            fail <- "YES"
        }
        temp <- data.frame(Variable=NULL,Value=NULL)
        for(i in 1:(length(Climate_List))+1){
            if(summary(lm)$coefficients[i,4]>0.05){
            Var<-(names(summary(lm)$coefficients[,4][i]))
            Val<-summary(lm)$coefficients[Var,4]
            temp[nrow(temp)+1,"Variable"]<-Var
            temp[nrow(temp),"Value"]<-Val
            } 
    }
        if(any(temp$Value>0.05)==TRUE){
        x <- temp[which.MEAN(temp$Value),]
        x <- x[,1]
        Climate_List<-Climate_List[names(Climate_List) != x]
        }   
    }
    if(fail == "NO"){
    lm_list_MEAN[[W]] <- lm
    bpmodel <- bptest(lm)
    Breusch_Pagan_list_Pheno[[W]] <- bpmodel
    }
    if(bpmodel$p.value<0.05){
    BP_table_PHENO[nrow(BP_table_PHENO)+1,"Heteroscedastic"]<-print("True")
        }else   {
    BP_table_PHENO[nrow(BP_table_PHENO)+1,"Heteroscedastic"]<-print("False")
    }
}
#MAX temp
lm_list_MAX <- list()
BP_table_PHENO<-data.frame("Heteroscedastic"=NULL)
Breusch_Pagan_list_Pheno<-list()
for(W in names(Phenomerged[,3:14])){
    Climate_List<-list(Phenomerged$MAXEarly,Phenomerged$MAXLate,Phenomerged$PRECIP_EARLY,Phenomerged$PRECIP_LATE)
    names(Climate_List)<-c("MAXEarly","MAXLate","PRECIP_EARLY","PRECIP_LATE")
    for(pp in seq(1,5)){
        a <- length(Climate_List)
        b <- a - 1
        c <- b - 1
        d <- c - 1
        ax <- names(Climate_List[a])
        bx <- names(Climate_List[b])
        cx <- names(Climate_List[c])
        dx <- names(Climate_List[d])
        if(a>=1 && b>=1 && c>=1 && d>= 1 ){
            lm <- lm(as.matrix(Phenomerged[,W]) ~ Climate_List[[a]]+Climate_List[[b]]+Climate_List[[c]]+Climate_List[[d]],data=Phenomerged)
            names(lm$coefficients)<-c("Intercept",ax,bx,cx,dx)
            fail <- "NO"
        }else if(a>=1 && b>=1 && c>=1){
            lm <- lm(as.matrix(Phenomerged[,W]) ~ Climate_List[[a]]+Climate_List[[b]]+Climate_List[[c]])
            names(lm$coefficients)<-c("Intercept",ax,bx,cx)
            fail <- "NO"
        }else if(a>=1 && b>=1){
            lm <- lm(as.matrix(Phenomerged[,W]) ~ Climate_List[[a]]+Climate_List[[b]])
            names(lm$coefficients)<-c("Intercept",ax,bx)
            fail <- "NO"
        }else if(a>=1){
            lm <- lm(as.matrix(Phenomerged[,W]) ~ Climate_List[[a]])
            names(lm$coefficients)<-c("Intercept",ax)
            fail <- "NO"
        }else {
            fail <- "YES"
        }
        temp <- data.frame(Variable=NULL,Value=NULL)
        for(i in 1:(length(Climate_List))+1){
            if(summary(lm)$coefficients[i,4]>0.05){
            Var<-(names(summary(lm)$coefficients[,4][i]))
            Val<-summary(lm)$coefficients[Var,4]
            temp[nrow(temp)+1,"Variable"]<-Var
            temp[nrow(temp),"Value"]<-Val
            } 
    }
        if(any(temp$Value>0.05)==TRUE){
        x <- temp[which.max(temp$Value),]
        x <- x[,1]
        Climate_List<-Climate_List[names(Climate_List) != x]
        }   
    }
    if(fail == "NO"){
    lm_list_MAX[[W]] <- lm
    bpmodel <- bptest(lm)
    Breusch_Pagan_list_Pheno[[W]] <- bpmodel
    }
    if(bpmodel$p.value<0.05){
    BP_table_PHENO[nrow(BP_table_PHENO)+1,"Heteroscedastic"]<-print("True")
        }else   {
    BP_table_PHENO[nrow(BP_table_PHENO)+1,"Heteroscedastic"]<-print("False")
    }
}
#MIN
lm_list_Min <- list()
BP_table_PHENO<-data.frame("Heteroscedastic"=NULL)
Breusch_Pagan_list_Pheno<-list()
for(W in names(Phenomerged[,3:14])){
    Climate_List<-list(Phenomerged$MinEARLY,Phenomerged$MinLate,Phenomerged$PRECIP_EARLY,Phenomerged$PRECIP_LATE)
    names(Climate_List)<-c("MinEarly","MinLate","PRECIP_EARLY","PRECIP_LATE")
    for(pp in seq(1,5)){
        a <- length(Climate_List)
        b <- a - 1
        c <- b - 1
        d <- c - 1
        ax <- names(Climate_List[a])
        bx <- names(Climate_List[b])
        cx <- names(Climate_List[c])
        dx <- names(Climate_List[d])
        if(a>=1 && b>=1 && c>=1 && d>= 1 ){
            lm <- lm(as.matrix(Phenomerged[,W]) ~ Climate_List[[a]]+Climate_List[[b]]+Climate_List[[c]]+Climate_List[[d]],data=Phenomerged)
            names(lm$coefficients)<-c("Intercept",ax,bx,cx,dx)
            fail <- "NO"
        }else if(a>=1 && b>=1 && c>=1){
            lm <- lm(as.matrix(Phenomerged[,W]) ~ Climate_List[[a]]+Climate_List[[b]]+Climate_List[[c]])
            names(lm$coefficients)<-c("Intercept",ax,bx,cx)
            fail <- "NO"
        }else if(a>=1 && b>=1){
            lm <- lm(as.matrix(Phenomerged[,W]) ~ Climate_List[[a]]+Climate_List[[b]])
            names(lm$coefficients)<-c("Intercept",ax,bx)
            fail <- "NO"
        }else if(a>=1){
            lm <- lm(as.matrix(Phenomerged[,W]) ~ Climate_List[[a]])
            names(lm$coefficients)<-c("Intercept",ax)
            fail <- "NO"
        }else {
            fail <- "YES"
        }
        temp <- data.frame(Variable=NULL,Value=NULL)
        for(i in 1:(length(Climate_List))+1){
            if(summary(lm)$coefficients[i,4]>0.05){
            Var<-(names(summary(lm)$coefficients[,4][i]))
            Val<-summary(lm)$coefficients[Var,4]
            temp[nrow(temp)+1,"Variable"]<-Var
            temp[nrow(temp),"Value"]<-Val
            } 
    }
        if(any(temp$Value>0.05)==TRUE){
        x <- temp[which.Min(temp$Value),]
        x <- x[,1]
        Climate_List<-Climate_List[names(Climate_List) != x]
        }   
    }
    if(fail == "NO"){
    lm_list_Min[[W]] <- lm
    bpmodel <- bptest(lm)
    Breusch_Pagan_list_Pheno[[W]] <- bpmodel
    }
    if(bpmodel$p.value<0.05){
    BP_table_PHENO[nrow(BP_table_PHENO)+1,"Heteroscedastic"]<-print("True")
        }else   {
    BP_table_PHENO[nrow(BP_table_PHENO)+1,"Heteroscedastic"]<-print("False")
    }
}
###############################################################
#Run Linear Regression Models: Phenology by Distance from Edge#
###############################################################
distance_lm_list=list()
for(i in names(ClayFlower[,5:10])){
    name = i
    distance_lm_list[[name]]<-lm(as.matrix(ClayFlower[,i])~DistanceFromEdge,data=ClayFlower)
}
clay_lm_list=list()
for(i in names(ClayFlower[,5:10])){
    name = i
    clay_lm_list[[name]]<-lm(as.matrix(ClayFlower[,i])~Clay_Percent,data=ClayFlower)
}

################################
#Analyze Rate of Pool Occupancy#
################################
#Meadowfoam 
lm_list <- list()
BP_table_PHENO<-data.frame("Heteroscedastic"=NULL)
Breusch_Pagan_list_Pheno<-list()
for(W in names(Phenomerged[,15:25])){
    Climate_List<-list(Phenomerged$GDD_WINTER,Phenomerged$PRECIP_WINTER,Phenomerged$MFStart,Phenomerged$MFEnd)
    names(Climate_List)<-c("GDD_WINTER","PRECIP_WINTER","MFSTART","MFEND")
    for(pp in seq(1,5)){
        a <- length(Climate_List)
        b <- a - 1
        c <- b - 1
        d <- c - 1
        ax <- names(Climate_List[a])
        bx <- names(Climate_List[b])
        cx <- names(Climate_List[c])
        dx <- names(Climate_List[d])
        if(a>=1 && b>=1 && c>=1 && d>= 1 ){
            lm <- lm(as.matrix(Phenomerged[,W]) ~ Climate_List[[a]]+Climate_List[[b]]+Climate_List[[c]]+Climate_List[[d]],data=Phenomerged)
            names(lm$coefficients)<-c("Intercept",ax,bx,cx,dx)
            fail <- "NO"
        }else if(a>=1 && b>=1 && c>=1){
            lm <- lm(as.matrix(Phenomerged[,W]) ~ Climate_List[[a]]+Climate_List[[b]]+Climate_List[[c]])
            names(lm$coefficients)<-c("Intercept",ax,bx,cx)
            fail <- "NO"
        }else if(a>=1 && b>=1){
            lm <- lm(as.matrix(Phenomerged[,W]) ~ Climate_List[[a]]+Climate_List[[b]])
            names(lm$coefficients)<-c("Intercept",ax,bx)
            fail <- "NO"
        }else if(a>=1){
            lm <- lm(as.matrix(Phenomerged[,W]) ~ Climate_List[[a]])
            names(lm$coefficients)<-c("Intercept",ax)
            fail <- "NO"
        }else {
            fail <- "YES"
        }
        temp <- data.frame(Variable=NULL,Value=NULL)
        for(i in 1:(length(Climate_List))+1){
            if(summary(lm)$coefficients[i,4]>0.05){
                Var<-(names(summary(lm)$coefficients[,4][i]))
                Val<-summary(lm)$coefficients[Var,4]
                temp[nrow(temp)+1,"Variable"]<-Var
                temp[nrow(temp),"Value"]<-Val
            } 
        }
        if(any(temp$Value>0.05)==TRUE){
            x <- temp[which.max(temp$Value),]
            x <- x[,1]
            Climate_List<-Climate_List[names(Climate_List) != x]
        }   
    }
    if(fail == "NO"){
        lm_list[[W]] <- lm
        bpmodel <- bptest(lm)
        Breusch_Pagan_list_Pheno[[W]] <- bpmodel
    }
    if(bpmodel$p.value<0.05){
        BP_table_PHENO[nrow(BP_table_PHENO)+1,"Heteroscedastic"]<-print("True")
    }else   {
        BP_table_PHENO[nrow(BP_table_PHENO)+1,"Heteroscedastic"]<-print("False")
    }
}
#Whitetip Clover
for(W in names(Phenomerged[,26:36])){
    Climate_List<-list(Phenomerged$GDD_WINTER,Phenomerged$PRECIP_WINTER,Phenomerged$WTCStart,Phenomerged$WTCEnd)
    names(Climate_List)<-c("GDD_WINTER","PRECIP_WINTER","WTSTART","WTEND")
    for(pp in seq(1,5)){
        a <- length(Climate_List)
        b <- a - 1
        c <- b - 1
        d <- c - 1
        ax <- names(Climate_List[a])
        bx <- names(Climate_List[b])
        cx <- names(Climate_List[c])
        dx <- names(Climate_List[d])
        if(a>=1 && b>=1 && c>=1 && d>= 1 ){
            lm <- lm(as.matrix(Phenomerged[,W]) ~ Climate_List[[a]]+Climate_List[[b]]+Climate_List[[c]]+Climate_List[[d]],data=Phenomerged)
            names(lm$coefficients)<-c("Intercept",ax,bx,cx,dx)
            fail <- "NO"
        }else if(a>=1 && b>=1 && c>=1){
            lm <- lm(as.matrix(Phenomerged[,W]) ~ Climate_List[[a]]+Climate_List[[b]]+Climate_List[[c]])
            names(lm$coefficients)<-c("Intercept",ax,bx,cx)
            fail <- "NO"
        }else if(a>=1 && b>=1){
            lm <- lm(as.matrix(Phenomerged[,W]) ~ Climate_List[[a]]+Climate_List[[b]])
            names(lm$coefficients)<-c("Intercept",ax,bx)
            fail <- "NO"
        }else if(a>=1){
            lm <- lm(as.matrix(Phenomerged[,W]) ~ Climate_List[[a]])
            names(lm$coefficients)<-c("Intercept",ax)
            fail <- "NO"
        }else {
            fail <- "YES"
        }
        temp <- data.frame(Variable=NULL,Value=NULL)
        for(i in 1:(length(Climate_List))+1){
            if(summary(lm)$coefficients[i,4]>0.05){
                Var<-(names(summary(lm)$coefficients[,4][i]))
                Val<-summary(lm)$coefficients[Var,4]
                temp[nrow(temp)+1,"Variable"]<-Var
                temp[nrow(temp),"Value"]<-Val
            } 
        }
        if(any(temp$Value>0.05)==TRUE){
            x <- temp[which.max(temp$Value),]
            x <- x[,1]
            Climate_List<-Climate_List[names(Climate_List) != x]
        }   
    }
    if(fail == "NO"){
        lm_list[[W]] <- lm
        bpmodel <- bptest(lm)
        Breusch_Pagan_list_Pheno[[W]] <- bpmodel
    }
    if(bpmodel$p.value<0.05){
        BP_table_PHENO[nrow(BP_table_PHENO)+1,"Heteroscedastic"]<-print("True")
    }else   {
        BP_table_PHENO[nrow(BP_table_PHENO)+1,"Heteroscedastic"]<-print("False")
    }
}


#########
#Figures#
#########

#Figure 1: Phenology Plots
av_list <- list()
for(i in 1:length(lm_list)){
    name = names(lm_list[i])
    av_list[[name]] <- avPlots(lm_list[[i]])
}
mfstart<-ggplot(data=as.data.frame(av_list[[1]][1]),aes(x=av_list[[1]][[1]][,1],y=av_list[[1]][[1]][,2]))+geom_point()+theme_light()+geom_smooth(method="lm", se=TRUE,color="black")+xlab("")+ylab("")+theme(axis.text.x = element_text(size=15),axis.text.y=element_text(size=15))
mfstart_b<-ggplot(data=as.data.frame(av_list[[1]][2]),aes(x=av_list[[1]][[2]][,1],y=av_list[[1]][[2]][,2]))+geom_point()+theme_light()+geom_smooth(method="lm", se=TRUE,color="black")+xlab("")+ylab("")+theme(axis.text.x = element_text(size=15),axis.text.y=element_text(size=15))
mfstart_c<-ggplot(data=as.data.frame(av_list[[1]][3]),aes(x=av_list[[1]][[3]][,1],y=av_list[[1]][[3]][,2]))+geom_point()+theme_light()+geom_smooth(method="lm", se=TRUE,color="black")+xlab("")+ylab("")+theme(axis.text.x = element_text(size=15),axis.text.y=element_text(size=15))

mfend<-ggplot(data=as.data.frame(av_list[[2]][1]),aes(x=av_list[[2]][[1]][,1],y=av_list[[2]][[1]][,2]))+geom_point()+theme_light()+geom_smooth(method="lm", se=TRUE,color="black")+xlab("")+ylab("")+theme(axis.text.x = element_text(size=15),axis.text.y=element_text(size=15))
mfend_b<-ggplot(data=as.data.frame(av_list[[2]][2]),aes(x=av_list[[2]][[2]][,1],y=av_list[[2]][[2]][,2]))+geom_point()+theme_light()+geom_smooth(method="lm", se=TRUE,color="black")+xlab("")+ylab("")+theme(axis.text.x = element_text(size=15),axis.text.y=element_text(size=15))

mfpeak<-ggplot(data=as.data.frame(av_list[[3]][1]),aes(x=av_list[[3]][[1]][,1],y=av_list[[3]][[1]][,2]))+geom_point()+theme_light()+geom_smooth(method="lm", se=TRUE,color="black")+xlab("")+ylab("")+theme(axis.text.x = element_text(size=15),axis.text.y=element_text(size=15))
mfpeak_b<-ggplot(data=as.data.frame(av_list[[3]][2]),aes(x=av_list[[3]][[2]][,1],y=av_list[[3]][[2]][,2]))+geom_point()+theme_light()+geom_smooth(method="lm", se=TRUE,color="black")+xlab("")+ylab("")+theme(axis.text.x = element_text(size=15),axis.text.y=element_text(size=15))

wtstart<-ggplot(data=as.data.frame(av_list[[4]][1]),aes(x=av_list[[4]][[1]][,1],y=av_list[[4]][[1]][,2]))+geom_point()+theme_light()+geom_smooth(method="lm", se=TRUE,color="black")+xlab("")+ylab("")+theme(axis.text.x = element_text(size=15),axis.text.y=element_text(size=15))
wtstart_b<-ggplot(data=as.data.frame(av_list[[4]][2]),aes(x=av_list[[4]][[2]][,1],y=av_list[[4]][[2]][,2]))+geom_point()+theme_light()+geom_smooth(method="lm", se=TRUE,color="black")+xlab("")+ylab("")+theme(axis.text.x = element_text(size=15),axis.text.y=element_text(size=15))
wtstart_c<-ggplot(data=as.data.frame(av_list[[4]][3]),aes(x=av_list[[4]][[3]][,1],y=av_list[[4]][[3]][,2]))+geom_point()+theme_light()+geom_smooth(method="lm", se=TRUE,color="black")+xlab("")+ylab("")+theme(axis.text.x = element_text(size=15),axis.text.y=element_text(size=15))

wtend<-ggplot(data=as.data.frame(av_list[[5]][1]),aes(x=av_list[[5]][[1]][,1],y=av_list[[5]][[1]][,2]))+geom_point()+theme_light()+geom_smooth(method="lm", se=TRUE,color="black")+xlab("")+ylab("")+theme(axis.text.x = element_text(size=15),axis.text.y=element_text(size=15))
wtend_b<-ggplot(data=as.data.frame(av_list[[5]][2]),aes(x=av_list[[5]][[2]][,1],y=av_list[[5]][[2]][,2]))+geom_point()+theme_light()+geom_smooth(method="lm", se=TRUE,color="black")+xlab("")+ylab("")+theme(axis.text.x = element_text(size=15),axis.text.y=element_text(size=15))

wtpeak<-ggplot(data=as.data.frame(av_list[[6]][1]),aes(x=av_list[[6]][[1]][,1],y=av_list[[6]][[1]][,2]))+geom_point()+theme_light()+geom_smooth(method="lm", se=TRUE,color="black")+xlab("")+ylab("")+theme(axis.text.x = element_text(size=15),axis.text.y=element_text(size=15))
wtpeak_b<-ggplot(data=as.data.frame(av_list[[6]][2]),aes(x=av_list[[6]][[2]][,1],y=av_list[[6]][[2]][,2]))+geom_point()+theme_light()+geom_smooth(method="lm", se=TRUE,color="black")+xlab("")+ylab("")+theme(axis.text.x = element_text(size=15),axis.text.y=element_text(size=15))
wtpeak_c<-ggplot(data=as.data.frame(av_list[[6]][3]),aes(x=av_list[[6]][[3]][,1],y=av_list[[6]][[3]][,2]))+geom_point()+theme_light()+geom_smooth(method="lm", se=TRUE,color="black")+xlab("")+ylab("")+theme(axis.text.x = element_text(size=15),axis.text.y=element_text(size=15))

cumabun<-ggplot(data=as.data.frame(av_list[[7]][1]),aes(x=av_list[[7]][[1]][,1],y=av_list[[7]][[1]][,2]))+geom_point()+theme_light()+geom_smooth(method="lm", se=TRUE,color="black")+xlab("")+ylab("")+theme(axis.text.x = element_text(size=15),axis.text.y=element_text(size=15))

peakabun<-ggplot(data=as.data.frame(av_list[[8]][1]),aes(x=av_list[[8]][[1]][,1],y=av_list[[8]][[1]][,2]))+geom_point()+theme_light()+geom_smooth(method="lm", se=TRUE,color="black")+xlab("")+ylab("")+theme(axis.text.x = element_text(size=15),axis.text.y=element_text(size=15))


#Plot Regression with year
pl<-ggplot(data=Phenomerged,aes(x=Year,y=PRECIP_LATE))+theme_light()+geom_smooth(method="lm",se=TRUE,color="black")+geom_point()+ylab("")+xlab("")+stat_poly_eq(use_label(c("R2","p")))+theme(axis.text.x=element_text(size=12),axis.text.y=element_text(size=12),text=element_text(family="Times New Roman", face="bold", size=12))
pe<-ggplot(data=Phenomerged,aes(x=Year,y=PRECIP_EARLY))+theme_light()+geom_smooth(method="lm",se=TRUE,color="black")+geom_point()+ylab("")+xlab("")+stat_poly_eq(use_label(c("R2","p")))+theme(axis.text.x=element_text(size=12),axis.text.y=element_text(size=12),text=element_text(family="Times New Roman", face="bold", size=12))
ge<-ggplot(data=Phenomerged,aes(x=Year,y=GDDEARLY))+theme_light()+geom_smooth(method="lm",se=TRUE,color="black")+geom_point()+ylab("")+xlab("")+stat_poly_eq(use_label(c("R2","p")))+theme(axis.text.x=element_text(size=12),axis.text.y=element_text(size=12),text=element_text(family="Times New Roman", face="bold", size=12))
gl<-ggplot(data=Phenomerged,aes(x=Year,y=GDDLATE))+theme_light()+geom_smooth(method="lm",se=TRUE,color="black")+geom_point()+ylab("")+xlab("")+stat_poly_eq(use_label(c("R2","p")))+theme(axis.text.x=element_text(size=12),axis.text.y=element_text(size=12),text=element_text(family="Times New Roman", face="bold", size=12))
mfstart<-ggplot(data=Phenomerged,aes(x=Year,y=MFStart))+theme_light()+geom_smooth(method="lm",se=TRUE,color="black")+geom_point()+ylab("")+xlab("")+stat_poly_eq(use_label(c("R2","p")))+theme(axis.text.x=element_text(size=12),axis.text.y=element_text(size=12),text=element_text(family="Times New Roman", face="bold", size=12))
mfend<-ggplot(data=Phenomerged,aes(x=Year,y=MFEnd))+theme_light()+geom_smooth(method="lm",se=TRUE,color="black")+geom_point()+ylab("")+xlab("")+stat_poly_eq(use_label(c("R2","p")))+theme(axis.text.x=element_text(size=12),axis.text.y=element_text(size=12),text=element_text(family="Times New Roman", face="bold", size=12))
mfpeak<-ggplot(data=Phenomerged,aes(x=Year,y=MFPeak))+theme_light()+geom_smooth(method="lm",se=TRUE,color="black")+geom_point()+ylab("")+xlab("")+stat_poly_eq(use_label(c("R2","p")))+theme(axis.text.x=element_text(size=12),axis.text.y=element_text(size=12),text=element_text(family="Times New Roman", face="bold", size=12))
mflength<-ggplot(data=Phenomerged,aes(x=Year,y=MFLength))+theme_light()+geom_smooth(method="lm",se=TRUE,color="black")+geom_point()+ylab("")+xlab("")+stat_poly_eq(use_label(c("R2","p")))+theme(axis.text.x=element_text(size=12),axis.text.y=element_text(size=12),text=element_text(family="Times New Roman", face="bold", size=12))
wstart<-ggplot(data=Phenomerged,aes(x=Year,y=WTCStart))+theme_light()+geom_smooth(method="lm",se=TRUE,color="black")+geom_point()+ylab("")+xlab("")+stat_poly_eq(use_label(c("R2","p")))+theme(axis.text.x=element_text(size=12),axis.text.y=element_text(size=12),text=element_text(family="Times New Roman", face="bold", size=12))
wend<-ggplot(data=Phenomerged,aes(x=Year,y=WTCEnd))+theme_light()+geom_smooth(method="lm",se=TRUE,color="black")+geom_point()+ylab("")+xlab("")+stat_poly_eq(use_label(c("R2","p")))+theme(axis.text.x=element_text(size=12),axis.text.y=element_text(size=12),text=element_text(family="Times New Roman", face="bold", size=12))
wpeak<-ggplot(data=Phenomerged,aes(x=Year,y=WTCPeak))+theme_light()+geom_smooth(method="lm",se=TRUE,color="black")+geom_point()+ylab("")+xlab("")+stat_poly_eq(use_label(c("R2","p")))+theme(axis.text.x=element_text(size=12),axis.text.y=element_text(size=12),text=element_text(family="Times New Roman", face="bold", size=12))
wlength<-ggplot(data=Phenomerged,aes(x=Year,y=WTFLength))+theme_light()+geom_smooth(method="lm",se=TRUE,color="black")+geom_point()+ylab("")+xlab("")+stat_poly_eq(use_label(c("R2","p")))+theme(axis.text.x=element_text(size=12),axis.text.y=element_text(size=12),text=element_text(family="Times New Roman", face="bold", size=12))
mcum<-ggplot(data=Phenomerged,aes(x=Year,y=CummulativeAbundanceMF))+theme_light()+geom_smooth(method="lm",se=TRUE,color="black")+geom_point()+ylab("")+xlab("")+stat_poly_eq(use_label(c("R2","p")))+theme(axis.text.x=element_text(size=12),axis.text.y=element_text(size=12),text=element_text(family="Times New Roman", face="bold", size=12))
mpa<-ggplot(data=Phenomerged,aes(x=Year,y=PeakAbundanceMF))+theme_light()+geom_smooth(method="lm",se=TRUE,color="black")+geom_point()+ylab("")+xlab("")+stat_poly_eq(use_label(c("R2","p")))+theme(axis.text.x=element_text(size=12),axis.text.y=element_text(size=12),text=element_text(family="Times New Roman", face="bold", size=12))
wcum<-ggplot(data=Phenomerged,aes(x=Year,y=CummulativeAbundanceWTC))+theme_light()+geom_smooth(method="lm",se=TRUE,color="black")+geom_point()+ylab("")+xlab("")+stat_poly_eq(use_label(c("R2","p")))+theme(axis.text.x=element_text(size=12),axis.text.y=element_text(size=12),text=element_text(family="Times New Roman", face="bold", size=12))
wpa<-ggplot(data=Phenomerged,aes(x=Year,y=PeakAbundanceWTC))+theme_light()+geom_smooth(method="lm",se=TRUE,color="black")+geom_point()+ylab("")+xlab("")+stat_poly_eq(use_label(c("R2","p")))+theme(axis.text.x=element_text(size=12),axis.text.y=element_text(size=12),text=element_text(family="Times New Roman", face="bold", size=12))
mfavg<-ggplot(data=Phenomerged,aes(x=Year,y=Avearge_MF_plants))+theme_light()+geom_smooth(method="lm",se=TRUE,color="black")+geom_point()+ylab("")+xlab("")+stat_poly_eq(use_label(c("R2","p")))+theme(axis.text.x=element_text(size=12),axis.text.y=element_text(size=12),text=element_text(family="Times New Roman", face="bold", size=12))
mfavgOP<-ggplot(data=Phenomerged,aes(x=Year,y=Average_MF_plants_ONLYIFPRESENT))+theme_light()+geom_smooth(method="lm",se=TRUE,color="black")+geom_point()+ylab("")+xlab("")+stat_poly_eq(use_label(c("R2","p")))+theme(axis.text.x=element_text(size=12),axis.text.y=element_text(size=12),text=element_text(family="Times New Roman", face="bold", size=12))
mffavg<-ggplot(data=Phenomerged,aes(x=Year,y=Average_Flower_MF))+theme_light()+geom_smooth(method="lm",se=TRUE,color="black")+geom_point()+ylab("")+xlab("")+stat_poly_eq(use_label(c("R2","p")))+theme(axis.text.x=element_text(size=12),axis.text.y=element_text(size=12),text=element_text(family="Times New Roman", face="bold", size=12))
mffavgOP<-ggplot(data=Phenomerged,aes(x=Year,y=Average_Flower_MF_IFPRESENT))+theme_light()+geom_smooth(method="lm",se=TRUE,color="black")+geom_point()+ylab("")+xlab("")+stat_poly_eq(use_label(c("R2","p")))+theme(axis.text.x=element_text(size=12),axis.text.y=element_text(size=12),text=element_text(family="Times New Roman", face="bold", size=12))
mfsavg<-ggplot(data=Phenomerged,aes(x=Year,y=Average_Seed_MF))+theme_light()+geom_smooth(method="lm",se=TRUE,color="black")+geom_point()+ylab("")+xlab("")+stat_poly_eq(use_label(c("R2","p")))+theme(axis.text.x=element_text(size=12),axis.text.y=element_text(size=12),text=element_text(family="Times New Roman", face="bold", size=12))
mfsavgOP<-ggplot(data=Phenomerged,aes(x=Year,y=SEED_MF_IFPRESENT))+theme_light()+geom_smooth(method="lm",se=TRUE,color="black")+geom_point()+ylab("")+xlab("")+stat_poly_eq(use_label(c("R2","p")))+theme(axis.text.x=element_text(size=12),axis.text.y=element_text(size=12),text=element_text(family="Times New Roman", face="bold", size=12))
mfoccupy<-ggplot(data=Phenomerged,aes(x=Year,y=AVGplots_occurpied_MF))+theme_light()+geom_smooth(method="lm",se=TRUE,color="black")+geom_point()+ylab("")+xlab("")+stat_poly_eq(use_label(c("R2","p")))+theme(axis.text.x=element_text(size=12),axis.text.y=element_text(size=12),text=element_text(family="Times New Roman", face="bold", size=12))
wtoccupy<-ggplot(data=Phenomerged,aes(x=Year,y=Quadrat_AvgWTC))+theme_light()+geom_smooth(method="lm",se=TRUE,color="black")+geom_point()+ylab("")+xlab("")+stat_poly_eq(use_label(c("R2","p")))+theme(axis.text.x=element_text(size=12),axis.text.y=element_text(size=12),text=element_text(family="Times New Roman", face="bold", size=12))
wtavg<-ggplot(data=Phenomerged,aes(x=Year,y=Average_WTC_Plant))+theme_light()+geom_smooth(method="lm",se=TRUE,color="black")+geom_point()+ylab("")+xlab("")+stat_poly_eq(use_label(c("R2","p")))+theme(axis.text.x=element_text(size=12),axis.text.y=element_text(size=12),text=element_text(family="Times New Roman", face="bold", size=12))
wtavgOP<-ggplot(data=Phenomerged,aes(x=Year,y=Average_WTC_plant_IFPRESENT))+theme_light()+geom_smooth(method="lm",se=TRUE,color="black")+geom_point()+ylab("")+xlab("")+stat_poly_eq(use_label(c("R2","p")))+theme(axis.text.x=element_text(size=12),axis.text.y=element_text(size=12),text=element_text(family="Times New Roman", face="bold", size=12))
wtfavg<-ggplot(data=Phenomerged,aes(x=Year,y=Average_WTC_flower))+theme_light()+geom_smooth(method="lm",se=TRUE,color="black")+geom_point()+ylab("")+xlab("")+stat_poly_eq(use_label(c("R2","p")))+theme(axis.text.x=element_text(size=12),axis.text.y=element_text(size=12),text=element_text(family="Times New Roman", face="bold", size=12))
wtfavgOP<-ggplot(data=Phenomerged,aes(x=Year,y=Average_WTC_flower_IFPRESENT))+theme_light()+geom_smooth(method="lm",se=TRUE,color="black")+geom_point()+ylab("")+xlab("")+stat_poly_eq(use_label(c("R2","p")))+theme(axis.text.x=element_text(size=12),axis.text.y=element_text(size=12),text=element_text(family="Times New Roman", face="bold", size=12))
wtsavg<-ggplot(data=Phenomerged,aes(x=Year,y=Average_WTC_seed))+theme_light()+geom_smooth(method="lm",se=TRUE,color="black")+geom_point()+ylab("")+xlab("")+stat_poly_eq(use_label(c("R2","p")))+theme(axis.text.x=element_text(size=12),axis.text.y=element_text(size=12),text=element_text(family="Times New Roman", face="bold", size=12))
wtsavgOP<-ggplot(data=Phenomerged,aes(x=Year,y=Average_WTC_seed_IFPRESENT))+theme_light()+geom_smooth(method="lm",se=TRUE,color="black")+geom_point()+ylab("")+xlab("")+stat_poly_eq(use_label(c("R2","p")))+theme(axis.text.x=element_text(size=12),axis.text.y=element_text(size=12),text=element_text(family="Times New Roman", face="bold", size=12))

density<-ggarrange(mfavg,wtavg,mfavgOP,wtavgOP,mffavg,wtfavg,mffavgOP,wtfavgOP,mfsavg,wtsavg,mfsavgOP,wtsavgOP,ncol=2,nrow=6)
pheno<-ggarrange(mfstart,wstart,mfend,wend,mfpeak,wpeak,mflength,wlength,mcum,wcum,mpa,wpa,ncol=2,nrow=6)


#CIMIS Sites
early<-ggplot(data=tab_early,aes(x=SITE,y=Accumulated_PRECIP))+geom_boxplot()+xlab("")+ylab("Rainfall (mm)")+theme_light()+theme(text=element_text(family="Times New Roman", face="bold", size=12),axis.title.y=element_text(size=12))
late<-ggplot(data=tab_late,aes(x=SITE,y=Accumulated_PRECIP))+geom_boxplot()+xlab("")+ylab("Rainfall (mm)")+theme_light()+theme(text=element_text(family="Times New Roman", face="bold", size=12),axis.title.y=element_text(size=12))

#Distance From Edge 
wp1l<-ggplot(data=Book8[Book8$Pool==1 & Book8$Ordinance=="WE" & Book8$DistanceFromEdge<=8.1,],aes(x=DistanceFromEdge,y=WTCLength))+theme_light()+xlab("")+ylab("Flowering Length")+scale_x_continuous(breaks=seq(0,12,1))+theme(axis.text.x = element_text(size=12),axis.text.y = element_text(size=12),axis.title.y = element_text(size=12),text=element_text(family="Times New Roman", face="bold", size=12))+geom_line()
wp2l<-ggplot(data=Book8[Book8$Pool==2 & Book8$Ordinance=="WE" & Book8$DistanceFromEdge<=4.9,],aes(x=DistanceFromEdge,y=WTCLength))+theme_light()+xlab("")+ylab("Flowering Length")+scale_x_continuous(breaks=seq(0,12,1))+theme(axis.text.x = element_text(size=12),axis.text.y = element_text(size=12),axis.title.y = element_text(size=12),text=element_text(family="Times New Roman", face="bold", size=12))+geom_line()
wp3l<-ggplot(data=Book8[Book8$Pool==3 & Book8$Ordinance=="WE",],aes(x=DistanceFromEdge,y=WTCLength))+theme_light()+xlab("Distance From Edge (m)")+ylab("Flowering Length")+scale_x_continuous(breaks=seq(0,12,1))+theme(axis.text.x = element_text(size=12),axis.text.y = element_text(size=12),axis.title.y = element_text(size=12),axis.title.x=element_text(size=12),text=element_text(family="Times New Roman", face="bold", size=12))+geom_line()

p3<-ggplot(data=Book8[Book8$Pool==3 & Book8$Ordinance=="WE",],aes(x=DistanceFromEdge))+geom_line(aes(y=MFStart),color="blue")+geom_line(aes(y=MFEnd),color="red")+theme_light()+xlab("Distance From Edge (m)")+ylab("Julian Date")+scale_x_continuous(breaks=seq(0,12,1))+geom_ribbon(aes(ymin=MFStart,ymax=MFEnd), fill="lightgrey", alpha=0.5)+theme(axis.text.x = element_text(size=12),axis.text.y = element_text(size=12),axis.title.y = element_text(size=12),axis.title.x=element_text(size=12),text=element_text(family="Times New Roman", face="bold", size=12))

p3l<-ggplot(data=Book8[Book8$Pool==3 & Book8$Ordinance=="WE",],aes(x=DistanceFromEdge,y=MeanLength))+theme_light()+xlab("Distance From Edge (m)")+ylab("Flowering Length")+scale_x_continuous(breaks=seq(0,12,1))+theme(axis.text.x = element_text(size=12),axis.text.y = element_text(size=12),axis.title.y = element_text(size=12),axis.title.x = element_text(size=12),text=element_text(family="Times New Roman", face="bold", size=12))+geom_line()

wp3<-ggplot(data=Book8[Book8$Pool==3 & Book8$Ordinance=="WE",],aes(x=DistanceFromEdge))+geom_line(aes(y=WTCStart),color="blue")+geom_line(aes(y=WTCEnd),color="red")+theme_light()+xlab("Distance From Edge (m)")+ylab("Julian Date")+scale_x_continuous(breaks=seq(0,12,1))+geom_ribbon(aes(ymin=WTCStart,ymax=WTCEnd), fill="lightgrey", alpha=0.5)+theme(axis.text.x = element_text(size=12),axis.text.y = element_text(size=12),axis.title.y = element_text(size=12),axis.title.x = element_text(size=12),text=element_text(family="Times New Roman", face="bold", size=12))

wp2<-ggplot(data=Book8[Book8$Pool==2 & Book8$Ordinance=="NS" & Book8$DistanceFromEdge<=4.9,],aes(x=DistanceFromEdge))+geom_line(aes(y=WTCStart),color="blue")+geom_line(aes(y=WTCEnd),color="red")+theme_light()+xlab("")+ylab("Julian Date")+scale_x_continuous(breaks=seq(0,12,1))+geom_ribbon(aes(ymin=WTCStart,ymax=WTCEnd), fill="lightgrey", alpha=0.5)+theme(axis.text.x = element_text(size=12),axis.text.y = element_text(size=12),axis.title.y = element_text(size=12),text=element_text(family="Times New Roman", face="bold", size=12))

wp1<-ggplot(data=Book8[Book8$Pool==1 & Book8$Ordinance=="WE" & Book8$DistanceFromEdge<=8.1,],aes(x=DistanceFromEdge))+geom_line(aes(y=WTCStart),color="blue")+geom_line(aes(y=WTCEnd),color="red")+theme_light()+xlab("")+ylab("Julian Date")+scale_x_continuous(breaks=seq(0,12,1))+geom_ribbon(aes(ymin=WTCStart,ymax=WTCEnd), fill="lightgrey", alpha=0.5)+theme(axis.text.x = element_text(size=12),axis.text.y = element_text(size=12),axis.title.y = element_text(size=12),text=element_text(family="Times New Roman", face="bold", size=12))

p1l<-ggplot(data=Book8[Book8$Pool==1 & Book8$Ordinance=="WE",],aes(x=DistanceFromEdge,y=MeanLength))+theme_light()+xlab("")+ylab("Flowering Length")+scale_x_continuous(breaks=seq(0,12,1))+theme(axis.text.x = element_text(size=12),axis.text.y = element_text(size=12),axis.title.y = element_text(size=12),text=element_text(family="Times New Roman", face="bold", size=12))+geom_line()

p2l<-ggplot(data=Book8[Book8$Pool==2 & Book8$Ordinance=="NS",],aes(x=DistanceFromEdge,y=MeanLength))+theme_light()+xlab("")+ylab("Flowering Length")+scale_x_continuous(breaks=seq(0,12,1))+theme(axis.text.x = element_text(size=12),axis.text.y = element_text(size=12),axis.title.y = element_text(size=12),text=element_text(family="Times New Roman", face="bold", size=12))+geom_line()

p1<-ggplot(data=Book8[Book8$Pool==1 & Book8$Ordinance=="WE",],aes(x=DistanceFromEdge))+geom_line(aes(y=MFStart),color="blue")+geom_line(aes(y=MFEnd),color="red")+theme_light()+xlab("")+ylab("Julian Date")+scale_x_continuous(breaks=seq(0,12,1))+geom_ribbon(aes(ymin=MFStart,ymax=MFEnd), fill="lightgrey", alpha=0.5)+theme(axis.text.x = element_text(size=12),axis.text.y = element_text(size=12),axis.title.y = element_text(size=12),text=element_text(family="Times New Roman", face="bold", size=12))

p2<-ggplot(data=Book8[Book8$Pool==2 & Book8$Ordinance=="NS",],aes(x=DistanceFromEdge))+geom_line(aes(y=MFStart),color="blue")+geom_line(aes(y=MFEnd),color="red")+theme_light()+xlab("")+ylab("Julian Date")+scale_x_continuous(breaks=seq(0,12,1))+geom_ribbon(aes(ymin=MFStart,ymax=MFEnd), fill="lightgrey", alpha=0.5)+theme(axis.text.x = element_text(size=12),axis.text.y = element_text(size=12),axis.title.y = element_text(size=12),text=element_text(family="Times New Roman", face="bold", size=12))

pool<-ggarrange(p1,p1l,wp1,wp1l,p2,p2l,wp2,wp2l,p3,p3l,wp3,wp3l,ncol=4,nrow=3)

#Pool Boxplots
Phenomerged$Pool<-as.factor(Phenomerged$Pool)
wtstart<-ggplot(data=Phenomerged,aes(x=Pool,y=WTCStart))+geom_boxplot()+xlab("")+ylab("Julian Date")+theme_light()+theme(axis.text.x = element_text(size=15),axis.text.y=element_text(size=15),axis.title.y = element_text(size=15),plot.title = element_text(hjust = 0.5,size=15),text=element_text(family="Times New Roman", face="bold", size=12))+ggtitle("Floral Initiation Date")
wtend<-ggplot(data=Phenomerged,aes(x=Pool,y=WTCEnd))+geom_boxplot()+xlab("")+ylab("Julian Date")+theme_light()+theme(axis.text.x = element_text(size=15),axis.text.y=element_text(size=15),axis.title.y = element_text(size=15),plot.title = element_text(hjust = 0.5,size=15),text=element_text(family="Times New Roman", face="bold", size=12))+ggtitle("Floral Termination Date")
wtpeak<-ggplot(data=Phenomerged,aes(x=Pool,y=WTCPeak))+geom_boxplot()+xlab("")+ylab("Julian Date")+theme_light()+theme(axis.text.x = element_text(size=15),axis.text.y=element_text(size=15),axis.title.y = element_text(size=15),plot.title = element_text(hjust = 0.5,size=15),text=element_text(family="Times New Roman", face="bold", size=12))+ggtitle("Peak Flowering Date")
wtlength<-ggplot(data=Phenomerged,aes(x=Pool,y=WTFLength))+geom_boxplot()+xlab("")+ylab("Day")+theme_light()+theme(axis.text.x = element_text(size=15),axis.text.y=element_text(size=15),axis.title.y = element_text(size=15),plot.title = element_text(hjust = 0.5,size=15),text=element_text(family="Times New Roman", face="bold", size=12))+ggtitle("Flowering Length")

mfstart<-ggplot(data=Phenomerged,aes(x=Pool,y=MFStart))+geom_boxplot()+xlab("")+ylab("Julian Date")+theme_light()+theme(axis.text.x = element_text(size=15),axis.text.y=element_text(size=15),axis.title.y = element_text(size=15),plot.title = element_text(hjust = 0.5,size=15),text=element_text(family="Times New Roman", face="bold", size=12))+ggtitle("Floral Initiation Date")
mfend<-ggplot(data=Phenomerged,aes(x=Pool,y=MFEnd))+geom_boxplot()+xlab("")+ylab("Julian Date")+theme_light()+theme(axis.text.x = element_text(size=15),axis.text.y=element_text(size=15),axis.title.y = element_text(size=15),plot.title = element_text(hjust = 0.5,size=15),text=element_text(family="Times New Roman", face="bold", size=12))+ggtitle("Floral Termination Date")
mfpeak<-ggplot(data=Phenomerged,aes(x=Pool,y=MFPeak))+geom_boxplot()+xlab("")+ylab("Julian Date")+theme_light()+theme(axis.text.x = element_text(size=15),axis.text.y=element_text(size=15),axis.title.y = element_text(size=15),plot.title = element_text(hjust = 0.5,size=15),text=element_text(family="Times New Roman", face="bold", size=12))+ggtitle("Peak Flowering Date")
mflength<-ggplot(data=Phenomerged,aes(x=Pool,y=MFLength))+geom_boxplot()+xlab("")+ylab("Day")+theme_light()+theme(axis.text.x = element_text(size=15),axis.text.y=element_text(size=15),axis.title.y = element_text(size=15),plot.title = element_text(hjust = 0.5,size=15),text=element_text(family="Times New Roman", face="bold", size=12))+ggtitle("Flowering Length")

wtcum<-ggplot(data=Phenomerged,aes(x=Pool,y=CummulativeAbundanceWTC))+geom_boxplot()+xlab("")+ylab("Number of Flowers")+theme_light()+theme(axis.text.x = element_text(size=15),axis.text.y=element_text(size=15),axis.title.y = element_text(size=15),plot.title = element_text(hjust = 0.5,size=15),text=element_text(family="Times New Roman", face="bold", size=12))+ggtitle("Cumulative Abundance")
wtpabun<-ggplot(data=Phenomerged,aes(x=Pool,y=PeakAbundanceWTC))+geom_boxplot()+xlab("")+ylab("Number of Flowers")+theme_light()+theme(axis.text.x = element_text(size=15),axis.text.y=element_text(size=15),axis.title.y = element_text(size=15),plot.title = element_text(hjust = 0.5,size=15),text=element_text(family="Times New Roman", face="bold", size=12))+ggtitle("Peak Abundance")
mfcum<-ggplot(data=Phenomerged,aes(x=Pool,y=CummulativeAbundanceMF))+geom_boxplot()+xlab("")+ylab("Number of Flowers")+theme_light()+theme(axis.text.x = element_text(size=15),axis.text.y=element_text(size=15),axis.title.y = element_text(size=15),plot.title = element_text(hjust = 0.5,size=15),text=element_text(family="Times New Roman", face="bold", size=12))+ggtitle("Cumulative Abundance")
mfpabun<-ggplot(data=Phenomerged,aes(x=Pool,y=PeakAbundanceMF))+geom_boxplot()+xlab("")+ylab("Number of Flowers")+theme_light()+theme(axis.text.x = element_text(size=15),axis.text.y=element_text(size=15),axis.title.y = element_text(size=15),plot.title = element_text(hjust = 0.5,size=15),text=element_text(family="Times New Roman", face="bold", size=12))+ggtitle("Peak Abundance")

mfbox<-ggarrange(mfstart,mfend,mfpeak,mflength,mfcum,mfpabun,ncol=2,nrow=3)
wtbox<-ggarrange(wtstart,wtend,wtpeak,wtlength,wtcum,wtpabun,ncol=2,nrow=3)
