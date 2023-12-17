###################
#Download Packages#
###################
install.packages(c("nortest","readr","sjPlot","ggplot2","lmtest","AICcmodavg","DescTools"))
lapply(c("nortest","readr","sjPlot","ggplot2","lmtest","AICcmodavg","DescTools"), require, character.only = TRUE)
list.of.packages <- c("nortest","readr","sjPlot","ggplot2","lmtest","AICcmodavg","DescTools")
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)
##############
#Upload Files#
##############
survey <- read.csv("~/Documents/Github_Projects/Meadowfoam_WhitetipClover_Phenology/VERNALPOOLPHENOLOGYV2csv.csv")
merced_climate <- read.csv("~/Documents/Github_Projects/Meadowfoam_WhitetipClover_Phenology/MERCED_CLIMATE_2000_2022.csv")
climate_stations <- read.csv("~/Documents/Github_Projects/Meadowfoam_WhitetipClover_Phenology/CIMIS_SITES.csv")

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
growing_season<-GDH(growing_season,10,GDH="GDH_MERCED")
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
model<-aov(Accumulated_PRECIP~SITE,data=tab_sum)
DunnettTest(x=tab_sum$Accumulated_PRECIP,g=tab_sum$SITE)



################################
#Summarize Climate Variables####
################################
climate_vars<-list("GDDEARLY","GDDLATE","MEANEARLY","MEANLATE","MinEARLY","MinLate","MAXEarly","MAXLate","PRECIP_EARLY","PRECIP_LATE")
climate_summary <- list()
for(i in 1:length(climate_vars)){
    name = climate_vars[[i]][1]
    climate_summary[[name]]<-summary(Phenomerged[,name])
}

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
###############################################################
#Run Linear Regression Models: Phenology by Distance from Edge#
###############################################################
distance_lm_list=list()
for(i in names(ClayFlower[,6:11])){
    name = i
    distance_lm_list[[name]]<-lm(as.matrix(ClayFlower[,i])~DistanceFromEdge,data=ClayFlower)
}
clay_lm_list=list()
for(i in names(ClayFlower[,6:11])){
    name = i
    clay_lm_list[[name]]<-lm(as.matrix(ClayFlower[,i])~Clay_Percent,data=ClayFlower)
}

################################
#Analyze Rate of Pool Occupancy#
################################
#Meadowfoam 
mf_lm_list <- list()
for(i in names(VP_Population[,65:75])){
    name = i 
    mf_lm_list[[i]]<-lm(as.matrix(VP_Population[,i])~GDDTotalWinter + PRECIPTotalWinter + MFStart + MFEnd,data=VP_Population)
}
#Whitetip Clover
wc_lm_list <- list()
for(i in names(VP_Population[,76:86])){
    name = i 
    wc_lm_list[[i]]<-lm(as.matrix(VP_Population[,i])~GDDTotalWinter + PRECIPTotalWinter + MFStart + MFEnd,data=VP_Population)
}


