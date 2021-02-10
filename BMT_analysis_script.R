rm(list=ls())

#setwd("")

library(corrplot)
library(car)
library(lattice)
library(mgcv)
library(MuMIn)
library(MASS)
library(FactoMineR)
library(factoextra)
require(glmulti)
require(nnet)
require(effects)
require(ggplot2)
require(ggpubr)
require(AER)
require(dplyr)
require(tidyr)
require(modEvA)


##################
# load data
##################

sq<-read.csv("sq_data.csv", header=T)

  # data subset without 2 species of surgeonfishes (to obtain a single species per family)
	sq <- droplevels(sq[-which(sq$genus=="Zebrasoma" | sq$genus=="Acanthurus"), ]); levels(sq$genus)
	rownames(sq) <- NULL # to reset row names after removing rows
	sq %>% group_by(genus) %>% tally()
	head(sq)

# remove 2 extreme observations that are statistical outliers
	c(dim(sq[sq$duration <150,])[1], dim(sq)[1]) # 2 extreme observations
	rownames(sq[sq$duration >150,])
	sq1 <- sq[-c(652,663),] # removes two Pomacanthus sexstriatus
	sq1$clientType <- factor(sq1$clientType); is.factor(sq1$clientType); levels(sq1$clientType)


##################
# descriptive stats and stats for Fig. 1 (validation of response variable for comparative purposes with previous studies)
##################

# number of individuals per spedcies
	sq %>% group_by(genus) %>% tally()

# stats for Table 1

  sq %>% 
	group_by(genus) %>% 
	summarise(total_jolts = sum(jolts), total_duration = sum(duration), mean_duration = mean(duration),
	sd_duration = sd(duration), mean_size = mean(TL), sd_size = sd(TL)) %>% 
	as.data.frame()

# compute stats for Fig. 1

  statsF1 <- sq %>% 
	group_by(genus) %>% 
	summarise(total_jolts = sum(jolts), total_duration = sum(duration), mean_duration = mean(duration),
	sd_duration = sd(duration), mean_size = mean(TL), sd_size = sd(TL)) %>% 
	as.data.frame()
  statsF1$jolts_100sec <- statsF1$total_jolts/statsF1$total_duration

#compute stats for analysis to validate use of jolts (absolute count) rather than jolts per 100s

	sq$Fjolt <- as.factor(ifelse(sq$jolts > 0 , "1+", sq$jolts))
	nojolt <- as.data.frame(table(sq$genus, sq$Fjolt))
	nojoltW <- spread(nojolt, Var2, Freq)
	nojoltW$prop_nojolt <- nojoltW[,2] / (nojoltW[,2] + nojoltW[,3])
	statsF1$prop_nojolt <- nojoltW$prop_nojolt

	p <- ggplot(statsF1, aes(prop_nojolt, jolts_100sec))
	coef(lm(jolts_100sec ~ prop_nojolt, data = statsF1)) # Arothron included
	cor(statsF1$jolts_100sec, statsF1$prop_nojolt, method="pearson") #-0.797
	cor(statsF1$jolts_100sec, statsF1$prop_nojolt, method="spearman") #-0.835
	statsF1.1 <- droplevels(statsF1[-which(statsF1$genus=="Arothron"), ]); levels(statsF1.1$genus) # remove genus Arothron
	coef(lm(jolts_100sec ~ prop_nojolt, data = statsF1.1)) # Arothron excluded
	cor(statsF1.1$jolts_100sec, statsF1.1$prop_nojolt, method="pearson") #-0.773

  p + geom_point() + geom_point(colour = "blue", size = 3) + 
	geom_abline(intercept = 0.4628, slope = -0.5640, colour = "blue") +
	geom_abline(intercept = 0.2056, slope = -0.2118, colour = "red") +
	theme_bw()



##################
# PCAs, variable standardization and corr plots
##################

# fast-start

	fs <- sq1[,c(9:11)]
	head(fs)
	pcafs0 <- PCA(fs, scale.unit=T, ncp=5)
  	get_eigenvalue(pcafs0)
  	varFS <- get_pca_var(pcafs0)
  	head(varFS$coord, 4)
  	fviz_pca_var(pcafs0, col.var = "black")
	sq1$pcFS <- pcafs0$ind$coord[,1] ; head(sq1)

# mucus

	mucus <- sq1[,c(14,16)]
	head(mucus)
	pcmuc0  <- PCA(mucus, scale.unit=T, ncp=5)
  	get_eigenvalue(pcmuc0)
  	varMUC <- get_pca_var(pcmuc0)
  	head(varMUC$coord, 4)
  	plot(sq1[,14], sq1[,16])
  	fviz_pca_var(pcmuc0, col.var = "black")
	sq1$pcMUC <- pcmuc0$ind$coord[,1]; head(sq1)

# predictor standardization

	s_TL <- scale(sq1$TL, scale=T, center=T)
	s_duration <- scale(sq1$duration, scale=T, center=T)
	s_turnRate <- scale(sq1$turnRate, scale=T, center=T)
	s_gnaCm2 <- scale(sq1$gnaCm2, scale=T, center=T)
	s_parCm2 <- scale(sq1$parCm2, scale=T, center=T)
	s_mucusProt <- scale(sq1$mucusProt, scale=T, center=T)

# correlations plots

	sq1$clientType <- as.factor(sq1$clientType)
	pairs(~ TL + duration + escapeDist + velocity + acceleration + turnRate + clientType + mucusWeight + 
	  mucusProt + mucusCal + parCm2 + gnaCm2, data=sq1, panel=panel.smooth, main="Simple Scatterplot Matrix")

	ct1.1<-sq1[,c(3,8:12, 14:18)]; M1<-cor(ct1.1); corrplot(M1, method = "number", type="upper", tl.col = "black")
	ct2<-sq1[,c(3,8,12,14:18)]; M2<-cor(ct2); corrplot(M2, method = "number") # correlations with pcFS
	ct3<-sq1[,c(3,8,12,15,17,18,19,20)]; M3<-cor(ct3); corrplot(M3, method = "number") # correlations with pcFS and pcMUC
	ct4<-data.frame(s_TL, s_duration, s_turnRate, s_gnaCm2, s_parCm2, s_mucusProt, sq1$pcFS, sq1$pcMUC); M4<-cor(ct4); corrplot(M4, method = "number")




##############################
# Visualization of the raw data 
##############################

bwplot(jolts~clientType, data=sq1); aggregate(jolts~clientType, data=sq1, mean)
xyplot(log(jolts+1)~pcMUC,    type=c("p", "smooth"), data=sq1)
xyplot(log(jolts+1)~s_turnRate, type=c("p", "smooth"), data=sq1)
xyplot(log(jolts+1)~s_gnaCm2, type=c("p", "smooth"), data=sq1)
xyplot(log(jolts+1)~pcFS, type=c("p", "smooth"), data=sq1)
xyplot(log(jolts+1)~s_mucusProt,  type=c("p", "smooth"), data=sq1)
xyplot(log(jolts+1)~s_duration,  type=c("p", "smooth"), data=sq1)
xyplot(log(jolts+1)~s_TL, type=c("p", "smooth"), data=sq1)



##################
### full GLMs (to look at diagnostics)
##################

hist(sq1$jolts)

# here, we run a negative binomial GLMM and a Poisson GLMM to evaluate whether one deals with the zeros better than the other

### negative binomial error distribution

  t1 <- glm.nb(jolts ~ I(s_duration^2) + I(s_TL^2) + s_turnRate + pcMUC + s_gnaCm2  + s_parCm2 + s_mucusProt + 
	s_duration + s_TL + clientType + pcFS, data=sq1)

  # diagnostics
	summary(t1) 
      par(mfrow=c(2,2)); plot(t1); par(mfrow=c(1,1))
	table(as.integer(fitted(t1, type="response")), sq1$jolts)
	summary(t1) # Residual deviance:  902.27  on 1041  degrees of freedom
	dispersiontest(t1, trafo=NULL, alternative = "greater") # p-value = 0.925
	dispersiontest(t1, trafo=NULL, alternative = "less") # p-value = 0.074
	X2 <- residuals(t1,type="pearson")^2
	X2 <- sum(X2)
	X2 # 1024.132
	t1$df.residual #1041
	1-pchisq(X2,t1$df.residual) #0.64
  # outliers
	plot(t1)
	qqPlot(resid(t1))
	influenceIndexPlot(t1, vars=c("Cook", "hat"), id.n=3)
  # independence of residuals from linear predictors
	plot(residuals(t1) ~ s_duration)
	plot(residuals(t1) ~ sq1$pcMUC)
	plot(residuals(t1) ~ s_turnRate)
	plot(residuals(t1) ~ s_gnaCm2)
	plot(residuals(t1) ~ s_parCm2)
	plot(residuals(t1) ~ s_mucusProt)
	plot(residuals(t1) ~ sq1$pcFS)
	plot(residuals(t1) ~ s_TL)
	plot(residuals(t1) ~ sq1$clientType)
  # variance explaineed
	r.squaredGLMM(t1) # 0.296
	Dsquared(t1) # 0.277
	AIC(t1)


### poisson error distribution

  t1.1 <- glm(jolts ~ I(s_duration^2) + I(s_TL^2) + s_turnRate + pcMUC + s_gnaCm2  + s_parCm2 + s_mucusProt + 
	s_duration + s_TL + clientType + pcFS, data=sq1, family = poisson)

  # diagnostics
	summary(t1.1) 
      par(mfrow=c(2,2)); plot(t1.1); par(mfrow=c(1,1))
	table(as.integer(fitted(t1.1, type="response")), sq1$jolts)
	summary(t1.1) # Residual deviance:  902.27  on 1041  degrees of freedom
	dispersiontest(t1.1, trafo=NULL, alternative = "greater") # p-value = 0.9254
	dispersiontest(t1.1, trafo=NULL, alternative = "less") # p-value = 0.074
	X2 <- residuals(t1.1,type="pearson")^2
	X2 <- sum(X2)
	X2 # 1024.132
	t1.1$df.residual #1041
	1-pchisq(X2,t1.1$df.residual) #0.64
  # outliers
	plot(t1.1)
	qqPlot(resid(t1.1))
	influenceIndexPlot(t1.1, vars=c("Cook", "hat"), id.n=3)
  # independence of residuals from linear predictors
	plot(residuals(t1.1) ~ s_duration)
	plot(residuals(t1.1) ~ sq1$pcMUC)
	plot(residuals(t1.1) ~ s_turnRate)
	plot(residuals(t1.1) ~ s_gnaCm2)
	plot(residuals(t1.1) ~ s_parCm2)
	plot(residuals(t1.1) ~ s_mucusProt)
	plot(residuals(t1.1) ~ sq1$pcFS)
	plot(residuals(t1.1) ~ s_TL)
	plot(residuals(t1.1) ~ sq1$clientType)
  # variance explaineed
	r.squaredGLMM(t1.1) # 0.2523
	Dsquared(t1.1) # 0.277
	AIC(t1.1)

# N.B. Diagnostics are almost the same for both models. Variance explained is only slightly higher for the negative binomial model
#	     but the AIC is >2 units lower for the Poisson GLMM, so we retain this (simpler) model.



##################
### model selection with AIC (exploratory analysis with all predictors included)
##################

t1.aicc <- glmulti("jolts", c("I(s_duration^2)", "I(s_TL^2)", "pcMUC", "s_turnRate", "s_gnaCm2", "s_parCm2", "s_mucusProt",
        "s_duration", "pcFS", "s_TL", "clientType"), data = sq1, fitfunction = "glm", level = 1, 
	marginality = TRUE, method = "h", crit = "aicc", family = poisson(link = "log"), confsetsize = 1000)

	print(t1.aicc)
	plot(t1.aicc)
	tmp1 <- weightable(t1.aicc); tmp1 <- tmp1[tmp1$aicc <= min(tmp1$aicc) + 4.5,]; tmp1
	plot(t1.aicc, type="s")
	multiest <- data.frame(coef.glmulti(t1.aicc, select="all"))
	plot(multiest$Importance)

# model averaging with MuMIn

	a1 <- glm(jolts ~ clientType + I(s_duration^2) + I(s_TL^2) + pcMUC + s_turnRate + s_gnaCm2 + s_parCm2 + s_mucusProt + s_duration + s_TL, data=sq1, family=poisson)
	a2 <- glm(jolts ~ clientType + I(s_duration^2) + I(s_TL^2) + pcMUC + s_turnRate + s_gnaCm2 + s_mucusProt + s_duration + s_TL, data=sq1, family=poisson)
	a3 <- glm(jolts ~ clientType + I(s_duration^2) + I(s_TL^2) + pcMUC + s_turnRate + s_gnaCm2 + s_parCm2 + s_mucusProt + s_duration + pcFS + s_TL, data=sq1, family=poisson)
	a4 <- glm(jolts ~ clientType + I(s_duration^2) + I(s_TL^2) + pcMUC + s_turnRate + s_gnaCm2 + s_mucusProt + s_duration + pcFS + s_TL, data=sq1, family=poisson)
	a5 <- glm(jolts ~ clientType + I(s_duration^2) + I(s_TL^2) + pcMUC + s_turnRate + s_gnaCm2 + s_parCm2 + s_duration + s_TL, data=sq1, family=poisson)
	a6 <- glm(jolts ~ clientType + I(s_duration^2) + I(s_TL^2) + pcMUC + s_turnRate + s_gnaCm2 + s_parCm2 + s_duration + pcFS + s_TL, data=sq1, family=poisson)
	m.avg <- model.avg(a1, a2, a3, a4, a5, a6)
	summary(m.avg)

### Variance explained (post model selection)

t.final <- glm(jolts ~ I(s_duration^2) + I(s_TL^2) + pcMUC + s_turnRate + s_gnaCm2  + s_parCm2 + 
	s_duration + clientType + s_TL, data=sq1, family=poisson)

	plot(t.final)
	table(as.integer(fitted(t.final, type="response")), sq1$jolts)
	r.squaredGLMM(t.final) # 0.286
	Dsquared(t.final) #0.273

### Plots of model predictions


	plot(effect("s_duration", t.final), colors = "blue", main="", axes=list(y=list(type="response", lab="Jolts", lim=c(0, 9))))
	plot(effect("s_TL", t.final), colors = "blue", main="", axes=list(y=list(type="response", lab="Jolts", lim=c(0, 0.7))))
	plot(effect("pcMUC", t.final), colors = "blue", main="", axes=list(y=list(type="response", lab="Jolts", lim=c(0, 4.5))))
	plot(effect("s_turnRate", t.final), colors = "blue", main="", axes=list(y=list(type="response", lab="Jolts", lim=c(0, 3.2))))
	plot(effect("s_gnaCm2", t.final), colors = "blue", main="", axes=list(y=list(type="response", lab="Jolts", lim=c(0, 0.9))))
	plot(effect("clientType", t.final), colors = "blue", main="", axes=list(y=list(type="response", lab="Jolts", lim=c(0, 2.7))))

### Plots of raw data

	qplot(duration, jolts, data = sq1, colour = I("blue"), size = I(2), alpha = I(0.2), position = position_jitter(w = 0, h = 0.25)) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), panel.border = element_rect(fill=NA, colour = "black"))
	qplot(pcMUC, jolts, data = sq1, colour = I("blue"), size = I(2), alpha = I(0.2), position = position_jitter(w = 0, h = 0.25)) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), panel.border = element_rect(fill=NA, colour = "black"))
	qplot(turnRate, jolts, data = sq1, colour = I("blue"), size = I(2), alpha = I(0.2), position = position_jitter(w = 0, h = 0.25)) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), panel.border = element_rect(fill=NA, colour = "black"))
	qplot(gnaCm2, jolts, data = sq1, colour = I("blue"), size = I(2), alpha = I(0.2), position = position_jitter(w = 0, h = 0.25)) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), panel.border = element_rect(fill=NA, colour = "black"))
	qplot(TL, jolts, data = sq1, colour = I("blue"), size = I(2),alpha = I(0.2), position = position_jitter(w = 0.35, h = 0.25)) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), panel.border = element_rect(fill=NA, colour = "black"))
	qplot(clientType, jolts, data = sq1, colour = I("blue"), size = I(2),alpha = I(0.2), position = position_jitter(w = 0.35, h = 0.25)) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), panel.border = element_rect(fill=NA, colour = "black"))

	
	
	##################
	# Descriptive stats for fast-start escape responses (Tables A1, A2) and mucus characteristics (Table A3)
	##################
	
	FS<-read.csv("raw_data_fastStart_and_mucus.csv", header=TRUE)
  head(FS)
  FS$genus <- as.factor(FS$genus)
  FS <- droplevels(FS[-which(FS$genus=="Zebrasoma" | FS$genus=="Acanthurus"), ]); levels(FS$genus)
  rownames(FS) <- NULL # to reset row names after removing rows

	# count + mean and SD for TL, Desc, Amax, and Umax (Table A1)
	FS %>% 
	  group_by(genus) %>%
	  summarise(n=n(),
	            mean_TL = mean(TL, na.rm=TRUE), sd_TL = sd(TL, na.rm=TRUE),
	            mean_D = mean(Desc, na.rm=T), sd_D = sd(Desc, na.rm=T),
	            mean_U = mean(Umax, na.rm=T), sd_U = sd(Umax, na.rm=T),
	            mean_A = mean(Amax, na.rm=T), sd_A = sd(Amax, na.rm=T)) %>%
	  mutate(sd_TL=round(sd_TL,1), sd_D=round(sd_D,2)) %>%
	  arrange(mean_D, .by_group = TRUE)
	
	# count + mean and SD for TL and Turning rate	 (Table A2)
		FS %>% 
	    drop_na(F_turnRate) %>%
	    group_by(genus) %>%
	    summarise(n=n(),
	            mean_TL = mean(TL, na.rm=TRUE), sd_TL = sd(TL, na.rm=TRUE),
	            mean_FTR = mean(F_turnRate, na.rm=TRUE), sd_FTR = sd(F_turnRate, na.rm=TRUE)) %>%
		  mutate(sd_TL=round(sd_TL,1)) %>%
		  arrange(mean_FTR, .by_group = TRUE)

		# count + mean and SD for TL, mucus weight, protein, and caloric content (Table A3)
		FS %>% 
		  drop_na(relMucusM) %>%
		  group_by(genus) %>%
		  summarise(n=n(),
		            mean_TL = mean(TL, na.rm=TRUE), sd_TL = sd(TL, na.rm=TRUE),
		            mean_Mm = mean(relMucusM, na.rm=T), sd_Mm = sd(relMucusM, na.rm=T),
		            mean_Mp = mean(mucusProt, na.rm=T), sd_Mp = sd(mucusProt, na.rm=T),
		            mean_Mc = mean(calories, na.rm=T), sd_Mc = sd(calories, na.rm=T)) %>%
		  mutate(across(is.numeric, ~ round(., 3))) %>%
		  arrange(mean_Mm, .by_group = TRUE) %>%
		  print.data.frame()
	
##################
# Descriptive stats for ectoparasites (Table A4)
##################
		
	M<-read.csv("raw_data_parasites.csv", header=TRUE)
	head(M)
	M$genus <- as.factor(M$genus)
	M <- droplevels(M[-which(M$genus=="Zebrasoma" | M$genus=="Acanthurus"), ]); levels(M$genus)
	rownames(FS) <- NULL # to reset row names after removing rows
  M$otherParasite <- M$totParasite - M$gna
		
# count + mean and SD for TL, gnathiids, and other ectoparasites (Table A4)
	M %>% 
	  drop_na(gna) %>%
	  group_by(genus) %>%
	  summarise(n=n(),
	            mean_TL = mean(TL, na.rm=TRUE), sd_TL = sd(TL, na.rm=TRUE),
	            mean_gnaCm2 = mean(gnaCm2, na.rm=TRUE), sd_gnaCm2 = sd(gnaCm2, na.rm=TRUE),
	            mean_othCm2 = mean(otherParasiteCm2, na.rm=T), sd_othCm2 = sd(otherParasiteCm2, na.rm=T),
	            mean_gna = mean(gna, na.rm=T), sd_gna = sd(gna, na.rm=T),
	            mean_oth = mean(otherParasite, na.rm=T), sd_oth = sd(otherParasite, na.rm=T)) %>%
	  mutate(across(is.numeric, ~ round(., 3))) %>%
	  arrange(mean_gnaCm2, .by_group = TRUE) %>%
	  print.data.frame()
	