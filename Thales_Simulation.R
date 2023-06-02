
# Packages
library(mvtnorm)
source("Thales_Header.R")
library(plyr)
library(sandwich)
library(boot)

# Set the parameters of the simulation - positionn of IT, smaple sizes to test, SD in intake and the experimental diets
IT<-c(100, 100)
n<-c(5, 25, 50)
SD<-10
experiment<-round(exp(seq(log(1/8), log(8), length=20)), 2)

# Convert to composition
diet_comps<-cbind(1-experiment/(1+experiment), experiment/(1+experiment))
colnames(diet_comps)<-c("X", "Y")

# Get the predicted intakes of each nutrients on each diet under CDM
cdm<-CDM(diet_comps, IT)

# Get the total eaten on each diet under CDM
total<-cdm[,2] + cdm[,1]

# Create a VCV, for variance in food intake on each diet - 0 cov across diets (i.e., independent intakes)	
s2<-matrix(0, nrow=length(experiment), ncol=length(experiment))
diag(s2)<-SD^2

# To hold estimated bias
accuracy<-NULL

# To hold coverage from CIs for angle with different methods
coverage_OLM<-NULL
coverage_DM<-NULL
coverage_ILM<-NULL
coverage_HLM<-NULL
coverage_HLMincIT<-NULL
coverage_BS<-NULL

# To hold standard errors
SE_sim<-NULL
SE_OLM<-NULL
SE_DM<-NULL
SE_ILM<-NULL
SE_HLM<-NULL
SE_HLMincIT<-NULL

# To hold mean sample variance of beta across diets
VAR_beta<-NULL

# Number of replicates
reps<-10000

# Loop for sample sizes
for(i in 1:length(n)){
	
	# Record bias in the angle
	accur<-NULL
	
	# Coverage as estimated by different methods for n[i]
	covers_OLM<-NULL
	covers_DM<-NULL
	covers_ILM<-NULL
	covers_HLM<-NULL
	covers_HLMincIT<-NULL
	covers_BS<-NULL
	
	# se for n[i]
	se_OLM<-NULL
	se_DM<-NULL
	se_ILM<-NULL
	se_HLM<-NULL
	se_HLMincIT<-NULL
	
	# sample variance for beta for n[i]
	var_beta<-NULL
	
	# Loop for Reps
	for(r in 1:reps){
	
		# Estimate the IT
		IT_data<-cbind(rnorm(n[i], IT[1], SD), rnorm(n[i], IT[2], SD))
		est_IT<-apply(IT_data, 2, mean)
		se_IT<-apply(IT_data, 2, sd) / sqrt(n[i]-1)
		
		# Get n intake data for each on each diet
		intakes<-rmvnorm(n[i], mean=total, sigma=s2)
		
		# Get the int of x and y for each as multiplication by proportion
		int_y<-NULL
		int_x<-NULL
		for(j in 1:length(experiment)){
			int_yj<-intakes[,j]*(experiment[j] / (1 + experiment[j]))
			int_xj<-intakes[,j] - int_yj
			int_y<-c(int_y, int_yj)
			int_x<-c(int_x, int_xj)
		}
		
		# Create the dataframe
		data<-data.frame(diet=as.factor(sort(rep(experiment, n[i]))))
		data$int_x<-int_x
		data$int_y<-int_y
		
		# Calculate the angles for each, using estimated IT
		data$beta<-calculate_Thales(data[,c("int_y", "int_x")], est_IT[c(2,1)])
		
		# Record sample variance in beta
		var_beta<-rbind(var_beta, ddply(data, .(diet), summarise, var(beta))[,2])
		
		# Calculate using true IT (out of interest for bias)
		data$beta2<-calculate_Thales(data[,c("int_y", "int_x")], IT[c(2,1)])
		
		# Fit the model, get estimates, SEs and CIs using the original method
		model<-lm(beta - 90 ~ 0 + diet, data=data)
		accur<-rbind(accur, model$coef)
		se_OLM<-rbind(se_OLM, sqrt(diag(vcov(model))))
		CIs<-confint(model)	
		covers_OLM<-rbind(covers_OLM, (CIs[,1] < 0 & CIs[,2] > 0))
		rm(CIs)
		
		# Get SE and CI using robust SE
		se_het<-sqrt(diag(vcovHC(model, type="HC4m")))
		se_HLM<-rbind(se_HLM, se_het)
		CIs<-cbind(model$coef - se_het * qt(0.975, df=n[i]-1), model$coef + se_het * qt(0.975, df=n[i]-1))
		covers_HLM<-rbind(covers_HLM, (CIs[,1] < 0 & CIs[,2] > 0))
		rm(CIs)
	
		# Now estimate SEs using the DM, and recalcualte CIs
		ses<-my_SE(int_x=data$int_x, int_y=data$int_y, diet=data$diet, IT_x=est_IT[1], IT_y=est_IT[2])
		mus<-ddply(data, .(diet), summarise, mean(beta)-90)[,2]
		se_DM<-rbind(se_DM, ses$se)
		# CI calculated using t-distribution (same as the confint function)
		CIs<-cbind(mus - ses$se * qt(0.975, df=n[i]-1), mus + ses$se * qt(0.975, df=n[i]-1)) 	
		covers_DM<-rbind(covers_DM, (CIs[,1] < 0 & CIs[,2] > 0))
		rm(CIs)

		# Now fit a lm for each diet independently
		se_ind<-NULL
		CIs<-NULL
		for(j in 1:length(experiment)){
			dat_j<-data[which(data$diet == experiment[j]),]$beta - 90
			lm_j<-lm(dat_j ~ 1)
			se_j<-sqrt(vcov(lm_j))
			se_ind<-c(se_ind, se_j)
			ci_j<-cbind(lm_j$coef - se_j * qt(0.975, df=n[i]-1), lm_j$coef + se_j * qt(0.975, df=n[i]-1))
			CIs<-rbind(CIs, ci_j)
		}
		# Save the se and CIs from this method
		se_ILM<-rbind(se_ILM, se_ind)
		covers_ILM<-rbind(covers_ILM, (CIs[,1] < 0 & CIs[,2] > 0))
		rm(CIs)

		# How does boot strapping do?
		CIs<-NULL
		for(j in 1:length(experiment)){
			dat_j<-data[which(data$diet == experiment[j]),]$beta - 90
			reps_j<-boot(data=dat_j, statistic=boot_mu, R=999)
			boot_j<-boot.ci(reps_j, type="stud")
			ci_j<-cbind(boot_j$stud[4], boot_j$stud[5])
			CIs<-rbind(CIs, ci_j)
		}
		# Save the CIs from this method
		covers_BS<-rbind(covers_BS, (CIs[,1] < 0 & CIs[,2] > 0))
		rm(CIs)

		# Repeat heteroskedastic model for known IT 
		model<-lm(beta2 - 90 ~ 0 + diet, data=data)
		
		# Get SE and CI using robust SE
		se_het<-sqrt(diag(vcovHC(model, type="HC4m")))
		se_HLMincIT<-rbind(se_HLMincIT, se_het)
		CIs<-cbind(model$coef - se_het * qt(0.975, df=n[i]-1), model$coef + se_het * qt(0.975, df=n[i]-1))
		covers_HLMincIT<-rbind(covers_HLMincIT, (CIs[,1] < 0 & CIs[,2] > 0))
		rm(CIs)

	}
	
	# Now we have all simulated data for n[i], we can get the simulated estimate of the sampling variance
	SE_sim<-rbind(SE_sim, apply(accur, 2, sd))
	
	# Now calculate the mean biases with known and estimated ITs
	accuracy<-rbind(accuracy, apply(accur, 2, mean))
	
	# Calculate the mean SE as estimated by the three methods
	SE_OLM<-rbind(SE_OLM, apply(se_OLM, 2, mean))
	SE_DM<-rbind(SE_DM, apply(se_DM, 2, mean, na.rm=T))
	SE_ILM<-rbind(SE_ILM, apply(se_ILM, 2, mean))
	SE_HLM<-rbind(SE_HLM, apply(se_HLM, 2, mean))
	SE_HLMincIT<-rbind(SE_HLMincIT, apply(se_HLMincIT, 2, mean))
	
	# Get coverage for CIs
	coverage_OLM<-rbind(coverage_OLM, apply(covers_OLM, 2, mean))
	coverage_DM<-rbind(coverage_DM, apply(covers_DM, 2, mean, na.rm=T))
	coverage_ILM<-rbind(coverage_ILM, apply(covers_ILM, 2, mean))
	coverage_HLM<-rbind(coverage_HLM, apply(covers_HLM, 2, mean))
	coverage_HLMincIT<-rbind(coverage_HLMincIT, apply(covers_HLMincIT, 2, mean))
	coverage_BS<-rbind(coverage_BS, apply(covers_BS, 2, mean))
	
	# Get mean sample variance across diets
	VAR_beta<-rbind(VAR_beta, apply(var_beta, 2, mean))

}

# Get all the column names labelled right
colnames(coverage_DM) <- colnames(coverage_OLM)
colnames(coverage_ILM) <- colnames(coverage_OLM)
colnames(coverage_BS) <- colnames(coverage_OLM)
colnames(SE_ILM) <- colnames(SE_HLM)
colnames(SE_DM) <- colnames(SE_HLM)
colnames(VAR_beta)<-colnames(SE_HLM)

# Save 
coverage<-as.data.frame(rbind(coverage_BS, coverage_DM,  coverage_HLM, coverage_HLMincIT, coverage_ILM, coverage_OLM))
coverage$n<-rep(n, 6)
coverage$method<-sort(rep(c("Bootstrap", "Delta Method", "Hetero. LM", "Hetero. LM (incl. IT)", "Ind. LM", "Orig. LM"), 3))
write.table(coverage, file="Thales_coverage.csv", sep=",", row.names=F, col.names=names(coverage))

SES<-as.data.frame(rbind(SE_DM, SE_HLM, SE_HLMincIT, SE_ILM, SE_OLM, SE_sim))
SES$n<-rep(n, 6)
SES$method<-sort(rep(c("Delta Method", "Hetero. LM", "Hetero. LM (incl. IT)", "Ind. LM", "Orig. LM", "Simulated"), 3))
write.table(SES, file="Thales_SES.csv", sep=",", row.names=F, col.names=names(SES))

BIAS<-as.data.frame(rbind(accuracy, VAR_beta))
BIAS$n<-rep(n, 2)
BIAS$statistic<-sort(rep(c("Bias_est", "var_beta"), 3))
write.table(BIAS, file="Thales_Bias.csv", sep=",", row.names=F, col.names=names(BIAS))

