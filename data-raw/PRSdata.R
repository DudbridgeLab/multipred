# Simulate polygenic risk scores (PRS) for six common diseases
# Model from table 1 of Dudbridge F (2020) Stat Methods Med Research
#
# T2D = Type-2 diabetes
# CAD = Coronary artery disease
# CD = Crohn's disease
# UC = Ulcerative colitis
# SCZ = schizophrenia
# RA = Rheumatoid Arthritis

library(mvtnorm)

# Reported area under ROC curve for published PRS
AUC = c(0.66,0.623,0.75,0.7,0.62,0.7)
names(AUC) = c("T2D","CAD","CD","UC","SCZ","RA")

# Reported prevalences
prevalence = c(0.102,0.0461,0.005,0.0025,0.01,0.01)
names(prevalence) = names(AUC)

# Liability variance explained by the PRS
# Equation 4 of Wray et al, PLoS Genet 2010
tau = -qnorm(prevalence)
i = dnorm(tau) / prevalence
v = -dnorm(tau) / (1-prevalence)
q = qnorm(AUC)
R2 = 2*q^2 / ((v-i)^2 + q^2*i*(i-tau) + v*(v-tau))
# R2 via observed scale from my paper
R2 = c(0.0856,0.0398,0.103,0.0553,0.0254,0.0732)
names(R2) = names(AUC)


# Variance-covariance matrix between PRS
VX = matrix(nrow=6,ncol=6)
VX[1,] = c(R2[1],0.0225,-0.0111,-0.0086,-0.00131,-0.038)
VX[2,] = c(VX[1,2],R2[2],0.0347,0.0191,0,-0.034)
VX[3,] = c(VX[1:2,3],R2[3],0.0409,0.00679,-0.00251)
VX[4,] = c(VX[1:3,4],R2[4],0.0048,0.00566)
VX[5,] = c(VX[1:4,5],R2[5],-0.00185)
VX[6,] = c(VX[1:5,6],R2[6])

# Genetic correlations between the traits
VL = matrix(nrow=6,ncol=6)
VL[1,] = c(1,0.384,-0.119,-0.125,-0.028,-0.048)
VL[2,] = c(VL[1,2],1,0.057,0.038,0,-0.063)
VL[3,] = c(VL[1:2,3],1,0.543,0.113,-0.029)
VL[4,] = c(VL[1:3,4],1,0.128,0.089)
VL[5,] = c(VL[1:4,5],1,-0.043)
VL[6,] = c(VL[1:5,6],1)

# SNP heritabilities
h2 = c(0.196,0.22,0.26,0.19,0.235,0.18)
VXh2 = diag(h2)
for(i in 1:6)
  for(j in 1:6)
    if (i!=j)
      VXh2[i,j] = VL[i,j]*sqrt(VXh2[i,i]*VXh2[j,j])

# simulate 10000 liabilities and PRS for these six diseases
PRS = rmvnorm(10000,sigma=rbind(cbind(VL,VX),cbind(VX,VX)))
risk = matrix(nrow=10000,ncol=6)
disease = matrix(nrow=10000,ncol=6)
for(i in 1:6) {
  risk[,i] = pnorm((-qnorm(prevalence[i])-PRS[,i+6])/sqrt(1-VX[i,i]),lower=F)
  disease[,i] = PRS[,i] > (-qnorm(prevalence[i]))
}

PRSdata = list(liability=PRS[,1:6],PRS=PRS[,7:12],risk=risk,disease=disease,VL=VL,VX=VX,h2=h2,prevalence=prevalence)
