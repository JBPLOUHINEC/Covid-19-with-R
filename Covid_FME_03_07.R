##==============================================================
## PART 1 : Systéme EDO  pour Covid er résolution 
##==============================================================

library(deSolve)
library(gdata)
library(ggplot2)
library(rootSolve)
library(coda)
library(data.table)
library(FME)
library(scales)

# Définition Function MySIRS :

mySIRS <- function(Time, State, pars) {
  with(as.list(c(State, pars)),{
    
    states<- matrix(State, nrow=NUM_COHORTS, ncol=NUM_STATES)
    
    beta <- CE/(CohortPopulations)
    
    Susceptible <- states[,1]
    InfectedAsymptomatic <- states[,2]
    InfectedMild <- states[,3]
    InfectedSerious <- states[,4]
    InfectedCritical <- states[,5]
    Recovered <-states[,6]
    Death <-states[,7]
    
    Lambda <- beta %*% InfectedAsymptomatic
    
    IR <- Lambda*Susceptible*(60-Time)*parameter2
    TR_as_m <- InfectedAsymptomatic/delays_r1
    R_as <- InfectedAsymptomatic /delays_1
    TR_m_ser <- InfectedMild /delays_2
    R_m <- InfectedMild  /delays_2
    TR_ser_crit <- Sigma_i_1*InfectedSerious /delays_r3
    R_ser <- InfectedSerious /delays_3
    TR_crit_death <-Sigma_i_2*InfectedCritical /delays_r4
    R_crit <- InfectedCritical /delays_r4
    R_death <- InfectedCritical *Phi
    
    dS_dt <- -IR
    dIa_dt <- IR-TR_as_m-R_as
    dIm_dt <- TR_as_m-TR_m_ser-R_m
    dIs_dt <- TR_m_ser-TR_ser_crit-R_ser
    dIc_dt <- TR_ser_crit-R_crit-R_death
    dR_dt <- R_as+R_m+R_ser+R_crit
    dDeath_dt <- R_death
    
    
    
    return (list(c(dS_dt, dIa_dt,dIm_dt,dIs_dt,dIc_dt,dR_dt,dDeath_dt)))
  })
}

# Définition des paramètres :
## initial "guess"

## 70 Parameters 

parameter2<- 0.009

### 25 Parameters Matrix

CE_11 <- 1.61
CE_12 <- 0.77
CE_13 <- 0.20
CE_14 <- 0.06
CE_15 <- 0.65
CE_21 <- 1.39
CE_22 <- 1.63
CE_23 <- 0.43
CE_24 <- 0.06
CE_25 <- 0.05
CE_31 <- 1.04
CE_32 <- 1.61
CE_33 <- 0.70
CE_34 <- 0.09
CE_35 <- 0.10
CE_41 <- 0.32
CE_42 <- 0.52
CE_43 <- 0.66
CE_44 <- 0.43
CE_45 <- 0.41
CE_51 <- 0.27
CE_52 <- 0.38
CE_53 <- 0.48
CE_54 <- 0.60
CE_55 <- 0.55

## 10 Parametres CohortPopulation

CohortPopulations_1 <- 625969
CohortPopulations_2 <- 631964
CohortPopulations_3 <- 705269
CohortPopulations_4 <- 688559
CohortPopulations_5 <- 750119
CohortPopulations_6 <- 655169
CohortPopulations_7 <- 570506
CohortPopulations_8 <- 341103
CohortPopulations_9 <- 146632
CohortPopulations_10 <- 73316

## 7 Parametres Délais
  
delays_r1_11 <- 14
delays_1_11 <- 5
delays_r2_11 <- 16
delays_2_11 <- 3
delays_r3_11 <- 12
delays_3_11 <-4
delays_r4_11 <-13

## 17 Parametres Taux Hospitalisations

Sigma_1_1 <-0.001
Sigma_2_1 <-0.0003
Sigma_3_1 <-0.011
Sigma_4_1 <-0.0339
Sigma_5_1 <-0.0429
Sigma_6_1 <-0.0819
Sigma_7_1 <-0.1179
Sigma_8_1 <-0.1659
Sigma_9_1 <-0.1599
Sigma_10_1 <-0.2

Sigma_1_2  <- 0.02
Sigma_4_2 <- 0.05
Sigma_5_2 <- 0.063
Sigma_6_2 <- 0.122
Sigma_7_2 <- 0.274
Sigma_8_2 <- 0.432
Sigma_9_2 <- 0.709
Sigma_10_2 <- 0.709

## 10 Parametres Taux Mortalité

Phi_1 <- 0.0001
Phi_2 <- 0.001
Phi_3 <-0.01
Phi_4 <-0.02
Phi_5 <-0.035
Phi_6 <-0.05
Phi_7 <-0.1
Phi_8 <-.5
Phi_9 <-1.0
Phi_10 <-1.0 

parameter <-5
CE <- matrix (c(CE_44
), nrow=10, ncol=10, byrow=TRUE)/parameter


CohortPopulations <- c(CohortPopulations_1, CohortPopulations_2, CohortPopulations_3,
                       CohortPopulations_4, CohortPopulations_5, CohortPopulations_6,
                       CohortPopulations_7, CohortPopulations_8, CohortPopulations_9,
                       CohortPopulations_10)

delays_r1<- c(delays_r1_11,delays_r1_11,delays_r1_11,delays_r1_11,delays_r1_11,
              delays_r1_11,delays_r1_11,delays_r1_11,delays_r1_11,delays_r1_11 )


delays_1<- c(delays_1_11, delays_1_11,delays_1_11,delays_1_11,delays_1_11,
             delays_1_11, delays_1_11, delays_1_11, delays_1_11,delays_1_11 )

delays_r2<- c(delays_r2_11, delays_r2_11, delays_r2_11, delays_r2_11,delays_r2_11,
              delays_r2_11, delays_r2_11,delays_r2_11,delays_r2_11,delays_r2_11 )

delays_2<- c(delays_2_11, delays_2_11, delays_2_11, delays_2_11,delays_2_11,
             delays_2_11, delays_2_11, delays_2_11,delays_2_11,delays_2_11)

delays_r3<- c(delays_r3_11, delays_r3_11, delays_r3_11, delays_r3_11,delays_r3_11,
              delays_r3_11, delays_r3_11, delays_r3_11, delays_r3_11, delays_r3_11)

delays_3<- c(delays_3_11, delays_3_11,delays_3_11, delays_3_11,delays_3_11,
             delays_3_11, delays_3_11, delays_3_11, delays_3_11,delays_3_11 )

delays_r4<- c(delays_r4_11, delays_r4_11, delays_r4_11, delays_r4_11,delays_r4_11,
              delays_r4_11, delays_r4_11, delays_r4_11, delays_r4_11, delays_r4_11 )

Sigma_i_1 <- c(Sigma_1_1, Sigma_2_1, Sigma_3_1, Sigma_4_1, Sigma_5_1,
               Sigma_6_1, Sigma_7_1, Sigma_8_1, Sigma_9_1, Sigma_10_1 )

Sigma_i_2 <- c(Sigma_1_2, Sigma_1_2, Sigma_1_2, Sigma_4_2, Sigma_5_2,
               Sigma_6_2, Sigma_7_2, Sigma_8_2, Sigma_9_2, Sigma_10_2 )

Phi<- c(Phi_1,Phi_2,Phi_3, Phi_4, Phi_5, Phi_6, Phi_7, Phi_8, Phi_9, Phi_10 )


pars <- c( parameter2)


#### Valeurs Initiales



SusceptibleClasse1_0 <- 625968
SusceptibleClasse2_0 <- 631963
SusceptibleClasse3_0 <- 705268
SusceptibleClasse4_0 <- 688558
SusceptibleClasse5_0 <- 750118
SusceptibleClasse6_0 <- 655168
SusceptibleClasse7_0 <- 570506
SusceptibleClasse8_0 <- 341102
SusceptibleClasse9_0 <- 146631
SusceptibleClasse10_0 <- 73315
InfectedAsymptomaticClasse1_0 <- 10
InfectedAsymptomaticClasse2_0 <- 10
InfectedAsymptomaticClasse3_0 <- 10
InfectedAsymptomaticClasse4_0 <- 10
InfectedAsymptomaticClasse5_0 <- 10
InfectedAsymptomaticClasse6_0 <- 10
InfectedAsymptomaticClasse7_0 <- 10
InfectedAsymptomaticClasse8_0 <- 10
InfectedAsymptomaticClasse9_0 <- 10
InfectedAsymptomaticClasse10_0 <- 10
InfectedMildClasse1_0 <- 1
InfectedMildClasse2_0 <- 1
InfectedMildClasse3_0 <- 1
InfectedMildClasse4_0 <- 1
InfectedMildClasse5_0 <- 1
InfectedMildClasse6_0 <- 1
InfectedMildClasse7_0 <- 1
InfectedMildClasse8_0 <- 1
InfectedMildClasse9_0 <- 1
InfectedMildClasse10_0 <- 1
InfectedSeriousClasse1_0 <- 1
InfectedSeriousClasse2_0 <- 1
InfectedSeriousClasse3_0 <- 1
InfectedSeriousClasse4_0 <- 1
InfectedSeriousClasse5_0 <- 1
InfectedSeriousClasse6_0 <- 1
InfectedSeriousClasse7_0 <- 1
InfectedSeriousClasse8_0 <- 1
InfectedSeriousClasse9_0 <- 1
InfectedSeriousClasse10_0 <- 1
InfectedCriticalClasse1_0 <- 1
InfectedCriticalClasse2_0 <- 1
InfectedCriticalClasse3_0 <- 1
InfectedCriticalClasse4_0 <- 1
InfectedCriticalClasse5_0 <- 1
InfectedCriticalClasse6_0 <- 1
InfectedCriticalClasse7_0 <- 1
InfectedCriticalClasse8_0 <- 1
InfectedCriticalClasse9_0 <- 1
InfectedCriticalClasse10_0 <- 1
RecoveredClasse1_0 <- 1
RecoveredClasse2_0 <- 1
RecoveredClasse3_0 <- 1
RecoveredClasse4_0 <- 1
RecoveredClasse5_0 <- 1
RecoveredClasse6_0 <- 1
RecoveredClasse7_0 <- 1
RecoveredClasse8_0 <- 1
RecoveredClasse9_0 <- 1
RecoveredClasse10_0 <- 1
DeathClasse1_0 <- 1
DeathClasse2_0 <- 1
DeathClasse3_0 <- 1
DeathClasse4_0 <- 1
DeathClasse5_0 <- 1
DeathClasse6_0 <- 1
DeathClasse7_0 <- 1
DeathClasse8_0 <- 1
DeathClasse9_0 <- 1
DeathClasse10_0 <- 1







#### Structures des Paramètres



auxs <-NULL
#Valeurs Initiales

yini <-  c(SusceptibleClasse1=SusceptibleClasse1_0 , SusceptibleClasse2=SusceptibleClasse2_0,SusceptibleClasse3=SusceptibleClasse3_0,
           SusceptibleClasse4=SusceptibleClasse4_0, SusceptibleClasse5=SusceptibleClasse5_0,SusceptibleClasse6=SusceptibleClasse6_0,
           SusceptibleClasse7=SusceptibleClasse7_0, SusceptibleClasse8=SusceptibleClasse8_0,SusceptibleClasse9=SusceptibleClasse9_0,
           SusceptibleClasse10=SusceptibleClasse10_0,
           InfectedAsymptomaticClasse1=InfectedAsymptomaticClasse1_0 ,InfectedAsymptomaticClasse2=InfectedAsymptomaticClasse2_0 ,
           InfectedAsymptomaticClasse3=InfectedAsymptomaticClasse3_0 , InfectedAsymptomaticClasse4=InfectedAsymptomaticClasse4_0 ,
           InfectedAsymptomaticClasse5=InfectedAsymptomaticClasse5_0 ,InfectedAsymptomaticClasse6=InfectedAsymptomaticClasse6_0 ,
           InfectedAsymptomaticClasse7=InfectedAsymptomaticClasse7_0 , InfectedAsymptomaticClasse8=InfectedAsymptomaticClasse8_0 ,
           InfectedAsymptomaticClasse9=InfectedAsymptomaticClasse9_0 ,InfectedAsymptomaticClasse10=InfectedAsymptomaticClasse10_0 ,
           InfectedMildClasse1=InfectedMildClasse1_0,InfectedMildClasse2=InfectedMildClasse2_0,InfectedMildClasse3=InfectedMildClasse3_0,
           InfectedMildClasse4=InfectedMildClasse4_0, InfectedMildClasse5=InfectedMildClasse5_0, InfectedMildClasse6=InfectedMildClasse6_0,
           InfectedMildClasse7=InfectedMildClasse7_0, InfectedMildClasse8=InfectedMildClasse8_0, InfectedMildClasse9=InfectedMildClasse9_0,
           InfectedMildClasse10=InfectedMildClasse10_0,
           InfectedSeriousClasse1=InfectedSeriousClasse1_0,InfectedSeriousClasse2=InfectedSeriousClasse2_0,InfectedSeriousClasse3=InfectedSeriousClasse3_0,
           InfectedSeriousClasse4=InfectedSeriousClasse4_0, InfectedSeriousClasse5=InfectedSeriousClasse5_0, InfectedSeriousClasse6=InfectedSeriousClasse6_0,
           InfectedSeriousClasse7=InfectedSeriousClasse7_0, InfectedSeriousClasse8=InfectedSeriousClasse8_0, InfectedSeriousClasse9=InfectedSeriousClasse9_0,
           InfectedSeriousClasse10=InfectedSeriousClasse10_0,
           InfectedCriticalClasse1=InfectedCriticalClasse1_0, InfectedCriticalClasse2=InfectedCriticalClasse2_0, InfectedCriticalClasse3=InfectedCriticalClasse3_0,
           InfectedCriticalClasse4=InfectedCriticalClasse4_0, InfectedCriticalClasse5=InfectedCriticalClasse5_0, InfectedCriticalClasse6=InfectedCriticalClasse6_0,
           InfectedCriticalClasse7=InfectedCriticalClasse7_0, InfectedCriticalClasse8=InfectedCriticalClasse8_0, InfectedCriticalClasse9=InfectedCriticalClasse9_0,
           InfectedCriticalClasse10=InfectedCriticalClasse10_0,
           RecoveredClasse1=RecoveredClasse1_0, RecoveredClasse2=RecoveredClasse2_0, RecoveredClasse3=RecoveredClasse3_0, RecoveredClasse4=RecoveredClasse4_0,
           RecoveredClasse5=RecoveredClasse5_0,RecoveredClasse6=RecoveredClasse6_0,RecoveredClasse7=RecoveredClasse7_0,RecoveredClasse8=RecoveredClasse8_0,
           RecoveredClasse9=RecoveredClasse9_0,RecoveredClasse10=RecoveredClasse10_0,
           DeathClasse1=DeathClasse1_0, DeathClasse2=DeathClasse2_0, DeathClasse3=DeathClasse3_0, DeathClasse4=DeathClasse4_0, DeathClasse5=DeathClasse5_0,
           DeathClasse6=DeathClasse6_0, DeathClasse7=DeathClasse7_0, DeathClasse8=DeathClasse8_0, DeathClasse9=DeathClasse9_0, DeathClasse10=DeathClasse10_0
)

State <- c(SusceptibleClasse1=SusceptibleClasse1_0 , SusceptibleClasse2=SusceptibleClasse2_0,SusceptibleClasse3=SusceptibleClasse3_0,
           SusceptibleClasse4=SusceptibleClasse4_0, SusceptibleClasse5=SusceptibleClasse5_0,SusceptibleClasse6=SusceptibleClasse6_0,
           SusceptibleClasse7=SusceptibleClasse7_0, SusceptibleClasse8=SusceptibleClasse8_0,SusceptibleClasse9=SusceptibleClasse9_0,
           SusceptibleClasse10=SusceptibleClasse10_0,
           InfectedAsymptomaticClasse1=InfectedAsymptomaticClasse1_0 ,InfectedAsymptomaticClasse2=InfectedAsymptomaticClasse2_0 ,
           InfectedAsymptomaticClasse3=InfectedAsymptomaticClasse3_0 , InfectedAsymptomaticClasse4=InfectedAsymptomaticClasse4_0 ,
           InfectedAsymptomaticClasse5=InfectedAsymptomaticClasse5_0 ,InfectedAsymptomaticClasse6=InfectedAsymptomaticClasse6_0 ,
           InfectedAsymptomaticClasse7=InfectedAsymptomaticClasse7_0 , InfectedAsymptomaticClasse8=InfectedAsymptomaticClasse8_0 ,
           InfectedAsymptomaticClasse9=InfectedAsymptomaticClasse9_0 ,InfectedAsymptomaticClasse10=InfectedAsymptomaticClasse10_0 ,
           InfectedMildClasse1=InfectedMildClasse1_0,InfectedMildClasse2=InfectedMildClasse2_0,InfectedMildClasse3=InfectedMildClasse3_0,
           InfectedMildClasse4=InfectedMildClasse4_0, InfectedMildClasse5=InfectedMildClasse5_0, InfectedMildClasse6=InfectedMildClasse6_0,
           InfectedMildClasse7=InfectedMildClasse7_0, InfectedMildClasse8=InfectedMildClasse8_0, InfectedMildClasse9=InfectedMildClasse9_0,
           InfectedMildClasse10=InfectedMildClasse10_0,
           InfectedSeriousClasse1=InfectedSeriousClasse1_0,InfectedSeriousClasse2=InfectedSeriousClasse2_0,InfectedSeriousClasse3=InfectedSeriousClasse3_0,
           InfectedSeriousClasse4=InfectedSeriousClasse4_0, InfectedSeriousClasse5=InfectedSeriousClasse5_0, InfectedSeriousClasse6=InfectedSeriousClasse6_0,
           InfectedSeriousClasse7=InfectedSeriousClasse7_0, InfectedSeriousClasse8=InfectedSeriousClasse8_0, InfectedSeriousClasse9=InfectedSeriousClasse9_0,
           InfectedSeriousClasse10=InfectedSeriousClasse10_0,
           InfectedCriticalClasse1=InfectedCriticalClasse1_0, InfectedCriticalClasse2=InfectedCriticalClasse2_0, InfectedCriticalClasse3=InfectedCriticalClasse3_0,
           InfectedCriticalClasse4=InfectedCriticalClasse4_0, InfectedCriticalClasse5=InfectedCriticalClasse5_0, InfectedCriticalClasse6=InfectedCriticalClasse6_0,
           InfectedCriticalClasse7=InfectedCriticalClasse7_0, InfectedCriticalClasse8=InfectedCriticalClasse8_0, InfectedCriticalClasse9=InfectedCriticalClasse9_0,
           InfectedCriticalClasse10=InfectedCriticalClasse10_0,
           RecoveredClasse1=RecoveredClasse1_0, RecoveredClasse2=RecoveredClasse2_0, RecoveredClasse3=RecoveredClasse3_0, RecoveredClasse4=RecoveredClasse4_0,
           RecoveredClasse5=RecoveredClasse5_0,RecoveredClasse6=RecoveredClasse6_0,RecoveredClasse7=RecoveredClasse7_0,RecoveredClasse8=RecoveredClasse8_0,
           RecoveredClasse9=RecoveredClasse9_0,RecoveredClasse10=RecoveredClasse10_0,
           DeathClasse1=DeathClasse1_0, DeathClasse2=DeathClasse2_0, DeathClasse3=DeathClasse3_0, DeathClasse4=DeathClasse4_0, DeathClasse5=DeathClasse5_0,
           DeathClasse6=DeathClasse6_0, DeathClasse7=DeathClasse7_0, DeathClasse8=DeathClasse8_0, DeathClasse9=DeathClasse9_0, DeathClasse10=DeathClasse10_0
)


#Structure

START<-0; FINISH<-56; STEP<-1
times <- seq(START, FINISH, by=STEP)
NUM_COHORTS<-10 
NUM_STATES<-7

#Résolution du système d’équations différentielles
Automatisation_Test <- 69

#Résolution du système d’équations différentielles
out   <- ode(y=yini, times=times ,func=mySIRS, parms=pars, method="euler")
Data_Mod <- data.frame(out)
#Data_Mod <- data.frame(Data_Mod[1],Data_Mod[Automatisation_Test])



##Observation des résultats de simulation
summary(out)
print(out)
#plot(out)


##==============================================================
## PART 2 : Importation Données et calibration
##==============================================================

#Détinition conditions de calibration
TT    <- seq(0, 56, 1)

#Importation des données 

#dat <- fread("Covid_FME_3006.txt",select=c(Automatisation_Test))

Data <- data.frame(
  time = TT,
  fread("Covid_FME_3006.txt",select=c(2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19, 20,
                                    21,22,23,24,25,26,27,28,29,30,31,32,33,34,35,36,37,38,39,40,
                                    41,42,43,44,45,46,47,48,49,50,51,52,53,54,55,56,57,58,59,60,
                                    61,62,63,64,65,66,67,68,69,70,71)))

#Observation des données de calibration dans un graphe

#plot(TT, Data[,"InfectedSeriousClasse8"], ylim = c(0, 100), pch = 16, col = "red",
#     main = "Evolution nb Serious 70 ans ", xlab = "time", ylab = "Nb people")


##=========================
## Fitted with Analytical Solution
##=========================

#Construction Fonction test des paramètres

solvemySIR <- function(pars){
  stocks <- c(SusceptibleClasse1=SusceptibleClasse1_0 , SusceptibleClasse2=SusceptibleClasse2_0,SusceptibleClasse3=SusceptibleClasse3_0,
              SusceptibleClasse4=SusceptibleClasse4_0, SusceptibleClasse5=SusceptibleClasse5_0,SusceptibleClasse6=SusceptibleClasse6_0,
              SusceptibleClasse7=SusceptibleClasse7_0, SusceptibleClasse8=SusceptibleClasse8_0,SusceptibleClasse9=SusceptibleClasse9_0,
              SusceptibleClasse10=SusceptibleClasse10_0,
              InfectedAsymptomaticClasse1=InfectedAsymptomaticClasse1_0 ,InfectedAsymptomaticClasse2=InfectedAsymptomaticClasse2_0 ,
              InfectedAsymptomaticClasse3=InfectedAsymptomaticClasse3_0 , InfectedAsymptomaticClasse4=InfectedAsymptomaticClasse4_0 ,
              InfectedAsymptomaticClasse5=InfectedAsymptomaticClasse5_0 ,InfectedAsymptomaticClasse6=InfectedAsymptomaticClasse6_0 ,
              InfectedAsymptomaticClasse7=InfectedAsymptomaticClasse7_0 , InfectedAsymptomaticClasse8=InfectedAsymptomaticClasse8_0 ,
              InfectedAsymptomaticClasse9=InfectedAsymptomaticClasse9_0 ,InfectedAsymptomaticClasse10=InfectedAsymptomaticClasse10_0 ,
              InfectedMildClasse1=InfectedMildClasse1_0,InfectedMildClasse2=InfectedMildClasse2_0,InfectedMildClasse3=InfectedMildClasse3_0,
              InfectedMildClasse4=InfectedMildClasse4_0, InfectedMildClasse5=InfectedMildClasse5_0, InfectedMildClasse6=InfectedMildClasse6_0,
              InfectedMildClasse7=InfectedMildClasse7_0, InfectedMildClasse8=InfectedMildClasse8_0, InfectedMildClasse9=InfectedMildClasse9_0,
              InfectedMildClasse10=InfectedMildClasse10_0,
              InfectedSeriousClasse1=InfectedSeriousClasse1_0,InfectedSeriousClasse2=InfectedSeriousClasse2_0,InfectedSeriousClasse3=InfectedSeriousClasse3_0,
              InfectedSeriousClasse4=InfectedSeriousClasse4_0, InfectedSeriousClasse5=InfectedSeriousClasse5_0, InfectedSeriousClasse6=InfectedSeriousClasse6_0,
              InfectedSeriousClasse7=InfectedSeriousClasse7_0, InfectedSeriousClasse8=InfectedSeriousClasse8_0, InfectedSeriousClasse9=InfectedSeriousClasse9_0,
              InfectedSeriousClasse10=InfectedSeriousClasse10_0,
              InfectedCriticalClasse1=InfectedCriticalClasse1_0, InfectedCriticalClasse2=InfectedCriticalClasse2_0, InfectedCriticalClasse3=InfectedCriticalClasse3_0,
              InfectedCriticalClasse4=InfectedCriticalClasse4_0, InfectedCriticalClasse5=InfectedCriticalClasse5_0, InfectedCriticalClasse6=InfectedCriticalClasse6_0,
              InfectedCriticalClasse7=InfectedCriticalClasse7_0, InfectedCriticalClasse8=InfectedCriticalClasse8_0, InfectedCriticalClasse9=InfectedCriticalClasse9_0,
              InfectedCriticalClasse10=InfectedCriticalClasse10_0,
              RecoveredClasse1=RecoveredClasse1_0, RecoveredClasse2=RecoveredClasse2_0, RecoveredClasse3=RecoveredClasse3_0, RecoveredClasse4=RecoveredClasse4_0,
              RecoveredClasse5=RecoveredClasse5_0,RecoveredClasse6=RecoveredClasse6_0,RecoveredClasse7=RecoveredClasse7_0,RecoveredClasse8=RecoveredClasse8_0,
              RecoveredClasse9=RecoveredClasse9_0,RecoveredClasse10=RecoveredClasse10_0,
              DeathClasse1=DeathClasse1_0, DeathClasse2=DeathClasse2_0, DeathClasse3=DeathClasse3_0, DeathClasse4=DeathClasse4_0, DeathClasse5=DeathClasse5_0,
              DeathClasse6=DeathClasse6_0, DeathClasse7=DeathClasse7_0, DeathClasse8=DeathClasse8_0, DeathClasse9=DeathClasse9_0, DeathClasse10=DeathClasse10_0
  )
  
  
  auxs <-  c(parameter2= unname(pars["parameter2"])
              
  )
                                
                        
  
  return(data.frame(ode(y=yini, times=times ,func=mySIRS, parms=pars, method="euler")))    
}

## FITTING algorithm 
getCost <- function(P) {
  out <- solvemySIR(P)
  cost<-return(modCost(obs=Data, model=out))# residuals
  return(cost)
}

## Bornes pour la calibration


lower <- pars*0.999999

upper <-pars*1.000001

Fit<-modFit(p=pars,f=getCost,lower=lower,upper=upper)

###Résultats calibration
#optPar<-c(Fit$par)
#optMod<-solvemySIR(optPar)
#summary(Fit)

#plot(residuals(Fit))

