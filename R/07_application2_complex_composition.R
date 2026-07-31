############################################################
## Application 2:
## Complex Composition Example – COVID-19 Ordinal Scale
##
## Ordinal Composite Endpoint W:
## 0: death
## 1: hospitalized on mechanical ventilation/ECMO
## 2: hospitalized on supplemental oxygen
## 3: hospitalized not on supplemental oxygen
## 4: not hospitalized with symptoms and limitation in activity
## 5: not hospitalized with symptoms but no limitation in activity
## 6: not hospitalized without symptoms or limitation in activity
############################################################

source_local <- function(file) {
  if (file.exists(file)) {
    source(file)
  } else if (file.exists(file.path("R", file))) {
    source(file.path("R", file))
  } else {
    stop("Cannot find file: ", file)
  }
}

source_local("03_functions_vco_pmf_power.R")

library(pacman)
p_load(
  dplyr,
  tidyr,
  readr,
  stringr,
  magrittr,
  rvinecopulib,
  VineCopula # for BiCopTau2Par()
)


############################################################
# 1. Marginal component prevalences from TREAT NOW
############################################################
preval_all <- c(
  hospital_day15      = 0.014005602,
  usual_activity_day15 = 0.151260504,
  symptoms_any_day15  = 0.288515406,
  oxygen_day15        = 0.008403361,
  vent_day15          = 0.002801120
)

############################################################
# 2. Dependence assumptions
############################################################

Gauss_param_HC_TN <- BiCopTau2Par(
  family = 1,
  tau = 0.6
)

Gauss_param_LC_TN <- BiCopTau2Par(
  family = 1,
  tau = 0.2
)

##########
# 1. Hospitalized (hospital_day15==1) and ventilation (vent_day15==1) [oxygen can be 0 or 1]
p1 = as.numeric(preval_all['hospital_day15'])
p2 = as.numeric(preval_all['vent_day15'])                

# assume low corr 

pmf1 = matrix(BiCop_PMF(p1,p2,parC= Gauss_param_LC_TN,  # BiCop_PMF is in 03_functions_vco_pmf_power.R
                        specific_method='Gaussian',print=T),nrow=1)
colnames(pmf1) = sort(paste0('X',expand.grid(rep(list(0:1),2)) %>% 
                               apply(1,function(x) paste0(x,collapse=''))))
pmf1

pmf1[,'X11']
# X11 
# 0.0002677709 

##########
# 2. Hospitalized (hospital_day15==1) and oxygen (oxygen_day15==1) [vent can only be 0 or 1]
p1 = as.numeric(preval_all['hospital_day15'])
p2 = as.numeric(preval_all['oxygen_day15'])                

pmf2 = matrix(BiCop_PMF(p1,p2,parC=Gauss_param_HC_TN,#find_parC(p1,p2,rho=corr_mat[1,4]),
                        specific_method='Gaussian'),nrow=1)
colnames(pmf2) = sort(paste0('X',expand.grid(rep(list(0:1),2)) %>% 
                               apply(1,function(x) paste0(x,collapse=''))))
pmf2

pmf2[,'X11']
# X11 
# 0.00418656

# 3. Hospitalized (hospital_day15==1) and not oxygen (oxygen_day15==0)
pmf2[,'X10']
# X10 
# 0.009819042


############################################################
# hospital---symptoms---activity
############################################################
## hospital,activity | symptom
SetupInfo = VCO_Initialize(nComp = 3,
                           prev = as.numeric(preval_all[c('hospital_day15',
                                                          'symptoms_any_day15',
                                                          'usual_activity_day15')]),
                           vine_type = 'D',
                           method_list=c('Gaussian','Indep'),
                           rho_list = NULL,
                           param_list = list(c(Gauss_param_LC_TN,Gauss_param_HC_TN),
                                             0)
)
### get multivariate binary joint distribution
output = gen_jointpmf(SetupInfo) # gen_jointpmf is in 03_functions_vco_pmf_power.R
list_joint_pmf = output$list_joint_pmf
jointpmf = output$jointpmf
jointpmf

# 4. Not hospitalized (hospital_day15==0) with symptoms (symptoms_any_day15==1) and limitation in activity (usual_activity_day15==1)
jointpmf['X011',]
# 0.1243809

# 5. Not hospitalized (hospital_day15==0) with symptoms (symptoms_any_day15==1) but with no limitation in activity (usual_activity_day15==0)
jointpmf['X010',]
# 0.1558024

# 6. Not hospitalized (hospital_day15==0) without symptoms  (symptoms_any_day15==0) nor limitation in activity (usual_activity_day15==0)
jointpmf['X000',]
# 0.6828152


############################################################
# Construct ordinal composite distribution
############################################################

VC_ordinal_dist <- c(
  0,                    # death
  pmf1["X11"],          # hospitalized + ventilation
  pmf2["X11"],          # hospitalized + oxygen
  pmf2["X10"],          # hospitalized without oxygen
  jointpmf["X011", ],   # symptoms + activity limitation
  jointpmf["X010", ],   # symptoms without activity limitation
  jointpmf["X000", ]    # no symptoms/activity limitation
)

VC_ordinal_dist[7] <- 1 - sum(VC_ordinal_dist[1:6])

VC_ordinal_dist

############################################################
# 6. Sample size calculation
############################################################

theta_R <- log(1.75)

get_SS <- function(target_power0, theta_R0, n0, pr0, A0) {
  
  calc_power <- wh_power(
    theta_R = theta_R0,
    n = n0,
    pr = pr0,
    A = A0
  )
  
  calc_power - target_power0
}


# 80% power
uniroot(
  get_SS,
  lower = 10,
  upper = 1000,
  theta_R0 = theta_R,
  pr0 = VC_ordinal_dist,
  A0 = 1,
  target_power0 = 0.8
)$root


# 90% power
uniroot(
  get_SS,
  lower = 10,
  upper = 1000,
  theta_R0 = theta_R,
  pr0 = VC_ordinal_dist,
  A0 = 1,
  target_power0 = 0.9
)$root
