# ################################################################# #
#### LOAD LIBRARY AND DEFINE CORE SETTINGS                       ####
# ################################################################# #

### Clear memory
rm(list = ls())

### Load Apollo library
library(apollo)

### Initialise code
apollo_initialise()

### Set core controls
apollo_control = list(
  modelName       = "HCM",
  modelDescr      = "Hybrid choice model",
  indivID         = "ID",
  nCores          = 14, 
  outputDirectory = "output_HCM10"
)

# ################################################################# #
#### LOAD DATA AND APPLY ANY TRANSFORMATIONS                     ####
# ################################################################# #

### Loading data from package
### if data is to be loaded from a file (e.g. called data.csv), 
database = read.csv("dataset.csv",header=TRUE)
### for data dictionary, use ?apollo_drugChoiceData

# ################################################################# #
#### DEFINE MODEL PARAMETERS                                     ####
# ################################################################# #

### Vector of parameters, including any that are kept fixed in estimation
apollo_beta = c(asc_s      = 0,
                asc_c      = 0,
                
                b_owning_new            = 0,
                b_owning_old            = 0,
                b_frequency_daily       = 0,
                b_distance_5km          = 0,
                b_growing_rate          = 0,
                b_purchase_price        = 0,
                b_access_time           = 0,
                b_recharging_cost       = 0,
                b_benefit1_s            = 0,
                b_benefit2_s            = 0,
                b_recharge_time         = 0,
                b_renew_cost_c          = 0,
                
                lambda1             = 1, 
                lambda2             = 1,
                
                gamma_1             = 0,
                gamma_2             = 0,
                gamma_3             = 0,
                gamma2_1            = 0,
                gamma2_2            = 0,
                gamma2_3            = 0,
                gamma2_4            = 0,
                
                gamma_constant     = 0,
                gamma_male         = 0.5, 
                gamma_age          = 0.5, 
                gamma_edu          = 0.5, 
                gamma_job          = 0.5,
                gamma_income       = 0.5,

                gamma2_constant     = 0,
                gamma2_male         = 0.5, 
                gamma2_age          = 0.5, 
                gamma2_edu          = 0.5, 
                gamma2_job          = 0.5,
                gamma2_income       = 0.5,
                
                zeta_risk1       = 1,
                zeta_risk2       = 1,
                zeta_risk3       = 1,
                zeta_will1       = 1,
                zeta_will2       = 1,
                zeta_will3       = 1,
                zeta_will4       = 1,

                phi_risk1_1     =1, 
                phi_risk1_2     =1,
                phi_risk2_1     =1, 
                phi_risk2_2     =1,
                phi_risk3_1     =1, 
                phi_risk3_2     =1, 

                phi_will1_1     =1, 
                phi_will1_2     =1,
                phi_will2_1     =1, 
                phi_will2_2     =1,
                phi_will3_1     =1, 
                phi_will3_2     =1,
                phi_will4_1     =1, 
                phi_will4_2     =1,
                
                SD_eta1         =1,
                SD_eta2         =1,
                SD_risk1        =1,
                SD_risk2        =1,
                SD_risk3        =1,
                SD_will1        =1,
                SD_will2        =1,
                SD_will3        =1,
                SD_will4        =1
                )

### Vector with names (in quotes) of parameters to be kept fixed at their starting value in apollo_beta, use apollo_beta_fixed = c() if none
apollo_fixed = c("asc_c","zeta_risk1","zeta_will1","gamma_1","gamma2_1")

# ################################################################# #
#### DEFINE RANDOM COMPONENTS                                    ####
# ################################################################# #

### Set parameters for generating draws
apollo_draws = list(
  interDrawsType="sobol", 
  interNDraws=100,          
  interUnifDraws=c(),      
  interNormDraws=c("eta1","eta2","eta_risk1","eta_risk2","eta_risk3","eta_will1","eta_will2","eta_will3","eta_will4"), 
  
  intraDrawsType="",
  intraNDraws=0,          
  intraUnifDraws=c(),     
  intraNormDraws=c()      
)

### Create random parameters
apollo_randCoeff=function(apollo_beta, apollo_inputs){
  randcoeff = list()
  
  randcoeff[["RISK"]] = gamma_constant + gamma_male*X_male + gamma_age*X_age_35_or_more + gamma_edu*X_edu_college_or_more + gamma_job*X_job_office + gamma_income*X_income_9k_or_more + eta1 * SD_eta1
  randcoeff[["WILL"]] = gamma2_constant + gamma2_male*X_male + gamma2_age*X_age_35_or_more + gamma2_edu*X_edu_college_or_more + gamma2_job*X_job_office + gamma2_income*X_income_9k_or_more + eta2 * SD_eta2
  
  return(randcoeff)
}

# ################################################################# #
#### GROUP AND VALIDATE INPUTS                                   ####
# ################################################################# #

apollo_inputs = apollo_validateInputs()

# ################################################################# #
#### DEFINE MODEL AND LIKELIHOOD FUNCTION                        ####
# ################################################################# #

apollo_probabilities=function(apollo_beta, apollo_inputs, functionality="estimate"){

  ### Attach inputs and detach after function exit
  apollo_attach(apollo_beta, apollo_inputs)
  on.exit(apollo_detach(apollo_beta, apollo_inputs))
  
  ### Create list of probabilities P
  P = list()
  
  ### Likelihood of indicators
  ol_settings1 = list(outcomeOrdered = risk_risk1, 
                      V              = gamma_1 + zeta_risk1*RISK + eta_risk1 * SD_risk1, 
                      tau            = list(-phi_risk1_1-phi_risk1_2, -phi_risk1_1, phi_risk1_2, phi_risk1_1+phi_risk1_2), 
                      rows           = (S_av==1),
                      componentName  = "indic_risk1")
  ol_settings2 = list(outcomeOrdered = risk_risk2, 
                      V              = gamma_2 + zeta_risk2*RISK + eta_risk2 * SD_risk2, 
                      tau            = list(-phi_risk2_1-phi_risk2_2, -phi_risk2_1, phi_risk2_2, phi_risk2_1+phi_risk2_2), 
                      rows           = (S_av==1),
                      componentName  = "indic_risk2")
  ol_settings3 = list(outcomeOrdered = risk_risk3, 
                      V              = gamma_3 + zeta_risk3*RISK + eta_risk3 * SD_risk3, 
                      tau            = list(-phi_risk3_1-phi_risk3_2, -phi_risk3_1, phi_risk3_2, phi_risk3_1+phi_risk3_2), 
                      rows           = (S_av==1),
                      componentName  = "indic_risk3")
  ol_settings4 = list(outcomeOrdered = will_will1, 
                      V              = gamma2_1 + zeta_will1*WILL + eta_will1 * SD_will1, 
                      tau            = list(-phi_will1_1-phi_will1_2, -phi_will1_1, phi_will1_2, phi_will1_1+phi_will1_2), 
                      rows           = (S_av==1),
                      componentName  = "indic_will1")
  ol_settings5 = list(outcomeOrdered = will_will2, 
                      V              = gamma2_2 + zeta_will2*WILL + eta_will2 * SD_will2, 
                      tau            = list(-phi_will2_1-phi_will2_2, -phi_will2_1, phi_will2_2, phi_will2_1+phi_will2_2), 
                      rows           = (S_av==1),
                      componentName  = "indic_will2")
  ol_settings6 = list(outcomeOrdered = will_will3, 
                      V              = gamma2_3 + zeta_will3*WILL + eta_will3 * SD_will3, 
                      tau            = list(-phi_will3_1-phi_will3_2, -phi_will3_1, phi_will3_2, phi_will3_1+phi_will3_2), 
                      rows           = (S_av==1),
                      componentName  = "indic_will3")
  ol_settings7 = list(outcomeOrdered = will_will4, 
                      V              = gamma2_4 + zeta_will4*WILL + eta_will4 * SD_will4, 
                      tau            = list(-phi_will4_1-phi_will4_2, -phi_will4_1, phi_will4_2, phi_will4_1+phi_will4_2), 
                      rows           = (S_av==1),
                      componentName  = "indic_will4")

  P[["indic_risk1"]]     = apollo_ol(ol_settings1, functionality)
  P[["indic_risk2"]]     = apollo_ol(ol_settings2, functionality)  
  P[["indic_risk3"]]     = apollo_ol(ol_settings3, functionality)
  P[["indic_will1"]]     = apollo_ol(ol_settings4, functionality)
  P[["indic_will2"]]     = apollo_ol(ol_settings5, functionality)
  P[["indic_will3"]]     = apollo_ol(ol_settings6, functionality)
  P[["indic_will4"]]     = apollo_ol(ol_settings7, functionality)
  
  
  ### Likelihood of choices
  ### List of utilities: these must use the same names as in mnl_settings, order is irrelevant
  V = list()
  V[["swapping"]]  = asc_s  + b_owning_new * X_owning_new + b_owning_old * X_owning_old + b_frequency_daily * X_frequency_daily + b_distance_5km * X_distance_5km +  b_growing_rate * S_growing_rate + b_purchase_price * S_purchase_price + b_access_time * S_access_time + b_recharging_cost * S_swapping_cost + b_benefit1_s * S_benefit1 + b_benefit2_s * S_benefit2 + b_recharge_time * S_swapping_time + lambda1*RISK + lambda2*WILL
  V[["charging"]]  = asc_c +  b_growing_rate * C_growing_rate + b_purchase_price * C_purchase_price + b_access_time * C_access_time + b_recharging_cost * C_charging_cost + b_recharge_time * C_charging_time + b_renew_cost_c * C_renew_cost
  
  ### Define settings for BC model component
  mnl_settings = list(
    alternatives  = c(swapping=2, charging=1), 
    avail         = list(swapping=S_av, charging=C_av), 
    choiceVar     = choice,
    utilities     = V
  )
  
  ### Compute probabilities for MNL model component
  P[["choice"]] = apollo_mnl(mnl_settings, functionality)
  
  ### Likelihood of the whole model
  P = apollo_combineModels(P, apollo_inputs, functionality)
  
  ### Take product across observation for same individual
  P = apollo_panelProd(P, apollo_inputs, functionality)
  
  ### Average across inter-individual draws
  P = apollo_avgInterDraws(P, apollo_inputs, functionality)
  
  ### Prepare and return outputs of function
  P = apollo_prepareProb(P, apollo_inputs, functionality)
  return(P)
}

# ################################################################# #
#### MODEL ESTIMATION                                            ####
# ################################################################# #

### Optional: calculate LL before model estimation
apollo_llCalc(apollo_beta, apollo_probabilities, apollo_inputs)

### Estimate model
model <- apollo_estimate(apollo_beta, apollo_fixed, apollo_probabilities, apollo_inputs,estimate_settings = list(scaleAfterConvergence = FALSE))
# ################################################################# #
#### MODEL OUTPUTS                                               ####
# ################################################################# #

# ----------------------------------------------------------------- #
#---- FORMATTED OUTPUT (TO SCREEN)                               ----
# ----------------------------------------------------------------- #

apollo_modelOutput(model)

# ----------------------------------------------------------------- #
#---- FORMATTED OUTPUT (TO FILE, using model name)               ----
# ----------------------------------------------------------------- #

apollo_saveOutput(model)

# ################################################################# #
##### POST-PROCESSING                                            ####
# ################################################################# #

### Print outputs of additional diagnostics to new output file (remember to close file writing when complete)
apollo_sink()

# ----------------------------------------------------------------- #
#---- MODEL PREDICTIONS                                          ----
# ----------------------------------------------------------------- #

forecast <- apollo_prediction(model, apollo_probabilities, apollo_inputs,
                              prediction_settings=list(modelComponent="indic_quality"))

# ----------------------------------------------------------------- #
#---- CONDITIONALS AND UNCONDITIONALS                            ----
# ----------------------------------------------------------------- #

conditionals <- apollo_conditionals(model,apollo_probabilities,apollo_inputs)

summary(conditionals)

unconditionals <- apollo_unconditionals(model,apollo_probabilities,apollo_inputs)

mean(unconditionals[[1]])
sd(unconditionals[[1]])

# ----------------------------------------------------------------- #
#---- switch off writing to file                                 ----
# ----------------------------------------------------------------- #

apollo_sink()