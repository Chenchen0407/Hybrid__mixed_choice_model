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
  modelName       = "ML",
  modelDescr      = "Mixed logit model",
  indivID         = "ID",  
  nCores          = 14,
  outputDirectory = "output_ML5"
)

# ################################################################# #
#### LOAD DATA AND APPLY ANY TRANSFORMATIONS                     ####
# ################################################################# #

### Loading data from package
### if data is to be loaded from a file (e.g. called data.csv), 
database = read.csv("dataset.csv",header=TRUE)
### database = apollo_swissRouteChoiceData
### for data dictionary, use ?apollo_swissRouteChoiceData

# ################################################################# #
#### DEFINE MODEL PARAMETERS                                     ####
# ################################################################# #

### Vector of parameters, including any that are kept fixed in estimation
apollo_beta=c(asc_c      = 0,
              
              b_owning_new            = 0,
              b_owning_old            = 0,
              b_frequency_daily       = 0,
              b_distance_5km          = 0,
              b_recharging_cost       = 0,
              b_benefit1_s            = 0,
              b_benefit2_s            = 0,
              b_recharge_time         = 0,
              b_renew_cost_c          = 0,
              
              mu_asc_s    = 0,
              sigma_asc_s = 1,
              
              mu_growing_rate    = 0,
              sigma_growing_rate = 1,
              
              mu_purchase_price    = 0,
              sigma_purchase_price = 1,
              
              mu_access_time    = 0,
              sigma_access_time = 1)

### Vector with names (in quotes) of parameters to be kept fixed at their starting value in apollo_beta, use apollo_beta_fixed = c() if none
apollo_fixed = c("asc_c")

# ################################################################# #
#### DEFINE RANDOM COMPONENTS                                    ####
# ################################################################# #

### Set parameters for generating draws
apollo_draws = list(
  interDrawsType = "halton",
  interNDraws    = 500,
  interUnifDraws = c(),
  interNormDraws = c("draws_asc_s","draws_growing_rate","draws_purchase_price","draws_access_time"),
  intraDrawsType = "halton",
  intraNDraws    = 0,
  intraUnifDraws = c(),
  intraNormDraws = c()
)

### Create random parameters
apollo_randCoeff = function(apollo_beta, apollo_inputs){
  randcoeff = list()
  
  randcoeff[["asc_s"]] = mu_asc_s + sigma_asc_s * draws_asc_s
  randcoeff[["b_growing_rate"]] = mu_growing_rate + sigma_growing_rate * draws_growing_rate
  randcoeff[["b_purchase_price"]] = mu_purchase_price + sigma_purchase_price * draws_purchase_price
  randcoeff[["b_access_time"]] = mu_access_time + sigma_access_time * draws_access_time
  
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
  
  ### Function initialisation: do not change the following three commands
  ### Attach inputs and detach after function exit
  apollo_attach(apollo_beta, apollo_inputs)
  on.exit(apollo_detach(apollo_beta, apollo_inputs))
  
  ### Create list of probabilities P
  P = list()
  
  ### List of utilities: these must use the same names as in mnl_settings, order is irrelevant
  V = list()
  V[["swapping"]]  = asc_s  + b_owning_new * X_owning_new + b_owning_old * X_owning_old + b_frequency_daily * X_frequency_daily + b_distance_5km * X_distance_5km +  b_growing_rate * S_growing_rate + b_purchase_price * S_purchase_price + b_access_time * S_access_time + b_recharging_cost * S_swapping_cost + b_benefit1_s * S_benefit1 + b_benefit2_s * S_benefit2 + b_recharge_time * S_swapping_time
  V[["charging"]]  = asc_c +  b_growing_rate * C_growing_rate + b_purchase_price * C_purchase_price + b_access_time * C_access_time + b_recharging_cost * C_charging_cost + b_recharge_time * C_charging_time + b_renew_cost_c * C_renew_cost
  
  ### Define settings for BC model component
  mnl_settings = list(
    alternatives  = c(swapping=2, charging=1), 
    avail         = list(swapping=S_av, charging=C_av), 
    choiceVar     = choice,
    utilities     = V
  )
  
  ### Compute probabilities using MNL model
  P[["model"]] = apollo_mnl(mnl_settings, functionality)
  
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

model = apollo_estimate(apollo_beta, apollo_fixed,apollo_probabilities, apollo_inputs,apollo_control)

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
