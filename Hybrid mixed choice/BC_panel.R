# ################################################################# #
#### LOAD LIBRARY AND DEFINE CORE SETTINGS                       ####
# ################################################################# #

### Clear memory
rm(list = ls())

### Load Apollo library
library(apollo)
if(!require(car)) {
  install.packages("car")
}
library(car)

### Initialise code
apollo_initialise()

### Set core controls
apollo_control = list(
  modelName       = "BC_Panel_effects",
  modelDescr      = "Binary choice model on mode choice SP data using effects coding",
  indivID         = "ID",
  outputDirectory = "output_BC5"
)

# ################################################################# #
#### LOAD DATA AND APPLY ANY TRANSFORMATIONS                     ####
# ################################################################# #

### Loading data from package
### if data is to be loaded from a file (e.g. called data.csv), 
### the code would be: 
database = read.csv("dataset.csv",header=TRUE)
# 检查多重共线性
model <- lm(choice ~ X_future_purchase_s + X_frequency_daily + X_distance_5km + S_growing_rate + S_purchase_price + S_access_time + S_swapping_cost + S_benefit1 + S_benefit2 + S_benefit3 + S_swapping_time, data = database)
summary(model)

# ################################################################# #
#### DEFINE MODEL PARAMETERS                                     ####
# ################################################################# #

### Vector of parameters, including any that are kept fixed in estimation
apollo_beta=c(asc_c      = 0,
              asc_s      = 0,
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
              b_renew_cost_c          = 0)

### Vector with names (in quotes) of parameters to be kept fixed at their starting value in apollo_beta, use apollo_beta_fixed = c() if none
apollo_fixed = c("asc_c")

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
  
  ### Prepare and return outputs of function
  P = apollo_prepareProb(P, apollo_inputs, functionality)
  return(P)
}

# ################################################################# #
#### MODEL ESTIMATION                                            ####
# ################################################################# #

model = apollo_estimate(apollo_beta, apollo_fixed, 
                        apollo_probabilities, apollo_inputs)

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

### calculate value and standard error for base of effects coded parameter
apollo_deltaMethod(model,deltaMethod_settings = list(expression=c(b_no_frills="-b_wifi-b_food")))

# ----------------------------------------------------------------- #
#---- switch off writing to file                                 ----
# ----------------------------------------------------------------- #

apollo_sink()