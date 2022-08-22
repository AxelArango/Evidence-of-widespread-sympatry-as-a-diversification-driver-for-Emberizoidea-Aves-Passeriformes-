##masterscript
##biofeobears
library(phytools)
library(cladoRcpp)
library(rexpokit)
library(snow)
library(BioGeoBEARS)

###loading data into the workspaes
setwd("~/Documents/Sympatry on Emberizoidea_Masterscript/biogeobears")
np<-"~/Documents/Sympatry on Emberizoidea_Masterscript/biogeobears"
wd<-np("~/Documents/Sympatry on Emberizoidea_Masterscript/biogeobears")
extdata_dir = np(system.file("extdata", package="BioGeoBEARS"))
extdata_dir
list.files(extdata_dir)
geogfn ="~/Documents/Sympatry on Emberizoidea_Masterscript/Data/breeding_phylip_k8.txt"#geographic data in phylip format
trfn="~/Documents/Sympatry on Emberizoidea_Masterscript/Data/nodedtree.txt"#Emberizoidea phylogenetic tree with named nodes

# Look at the raw geography text file:
moref(geogfn)
moref(trfn)
moref(trfn)
tr = read.tree(trfn)
#plot(tr)
# Look at your geographic range data:
tipranges = getranges_from_LagrangePHYLIP(lgdata_fn=geogfn)
#tipranges

# Maximum range size observed:
max(rowSums(dfnums_to_numeric(tipranges@df)))
max_range_size=5
##
# Intitialize a default model (DEC model)
BioGeoBEARS_run_object = define_BioGeoBEARS_run()

# Give BioGeoBEARS the location of the phylogeny Newick file
BioGeoBEARS_run_object$trfn = trfn

# Give BioGeoBEARS the location of the geography text file
BioGeoBEARS_run_object$geogfn = geogfn

# Input the maximum range size
BioGeoBEARS_run_object$max_range_size = max_range_size

BioGeoBEARS_run_object$min_branchlength = 0.000001    # Min to treat tip as a direct ancestor (no speciation event)
BioGeoBEARS_run_object$include_null_range = F    # set to FALSE for e.g. DEC* model, DEC*+J, etc.

BioGeoBEARS_run_object$on_NaN_error = -1e50    # returns very low lnL if parameters produce NaN error (underflow check)
BioGeoBEARS_run_object$speedup = TRUE          # shorcuts to speed ML search; use FALSE if worried (e.g. >3 params)
BioGeoBEARS_run_object$use_optimx = "GenSA"    # if FALSE, use optim() instead of optimx()
BioGeoBEARS_run_object$num_cores_to_use = 4
BioGeoBEARS_run_object$force_sparse = FALSE 
BioGeoBEARS_run_object = readfiles_BioGeoBEARS_run(BioGeoBEARS_run_object)

BioGeoBEARS_run_object$return_condlikes_table = TRUE
BioGeoBEARS_run_object$calc_TTL_loglike_from_condlikes_table = TRUE
BioGeoBEARS_run_object$calc_ancprobs = TRUE    # get ancestral states from optim run

####prerun
BioGeoBEARS_run_object

BioGeoBEARS_run_object$BioGeoBEARS_model_object
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table

#check
check_BioGeoBEARS_run(BioGeoBEARS_run_object)
#dec
runslow = F
resfn = "Allrange_DEC_k8_unconstrained_v2.Rdata"
if (runslow)
{
  res = bears_optim_run(BioGeoBEARS_run_object)
  res    
  
  save(res, file=resfn)
  resDEC = res
} else {
  # Loads to "res"
  load(resfn)
  resDEC = res
}


##dec+j
BioGeoBEARS_run_object = readfiles_BioGeoBEARS_run(BioGeoBEARS_run_object)
BioGeoBEARS_run_object$return_condlikes_table = TRUE
BioGeoBEARS_run_object$calc_TTL_loglike_from_condlikes_table = TRUE
BioGeoBEARS_run_object$calc_ancprobs = TRUE    # get ancestral states from optim run

dstart = resDEC$outputs@params_table["d","est"]
estart = resDEC$outputs@params_table["e","est"]
jstart = 0.0001

# Input starting values for d, e
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["d","init"] = dstart
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["d","est"] = dstart
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["e","init"] = estart
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["e","est"] = estart

# Add j as a free parameter
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["j","type"] = "free"
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["j","init"] = jstart
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["j","est"] = jstart

check_BioGeoBEARS_run(BioGeoBEARS_run_object)

resfn = "all range_DEC+J_k8_unconstrained_v2.Rdata"
runslow = F
if (runslow)
{
  #sourceall("/Dropbox/_njm/__packages/BioGeoBEARS_setup/")
  
  res = bears_optim_run(BioGeoBEARS_run_object)
  res    
  
  save(res, file=resfn)
  
  resDECj = res
} else {
  # Loads to "res"
  load(resfn)
  resDECj = res
}
#####DIVA(like)

# Good default settings to get ancestral states
BioGeoBEARS_run_object$return_condlikes_table = TRUE
BioGeoBEARS_run_object$calc_TTL_loglike_from_condlikes_table = TRUE
BioGeoBEARS_run_object$calc_ancprobs = TRUE    # get ancestral states from optim run

# Set up DIVALIKE model
# Remove subset-sympatry
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["s","type"] = "fixed"
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["s","init"] = 0.0
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["s","est"] = 0.0

BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["ysv","type"] = "2-j"
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["ys","type"] = "ysv*1/2"
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["y","type"] = "ysv*1/2"
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["v","type"] = "ysv*1/2"

# Allow classic, widespread vicariance; all events equiprobable
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["mx01v","type"] = "fixed"
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["mx01v","init"] = 0.5
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["mx01v","est"] = 0.5

# No jump dispersal/founder-event speciation
# BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["j","type"] = "free"
# BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["j","init"] = 0.01
# BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["j","est"] = 0.01

check_BioGeoBEARS_run(BioGeoBEARS_run_object)

runslow = F
resfn = "All range_DIVALIKE_k8_unconstrained_v2.Rdata"
if (runslow)
{
  res = bears_optim_run(BioGeoBEARS_run_object)
  res    
  
  save(res, file=resfn)
  resDIVALIKE = res
} else {
  # Loads to "res"
  load(resfn)
  resDIVALIKE = res
}
####DIVA like +j
BioGeoBEARS_run_object = readfiles_BioGeoBEARS_run(BioGeoBEARS_run_object)
BioGeoBEARS_run_object$return_condlikes_table = TRUE
BioGeoBEARS_run_object$calc_TTL_loglike_from_condlikes_table = TRUE
BioGeoBEARS_run_object$calc_ancprobs = TRUE    # get ancestral states from optim run
dstart = resDIVALIKE$outputs@params_table["d","est"]
estart = resDIVALIKE$outputs@params_table["e","est"]
jstart = 0.0001

# Input starting values for d, e
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["d","init"] = dstart
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["d","est"] = dstart
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["e","init"] = estart
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["e","est"] = estart

# Remove subset-sympatry
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["s","type"] = "fixed"
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["s","init"] = 0.0
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["s","est"] = 0.0

BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["ysv","type"] = "2-j"
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["ys","type"] = "ysv*1/2"
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["y","type"] = "ysv*1/2"
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["v","type"] = "ysv*1/2"

# Allow classic, widespread vicariance; all events equiprobable
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["mx01v","type"] = "fixed"
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["mx01v","init"] = 0.5
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["mx01v","est"] = 0.5

# Add jump dispersal/founder-event speciation
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["j","type"] = "free"
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["j","init"] = jstart
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["j","est"] = jstart

# Under DIVALIKE+J, the max of "j" should be 2, not 3 (as is default in DEC+J)
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["j","min"] = 0.00001
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["j","max"] = 1.99999

check_BioGeoBEARS_run(BioGeoBEARS_run_object)

resfn = "All range_DIVALIKE+J_k8_unconstrained_v2.Rdata"
runslow = F
if (runslow)
{
  #sourceall("/Dropbox/_njm/__packages/BioGeoBEARS_setup/")
  
  res = bears_optim_run(BioGeoBEARS_run_object)
  res    
  
  save(res, file=resfn)
  
  resDIVALIKEj = res
} else {
  # Loads to "res"
  load(resfn)
  resDIVALIKEj = res
}

#BAYAREALIKE
BioGeoBEARS_run_object = define_BioGeoBEARS_run()
BioGeoBEARS_run_object$trfn = trfn
BioGeoBEARS_run_object$geogfn = geogfn
BioGeoBEARS_run_object$max_range_size = max_range_size
BioGeoBEARS_run_object$min_branchlength = 0.000001    # Min to treat tip as a direct ancestor (no speciation event)
BioGeoBEARS_run_object$include_null_range = TRUE    # set to FALSE for e.g. DEC* model, DEC*+J, etc.
# (For DEC* and other "*" models, please cite: Massana, Kathryn A.; Beaulieu, 
#  Jeremy M.; Matzke, Nicholas J.; O'Meara, Brian C. (2015). Non-null Effects of 
#  the Null Range in Biogeographic Models: Exploring Parameter Estimation in the 
#  DEC Model. bioRxiv,  http://biorxiv.org/content/early/2015/09/16/026914 )
# Also: search script on "include_null_range" for other places to change

# Set up a time-stratified analysis:
#BioGeoBEARS_run_object$timesfn = "timeperiods.txt"
#BioGeoBEARS_run_object$dispersal_multipliers_fn = "manual_dispersal_multipliers.txt"
#BioGeoBEARS_run_object$areas_allowed_fn = "areas_allowed.txt"
#BioGeoBEARS_run_object$areas_adjacency_fn = "areas_adjacency.txt"
#BioGeoBEARS_run_object$distsfn = "distances_matrix.txt"
# See notes on the distances model on PhyloWiki's BioGeoBEARS updates page.

# Speed options and multicore processing if desired
BioGeoBEARS_run_object$on_NaN_error = -1e50    # returns very low lnL if parameters produce NaN error (underflow check)
BioGeoBEARS_run_object$speedup = TRUE          # shorcuts to speed ML search; use FALSE if worried (e.g. >3 params)
BioGeoBEARS_run_object$use_optimx = "GenSA"    # if FALSE, use optim() instead of optimx()
BioGeoBEARS_run_object$num_cores_to_use = 4
BioGeoBEARS_run_object$force_sparse = FALSE    # force_sparse=TRUE causes pathology & isn't much faster at this scale

# This function loads the dispersal multiplier matrix etc. from the text files into the model object. Required for these to work!
# (It also runs some checks on these inputs for certain errors.)
BioGeoBEARS_run_object = readfiles_BioGeoBEARS_run(BioGeoBEARS_run_object)

# Divide the tree up by timeperiods/strata (uncomment this for stratified analysis)
#BioGeoBEARS_run_object = section_the_tree(inputs=BioGeoBEARS_run_object, make_master_table=TRUE, plot_pieces=FALSE)
# The stratified tree is described in this table:
#BioGeoBEARS_run_object$master_table

# Good default settings to get ancestral states
BioGeoBEARS_run_object$return_condlikes_table = TRUE
BioGeoBEARS_run_object$calc_TTL_loglike_from_condlikes_table = TRUE
BioGeoBEARS_run_object$calc_ancprobs = TRUE    # get ancestral states from optim run

# Set up BAYAREALIKE model
# No subset sympatry
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["s","type"] = "fixed"
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["s","init"] = 0.0
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["s","est"] = 0.0

# No vicariance
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["v","type"] = "fixed"
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["v","init"] = 0.0
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["v","est"] = 0.0

# No jump dispersal/founder-event speciation
# BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["j","type"] = "free"
# BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["j","init"] = 0.01
# BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["j","est"] = 0.01

# Adjust linkage between parameters
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["ysv","type"] = "1-j"
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["ys","type"] = "ysv*1/1"
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["y","type"] = "1-j"

# Only sympatric/range-copying (y) events allowed, and with 
# exact copying (both descendants always the same size as the ancestor)
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["mx01y","type"] = "fixed"
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["mx01y","init"] = 0.9999
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["mx01y","est"] = 0.9999

# Check the inputs
check_BioGeoBEARS_run(BioGeoBEARS_run_object)

runslow = F
resfn = "All range_BAYAREALIKE_k8_unconstrained_v2.Rdata"
if (runslow)
{
  res = bears_optim_run(BioGeoBEARS_run_object)
  res    
  
  save(res, file=resfn)
  resBAYAREALIKE = res
} else {
  # Loads to "res"
  load(resfn)
  resBAYAREALIKE = res
}
#BAYAREALIKE+j
BioGeoBEARS_run_object$on_NaN_error = -1e50    # returns very low lnL if parameters produce NaN error (underflow check)
BioGeoBEARS_run_object$speedup = TRUE          # shorcuts to speed ML search; use FALSE if worried (e.g. >3 params)
BioGeoBEARS_run_object$use_optimx = "GenSA"
BioGeoBEARS_run_object$num_cores_to_use = 4
BioGeoBEARS_run_object$force_sparse = FALSE    # force_sparse=TRUE causes pathology & isn't much faster at this scale
BioGeoBEARS_run_object = readfiles_BioGeoBEARS_run(BioGeoBEARS_run_object)
BioGeoBEARS_run_object$return_condlikes_table = TRUE
BioGeoBEARS_run_object$calc_TTL_loglike_from_condlikes_table = TRUE
BioGeoBEARS_run_object$calc_ancprobs = TRUE    # get ancestral states from optim run
dstart = resBAYAREALIKE$outputs@params_table["d","est"]
estart = resBAYAREALIKE$outputs@params_table["e","est"]
jstart = 0.0001
# Input starting values for d, e
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["d","init"] = dstart
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["d","est"] = dstart
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["e","init"] = estart
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["e","est"] = estart

# No subset sympatry
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["s","type"] = "fixed"
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["s","init"] = 0.0
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["s","est"] = 0.0

# No vicariance
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["v","type"] = "fixed"
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["v","init"] = 0.0
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["v","est"] = 0.0

# *DO* allow jump dispersal/founder-event speciation (set the starting value close to 0)
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["j","type"] = "free"
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["j","init"] = jstart
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["j","est"] = jstart

# Under BAYAREALIKE+J, the max of "j" should be 1, not 3 (as is default in DEC+J) or 2 (as in DIVALIKE+J)
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["j","max"] = 0.99999

# Adjust linkage between parameters
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["ysv","type"] = "1-j"
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["ys","type"] = "ysv*1/1"
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["y","type"] = "1-j"

# Only sympatric/range-copying (y) events allowed, and with 
# exact copying (both descendants always the same size as the ancestor)
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["mx01y","type"] = "fixed"
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["mx01y","init"] = 0.9999
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["mx01y","est"] = 0.9999

# NOTE (NJM, 2014-04): BAYAREALIKE+J seems to crash on some computers, usually Windows 
# machines. I can't replicate this on my Mac machines, but it is almost certainly
# just some precision under-run issue, when optim/optimx tries some parameter value 
# just below zero.  The "min" and "max" options on each parameter are supposed to
# prevent this, but apparently optim/optimx sometimes go slightly beyond 
# these limits.  Anyway, if you get a crash, try raising "min" and lowering "max" 
# slightly for each parameter:
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["d","min"] = 0.0000001
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["d","max"] = 4.9999999

BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["e","min"] = 0.0000001
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["e","max"] = 4.9999999

BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["j","min"] = 0.00001
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["j","max"] = 0.99999

check_BioGeoBEARS_run(BioGeoBEARS_run_object)

resfn = "All range_BAYAREALIKE+J_k8_unconstrained_v2.Rdata"
runslow = F
if (runslow)
{
  res = bears_optim_run(BioGeoBEARS_run_object)
  res    
  
  save(res, file=resfn)
  
  resBAYAREALIKEj = res
} else {
  # Loads to "res"
  load(resfn)
  resBAYAREALIKEj = res
}

resBAYAREALIKEj$optim_result$trace.mat


results_object = resBAYAREALIKEj
res$optim_result
#create tables
##test table
# Set up empty tables to hold the statistical results
restable = NULL
teststable = NULL

#######################################################
# Statistics -- DEC vs. DEC+J
#######################################################
# We have to extract the log-likelihood differently, depending on the 
# version of optim/optimx
LnL_2 = get_LnL_from_BioGeoBEARS_results_object(resDEC)
LnL_1 = get_LnL_from_BioGeoBEARS_results_object(resDECj)

numparams1 = 3
numparams2 = 2
stats = AICstats_2models(LnL_1, LnL_2, numparams1, numparams2)
stats

# DEC, null model for Likelihood Ratio Test (LRT)
res2 = extract_params_from_BioGeoBEARS_results_object(results_object=resDEC, returnwhat="table", addl_params=c("j"), paramsstr_digits=4)
# DEC+J, alternative model for Likelihood Ratio Test (LRT)
res1 = extract_params_from_BioGeoBEARS_results_object(results_object=resDECj, returnwhat="table", addl_params=c("j"), paramsstr_digits=4)

# The null hypothesis for a Likelihood Ratio Test (LRT) is that two models
# confer the same likelihood on the data. See: Brian O'Meara's webpage:
# http://www.brianomeara.info/tutorials/aic
# ...for an intro to LRT, AIC, and AICc

rbind(res2, res1)
tmp_tests = conditional_format_table(stats)

restable = rbind(restable, res2, res1)
teststable = rbind(teststable, tmp_tests)

#######################################################
# Statistics -- DIVALIKE vs. DIVALIKE+J
#######################################################
# We have to extract the log-likelihood differently, depending on the 
# version of optim/optimx
LnL_2 = get_LnL_from_BioGeoBEARS_results_object(resDIVALIKE)
LnL_1 = get_LnL_from_BioGeoBEARS_results_object(resDIVALIKEj)

numparams1 = 3
numparams2 = 2
stats = AICstats_2models(LnL_1, LnL_2, numparams1, numparams2)
stats

# DIVALIKE, null model for Likelihood Ratio Test (LRT)
res2 = extract_params_from_BioGeoBEARS_results_object(results_object=resDIVALIKE, returnwhat="table", addl_params=c("j"), paramsstr_digits=4)
# DIVALIKE+J, alternative model for Likelihood Ratio Test (LRT)
res1 = extract_params_from_BioGeoBEARS_results_object(results_object=resDIVALIKEj, returnwhat="table", addl_params=c("j"), paramsstr_digits=4)

rbind(res2, res1)
conditional_format_table(stats)

tmp_tests = conditional_format_table(stats)

restable = rbind(restable, res2, res1)
teststable = rbind(teststable, tmp_tests)

#######################################################
# Statistics -- BAYAREALIKE vs. BAYAREALIKE+J
#######################################################
# We have to extract the log-likelihood differently, depending on the 
# version of optim/optimx
LnL_2 = get_LnL_from_BioGeoBEARS_results_object(resBAYAREALIKE)
LnL_1 = get_LnL_from_BioGeoBEARS_results_object(resBAYAREALIKEj)

numparams1 = 3
numparams2 = 2
stats = AICstats_2models(LnL_1, LnL_2, numparams1, numparams2)
stats

# BAYAREALIKE, null model for Likelihood Ratio Test (LRT)
res2 = extract_params_from_BioGeoBEARS_results_object(results_object=resBAYAREALIKE, returnwhat="table", addl_params=c("j"), paramsstr_digits=4)
# BAYAREALIKE+J, alternative model for Likelihood Ratio Test (LRT)
res1 = extract_params_from_BioGeoBEARS_results_object(results_object=resBAYAREALIKEj, returnwhat="table", addl_params=c("j"), paramsstr_digits=4)

rbind(res2, res1)
conditional_format_table(stats)

tmp_tests = conditional_format_table(stats)

restable = rbind(restable, res2, res1)
teststable = rbind(teststable, tmp_tests)

#########################################################################
# ASSEMBLE RESULTS TABLES: DEC, DEC+J, DIVALIKE, DIVALIKE+J, BAYAREALIKE, BAYAREALIKE+J
#########################################################################
teststable$alt = c("DEC+J", "DIVALIKE+J", "BAYAREALIKE+J")
teststable$null = c("DEC", "DIVALIKE", "BAYAREALIKE")
row.names(restable) = c("DEC", "DEC+J", "DIVALIKE", "DIVALIKE+J", "BAYAREALIKE", "BAYAREALIKE+J")
restable = put_jcol_after_ecol(restable)
restable

# Look at the results!!
restable
teststable

#######################################################
# Save the results tables for later -- check for e.g.
# convergence issues
#######################################################

# Loads to "restable"
save(restable, file="restable_v1.Rdata")
load(file="restable_v1.Rdata")

# Loads to "teststable"
save(teststable, file="teststable_v1.Rdata")
load(file="teststable_v1.Rdata")

# Also save to text files
write.table(restable, file="restable.txt", quote=FALSE, sep="\t")
write.table(unlist_df(teststable), file="teststable.txt", quote=FALSE, sep="\t")

#######################################################
# Model weights of all six models
#######################################################
restable2 = restable

# With AICs:
AICtable = calc_AIC_column(LnL_vals=restable$LnL, nparam_vals=restable$numparams)
restable = cbind(restable, AICtable)
restable_AIC_rellike = AkaikeWeights_on_summary_table(restable=restable, colname_to_use="AIC")
restable_AIC_rellike = put_jcol_after_ecol(restable_AIC_rellike)
restable_AIC_rellike

# With AICcs -- factors in sample size
samplesize = length(tr$tip.label)
AICtable = calc_AICc_column(LnL_vals=restable$LnL, nparam_vals=restable$numparams, samplesize=samplesize)
restable2 = cbind(restable2, AICtable)
restable_AICc_rellike = AkaikeWeights_on_summary_table(restable=restable2, colname_to_use="AICc")
restable_AICc_rellike = put_jcol_after_ecol(restable_AICc_rellike)
restable_AICc_rellike

# Also save to text files
write.table(restable_AIC_rellike, file="restable_AIC_rellike.txt", quote=FALSE, sep="\t")
write.table(restable_AICc_rellike, file="restable_AICc_rellike.txt", quote=FALSE, sep="\t")

# Save with nice conditional formatting
write.table(conditional_format_table(restable_AIC_rellike), file="restable_AIC_rellike_formatted.txt", quote=FALSE, sep="\t")
write.table(conditional_format_table(restable_AICc_rellike), file="restable_AICc_rellike_formatted.txt", quote=FALSE, sep="\t")

############################
#BSM
############################
setwd("~/Documents/Sympatry on Emberizoidea_Masterscript/bsm")
#bestmodel (or the chosen one)
model_name = "BAYArea"
res = resBAYAREALIKE

#pdffn = paste0("Emberizoidea", model_name, "_v1.pdf")
#pdf(pdffn, width=6, height=6)

#analysis_titletxt = paste0(model_name, " on Emberizoidea")

# Setup
results_object = res
scriptdir = np(system.file("extdata/a_scripts", package="BioGeoBEARS"))

#plot States
#res2 = plot_BioGeoBEARS_results(results_object, analysis_titletxt, addl_params=list("j"), plotwhat="text", label.offset=0.45, tipcex=0.7, statecex=0.7, splitcex=0.6, titlecex=0.8, plotsplits=TRUE, cornercoords_loc=scriptdir, include_null_range=TRUE, tr=tr, tipranges=tipranges)

#plot Pie chart
#plot_BioGeoBEARS_results(results_object, analysis_titletxt, addl_params=list("j"), plotwhat="pie", label.offset=0.45, tipcex=0.7, statecex=0.7, splitcex=0.6, titlecex=0.8, plotsplits=TRUE, cornercoords_loc=scriptdir, include_null_range=TRUE, tr=tr, tipranges=tipranges)

#dev.off()  # Turn off PDF
#cmdstr = paste("open ", pdffn, sep="")
#system(cmdstr) # Plot it

#######################################################
# Stochastic mapping
#######################################################
clado_events_tables = NULL
ana_events_tables = NULL
lnum = 0

#######################################################
#######################################################
BSM_inputs_fn = "BSM_inputs_file.Rdata"
runInputsSlow = T
if (runInputsSlow)
{
  stochastic_mapping_inputs_list = get_inputs_for_stochastic_mapping(res=res)
  save(stochastic_mapping_inputs_list, file=BSM_inputs_fn)
} else {
  # Loads to "stochastic_mapping_inputs_list"
  load(BSM_inputs_fn)
} # END if (runInputsSlow)

# Check inputs (doesn't work the same on unconstr)
names(stochastic_mapping_inputs_list)
stochastic_mapping_inputs_list$phy2
stochastic_mapping_inputs_list$COO_weights_columnar
#stochastic_mapping_inputs_list$unconstr
set.seed(seed=as.numeric(Sys.time()))

runBSMslow = T
if (runBSMslow == TRUE)
{
  # Saves to: RES_clado_events_tables.Rdata
  # Saves to: RES_ana_events_tables.Rdata
  BSM_output = runBSM(res, stochastic_mapping_inputs_list=stochastic_mapping_inputs_list, maxnum_maps_to_try=1000, nummaps_goal=1000, maxtries_per_branch=40000, save_after_every_try=TRUE, savedir=getwd(), seedval=12345, wait_before_save=0.01)
  
  RES_clado_events_tables = BSM_output$RES_clado_events_tables
  RES_ana_events_tables = BSM_output$RES_ana_events_tables
} else {
  # Load previously saved...
  
  # Loads to: RES_clado_events_tables
  load(file="RES_clado_events_tables.Rdata")
  # Loads to: RES_ana_events_tables
  load(file="RES_ana_events_tables.Rdata")
  BSM_output = NULL
  BSM_output$RES_clado_events_tables = RES_clado_events_tables
  BSM_output$RES_ana_events_tables = RES_ana_events_tables
} # END if (runBSMslow == TRUE)

# Extract BSM output
clado_events_tables = BSM_output$RES_clado_events_tables
ana_events_tables = BSM_output$RES_ana_events_tables
#head(clado_events_tables[[1]])
#head(ana_events_tables[[1]])
#length(clado_events_tables)
#length(ana_events_tables)

include_null_range = TRUE
areanames = names(tipranges@df)
areas = areanames
max_range_size = 5


states_list_0based = rcpp_areas_list_to_states_list(areas=areas, maxareas=max_range_size, include_null_range=include_null_range)

colors_list_for_states = get_colors_for_states_list_0based(areanames=areanames, states_list_0based=states_list_0based, max_range_size=max_range_size, plot_null_range=TRUE)

############################################
# Setup for painting a single stochastic map
############################################
scriptdir = np(system.file("extdata/a_scripts", package="BioGeoBEARS"))
stratified = FALSE
clado_events_table = clado_events_tables[[1]]
ana_events_table = ana_events_tables[[1]]

cols_to_get = names(clado_events_table[,-ncol(clado_events_table)])
colnums = match(cols_to_get, names(ana_events_table))
 ana_events_table_cols_to_add = ana_events_table[,colnums]
 anagenetic_events_txt_below_node = rep("none", nrow(ana_events_table_cols_to_add))
 ana_events_table_cols_to_add = cbind(ana_events_table_cols_to_add, anagenetic_events_txt_below_node)
 rows_to_get_TF = ana_events_table_cols_to_add$node <= length(tr$tip.label)
 master_table_cladogenetic_events = rbind(ana_events_table_cols_to_add[rows_to_get_TF,], clado_events_table)

############################################
# Open a PDF fo single stochastic map
############################################
#pdffn = paste0(model_name, "_single_stochastic_map_n1.pdf")
#pdf(file=pdffn, width=6, height=6)

#Convert the BSM into a modified res object
master_table_cladogenetic_events = clado_events_tables[[1]]
resmod = stochastic_map_states_into_res(res=res, master_table_cladogenetic_events=master_table_cladogenetic_events, stratified=stratified)

plot_BioGeoBEARS_results(results_object=resmod, analysis_titletxt="Stochastic map", addl_params=list("j"), label.offset=0.5, plotwhat="text", cornercoords_loc=scriptdir, root.edge=TRUE, colors_list_for_states=colors_list_for_states, skiptree=FALSE, show.tip.label=TRUE)

 #Paint on the branch states
#paint_stochastic_map_branches(res=resmod, master_table_cladogenetic_events=master_table_cladogenetic_events, colors_list_for_states=colors_list_for_states, lwd=5, lty=par("lty"), root.edge=TRUE, stratified=stratified)

#plot_BioGeoBEARS_results(results_object=resmod, analysis_titletxt="Stochastic map", addl_params=list("j"), plotwhat="text", cornercoords_loc=scriptdir, root.edge=TRUE, colors_list_for_states=colors_list_for_states, skiptree=TRUE, show.tip.label=TRUE)

############################################
# Close PDF
############################################
#dev.off()
#cmdstr = paste("open ", pdffn, sep="")
#system(cmdstr)

#######################################################
# Summarize stochastic map tables
#######################################################
length(clado_events_tables)
length(ana_events_tables)

head(clado_events_tables[[1]][,-20])
tail(clado_events_tables[[1]][,-20])

head(ana_events_tables[[1]])
tail(ana_events_tables[[1]])

areanames = names(tipranges@df)
actual_names = areanames
actual_names

 #Get the dmat and times (if any)
dmat_times = get_dmat_times_from_res(res=res, numstates=NULL)
dmat_times

# Extract BSM output
clado_events_tables = BSM_output$RES_clado_events_tables
ana_events_tables = BSM_output$RES_ana_events_tables

# Simulate the source areas
BSMs_w_sourceAreas = simulate_source_areas_ana_clado(res, clado_events_tables, ana_events_tables, areanames)
clado_events_tables = BSMs_w_sourceAreas$clado_events_tables
ana_events_tables = BSMs_w_sourceAreas$ana_events_tables

# Count all anagenetic and cladogenetic events
counts_list = count_ana_clado_events(clado_events_tables, ana_events_tables, areanames, actual_names)

summary_counts_BSMs = counts_list$summary_counts_BSMs
print(conditional_format_table(summary_counts_BSMs))

# Histogram of event counts
hist_event_counts(counts_list, pdffn=paste0(model_name, "_histograms_of_event_counts_v2.pdf"))

#######################################################
# Print counts to files
#######################################################
tmpnames = names(counts_list)
cat("\n\nWriting tables* of counts to tab-delimited text files:\n(* = Tables have dimension=2 (rows and columns). Cubes (dimension 3) and lists (dimension 1) will not be printed to text files.) \n\n")
for (i in 1:length(tmpnames))
{
  cmdtxt = paste0("item = counts_list$", tmpnames[i])
  eval(parse(text=cmdtxt))
  
  # Skip cubes
  if (length(dim(item)) != 2)
  {
    next()
  }
  
  outfn = paste0(tmpnames[i], ".txt")
  if (length(item) == 0)
  {
    cat(outfn, " -- NOT written, *NO* events recorded of this type", sep="")
    cat("\n")
  } else {
    cat(outfn)
    cat("\n")
    write.table(conditional_format_table(item), file=outfn, quote=FALSE, sep="\t", col.names=TRUE, row.names=TRUE)
  } # END if (length(item) == 0)
} # END for (i in 1:length(tmpnames))
cat("...done.\n")

#######################################################
# Check that ML ancestral state/range probabilities and
# the mean of the BSMs approximately line up
#######################################################
#library(MultinomialCI)    # For 95% CIs on BSM counts
#check_ML_vs_BSM(res, clado_events_tables, model_name, tr=NULL, plot_each_node=FALSE, linreg_plot=TRUE, MultinomialCI=TRUE)

#########Transition rates measures
#calculate dispertion and extirpation rate
parprobx<-read.csv("~/Documents/Sympatry on Emberizoidea_Masterscript/Data/parprobx.csv",header=T)
familix<-unique(as.character(na.omit(parprobx$fam)))
event="e"
floxi<-data.frame()
for(k in 1:length(familix)){
  test3<-parprobx[which(parprobx$fam==familix[k]),]
  clex<-vector()
  for (i in 1:length(ana_events_tables)){
    anatest1<-ana_events_tables[[i]]
    anatest2<-anatest1[which(anatest1$event_type==event),]
    mergan<-merge(test3,anatest2,by="node",all=F)
    merx<-length(unique(mergan$node))
    if(event=="e"){clex[[i]]<-merx+1}else{clex[[i]]<-merx}
  }
  meanx<-mean(clex)
  sdx<-sd(clex)
  nevents<-length(test3$node)+1
  srate<-meanx/nevents
  ntips<-length(which(test3$node.type=="tip"))
  ratefam<-as.data.frame(cbind(ntips,nevents,meanx,sdx,srate,unique(test3$fam)));names(ratefam)<-c("ntips","nevents","mean","sd","rate","family")
  floxi<-rbind(floxi,ratefam)
}



extinctionrate<-floxi

event="d"
floxi<-data.frame()
for(k in 1:length(familix)){
  test3<-parprobx[which(parprobx$fam==familix[k]),]
  clex<-vector()
  for (i in 1:length(ana_events_tables)){
    anatest1<-ana_events_tables[[i]]
    anatest2<-anatest1[which(anatest1$event_type==event),]
    mergan<-merge(test3,anatest2,by="node",all=F)
    merx<-length(unique(mergan$node))
    if(event=="e"){clex[[i]]<-merx+1}else{clex[[i]]<-merx}
  }
  meanx<-mean(clex)
  sdx<-sd(clex)
  nevents<-length(test3$node)+1
  srate<-meanx/nevents
  ntips<-length(which(test3$node.type=="tip"))
  ratefam<-as.data.frame(cbind(ntips,nevents,meanx,sdx,srate,unique(test3$fam)));names(ratefam)<-c("ntips","nevents","mean","sd","rate","family")
  floxi<-rbind(floxi,ratefam)
}


dispertionrate<-floxi

write.csv(extinctionrate,"~/Documents/Sympatry on Emberizoidea_Masterscript/Data/Bextinction rate_bayarea.csv",row.names = F)
write.csv(dispertionrate,"~/Documents/Sympatry on Emberizoidea_Masterscript/Data/Bdispertion rate_bayarea.csv",row.names = F)
extinctionrate
dispertionrate
#Diversification measures
require(ape)
###load tree and taxonomy
mcctree<-read.tree("~/Documents/Sympatry on Emberizoidea_Masterscript/Data/Breeding-phylogeny_MCC.tre")
comperx<-read.csv("~/Documents/Sympatry on Emberizoidea_Masterscript/Data//BAMM-MMC_taxonomy.csv",header=T)
comperx2<-comperx[,-2]
#calculate Diversification Rate (DR)
DR_statistic <- function(x, return.mean = FALSE){
  rootnode <- length(x$tip.label) + 1
  sprates <- numeric(length(x$tip.label))
  for (i in 1:length(sprates)){
    node <- i
    index <- 1
    qx <- 0
    while (node != rootnode){
      el <- x$edge.length[x$edge[,2] == node]
      node <- x$edge[,1][x$edge[,2] == node]			
      qx <- qx + el* (1 / 2^(index-1))			
      index <- index + 1
    }
    sprates[i] <- 1/qx
  }
  if (return.mean){
    return(mean(sprates))		
  }else{
    names(sprates) <- x$tip.label
    return(sprates)
  }
}


DR<-DR_statistic(mcctree)
DRdf<-as.data.frame(DR)
drnames<-row.names(DRdf)
DRx<-data.frame(drnames,DRdf$DR)
colnames(DRx)<-c("sp","DR")




###load clads ClaDS results
require("RPANDA")
load("~/Documents/Sympatry on Emberizoidea_Masterscript/Diversification/ClaDS/clads_output.RData")
cladsx<-data.frame(drnames,CladsOutput$lambdatip_map)
colnames(cladsx)<-c("sp","ClaDS")
DivRat<-merge(DRx,cladsx,by="sp")
cor.test(DivRat$DR,DivRat$ClaDS)

#load BAMM results
library(BAMMtools)
library(phytools)
setwd("~/Documents/Sympatry on Emberizoidea_Masterscript/Diversification/bamm/")
mcmcout <- read.csv("BAMM_mcmc_emberizo_1.txt", header=T)
plot(mcmcout$logLik ~ mcmcout$generation)
burnstart <- floor(0.1 * nrow(mcmcout))
postburn <- mcmcout[burnstart:nrow(mcmcout), ]


embtree<-mcctree
postfile<-"BAMM_mcmc_emberizo_1.txt"
eventsx<-"BAMM_event_emberizo_1.txt"
###getting the  rate shifts
bfmat<-computeBayesFactors(postfile, expectedNumberOfShifts = 1, burnin = 0.25)
edata<-getEventData(embtree, eventsx, burnin = 0.25)
msc.set <- maximumShiftCredibility(edata, maximize='product')
msc.config <- subsetEventData(edata, index = msc.set$sampleindex)
###obtaning lambda per tip
bammlambda<-msc.config$meanTipLambda
bammlambdax<-data.frame(drnames,bammlambda)
colnames(bammlambdax)<-c("sp","BAMM")
#correlation between clads and BAMM
DivRat<-merge(DivRat,bammlambdax,by="sp")
cor.test(DivRat$DR,DivRat$BAMM)
cor.test(DivRat$BAMM,DivRat$ClaDS)

#bamm 2 using expected 10
#setwd("~/Desktop/mcmcs/mcm_post/")
#mcmcout <- read.csv("BAMM_mcmc_emberizo_test.txt", header=T)
#plot(mcmcout$logLik ~ mcmcout$generation)
#burnstart <- floor(0.1 * nrow(mcmcout))
#postburn <- mcmcout[burnstart:nrow(mcmcout), ]
#postfile<-"BAMM_mcmc_emberizo_test.txt"
#eventsx<-"BAMM_event_emberizo_test.txt"
#bfmat<-computeBayesFactors(postfile, expectedNumberOfShifts = 10, burnin = 0.25)
#edata<-getEventData(embtree, eventsx, burnin = 0.25)
#msc.set <- maximumShiftCredibility(edata, maximize='product')
#msc.config <- subsetEventData(edata, index = msc.set$sampleindex)

##Speciation rates mean by family
DivRatF<-merge(DivRat,comperx2,by.x="sp",by.y = "species")
Fdr<-tapply(DivRatF$DR,DivRatF$family,mean)
Fdrdf<-data.frame(row.names(Fdr),Fdr,row.names = NULL);colnames(Fdrdf)<-c("family","DR")
Fclads<-tapply(DivRatF$ClaDS,DivRatF$family,mean)
Fcladsdf<-data.frame(row.names(Fclads),Fclads,row.names = NULL);colnames(Fcladsdf)<-c("family","ClaDS")
Fbamm<-tapply(DivRatF$BAMM,DivRatF$family,mean)
Fbammdf<-data.frame(row.names(Fbamm),Fbamm,row.names = NULL);colnames(Fbammdf)<-c("family","BAMM")

Famrates<-merge(Fdrdf,Fcladsdf,by="family")
Famrates<-merge(Famrates,Fbammdf,by="family")
Famrates
#correlations between measures
cor.test(Famrates$BAMM,Famrates$ClaDS)
cor.test(Famrates$DR,Famrates$ClaDS)
cor.test(Famrates$DR,Famrates$BAMM)

########calculate method of moments
library(caper)
library(geiger)
setwd("~/Documents/Sympatry on Emberizoidea_Masterscript/Diversification/rates/")
missingx3<-read.csv("missing_families.csv",header=T)
families<-as.character(na.omit(unique(missingx3$families)))
msrates<-data.frame()
treez<-embtree
for (i in 1:length(families)){
  separt<-parprobx[which(parprobx$fam==families[i] & parprobx$node.type=="tip"),]
  separtcomp<-comparative.data(treez,separt,names.col=label)
  prtx<-prt(separtcomp$phy)
  timex<-prtx$time_bp[which(prtx$node.type=="root")]
  msx<-bd.ms(phy = separtcomp$phy,time= timex, missing =missingx3$missing[i] )
  msrate<-data.frame(cbind(families[i],msx));names(msrate)<-c("family","ms_rate")
  msrates<-rbind(msrates,msrate)
  i<-i+1
  
}

#comparations between all measurements
write.csv(msrates,"msrates.csv",row.names=F)
msrates<-read.csv("msrates.csv",header=T)
msratespurg<-msrates[which(msrates$ms_rate>=0),]
msratespurg
msratescomp<-merge(msratespurg,Famrates,by="family",all = F)
msratescomp

cor.test(DivRat$DR,DivRat$ClaDS)
cor.test(DivRat$DR,DivRat$BAMM)
cor.test(DivRat$BAMM,DivRat$ClaDS)
cor.test(Famrates$BAMM,Famrates$ClaDS)
cor.test(Famrates$DR,Famrates$ClaDS)
cor.test(Famrates$DR,Famrates$BAMM)
cor.test(msratescomp$ms_rate,msratescomp$DR)
cor.test(msratescomp$ms_rate,msratescomp$ClaDS)
cor.test(msratescomp$ms_rate,msratescomp$BAMM)
#####

####Creating a 1sp per family phylogenetic tree for PGLSs
library(geiger)
library(caper)
library(BioGeoBEARS)
library(ape)

spfam<-data.frame(DivRatF$sp,DivRatF$family);names(spfam)<-c("sp","family")
families<-as.character(unique(spfam$family))
embtree
speciestokeep<-data.frame()
for (i in 1:length(families)){
  tolook<-families[i]
  lookedfam<-DivRatF[which(DivRatF$family==tolook),]
  compfam<-comparative.data(embtree,lookedfam,names.col = "sp")
  prtphy<-prt(compfam$phy)
  numy<-log(prtphy$time_bp)
  numx<-sort(numy)
  numx2<-numx[1]
  newsp<-prtphy$label[which(log(prtphy$time_bp)==numx2)]
  spx<-data.frame(newsp,families[i]);names(spx)<-c("sp","family")
  speciestokeep<-rbind(speciestokeep,spx)}

speciestokeepx<-subset(speciestokeep, !duplicated(subset(speciestokeep, select=c(family))))
length(speciestokeepx$family)
#######PGLSs for diversification rates and transition rates
##still runs, but is better to run the dedicated script for that.
##########loads the transition rates
#spkee<-read.csv("~/Documents/masterscript/results/deprecated/pglsready_d.csv",header=T)
#spkee<-spkee[,1:2]
extinctionrate<-read.csv("~/Documents/Sympatry on Emberizoidea_Masterscript/Data/Extinction rate_bayarea.csv",header=T)
dispertionrate<-read.csv("~/Documents/Sympatry on Emberizoidea_Masterscript/Data/dispertion rate_bayarea.csv",header=T)
Famrates<-Fbammdf
#speciestokeepx<-spkee
exdiv<-merge(extinctionrate,Famrates,by="family")
disdiv<-merge(dispertionrate,Famrates,by="family")
#############pgls extinction
library(phytools)
library(caper)
###merges the extinction rate with the speciation rate
exdiv<-merge(extinctionrate,Famrates,by="family")
famratx<-data.frame(exdiv$family,exdiv$BAMM,exdiv$rate);names(famratx)<-c("family","BAMM","extinction_rate")
speciesmatrix<-merge(speciestokeepx,famratx,by="family")
speciesmatrix
#
####creates a comparative object so it can run the pgls
comprates<-comparative.data(embtree,speciesmatrix,names.col="sp",vcv.dim = T)
uniquephy<-comprates$phy
uniquephy$tip.label<-as.character(comprates$data$family)
#####BAMM~eTR pgls
ratepgls<-pgls(BAMM~extinction_rate,comprates,lambda = "ML")
summary(ratepgls)
#plot(log(comprates$data$extinction_rate),log(comprates$data$BAMM),ylab = "Diversification Rate (BAMM)",xlab = "extinction rate")abline(ratepgls,col="red")
#########pgls dispertion
disdiv<-merge(dispertionrate,Famrates,by="family")
famratx2<-data.frame(disdiv$family,disdiv$BAMM,disdiv$rate);names(famratx2)<-c("family","BAMM","dispertion_rate")
speciesmatrix2<-merge(speciestokeepx,famratx2,by="family")
comprates2<-comparative.data(embtree,speciesmatrix2,names.col="sp",vcv.dim = T)
##########BAMM~dTR pgls
ratepgls2<-pgls(BAMM~dispertion_rate,comprates2,lambda = "ML")
summary(ratepgls2)

