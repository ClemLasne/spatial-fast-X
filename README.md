# spatial-fast-X
Codes for 'The relative contributions of the X chromosome and autosomes to local adaptation'

```ruby
# Environmental transition : abrupt
# dominance model: reversal
# deterministic model

########################### PARAMETERS ############################
# HABITATS
H = 100
# Number of generations per run
G = 10000
# sex-specific migration coefficients between two adjacent patches
mm = 0.1
mf = 0.1
# sex-specific selection coefficient against the deleterious allele
sm = 0.05
sf = 0.05
# dominance coefficient of the deleterious allele
h = 0.5

############ AUTOSOMAL AND X-LINKED GENOTYPE FREQUENCIES ############

###### AUTOSOMAL LOCUS ######
# vectors to hold genotype frequencies in FEMALES
fAA_a = rep(0,H)
fAB_a = rep(0,H)
fBB_a = rep(0,H)
# vectors to hold genotype frequencies in MALES
gAA_a = rep(0,H)
gAB_a = rep(0,H)
gBB_a = rep(0,H)

# Initial genotype frequencies at autosomal locus in FEMALES and MALES
for (i in 1:(H/2) ){
  fAA_a[i] = 0.25
  fAB_a[i] = 0.5
  fBB_a[i] = 0.25
  gAA_a[i] = 0.25
  gAB_a[i] = 0.5
  gBB_a[i] = 0.25
}

for (i in (H/2+1):H){
  fAA_a[i] = 0.25
  fAB_a[i] = 0.5
  fBB_a[i] = 0.25
  gAA_a[i] = 0.25
  gAB_a[i] = 0.5
  gBB_a[i] = 0.25
}

###### X-LINKED LOCUS ######
# vectors to hold genotype frequencies in FEMALES
fAA_x = rep(0,H)
fAB_x = rep(0,H)
fBB_x = rep(0,H)
# vectors to hold genotype frequencies in MALES
gA_x = rep(0,H)
gB_x = rep(0,H)

# Initial genotype frequencies at X-linked locus in FEMALES and MALES
for (i in 1:(H/2) ){
  fAA_x[i] = 0.25
  fAB_x[i] = 0.5
  fBB_x[i] = 0.25
  gA_x[i] = 0.5
  gB_x[i] = 0.5
}

for (i in (H/2+1):H){
  fAA_x[i] = 0.25
  fAB_x[i] = 0.5
  fBB_x[i] = 0.25
  gA_x[i] = 0.5
  gB_x[i] = 0.5
}

################################
######### 1. MIGRATION #########
################################

# vectors to hold genotype frequencies AFTER SEX-SPECIFIC MIGRATION

###### AUTOSOMAL LOCUS #######
# FEMALES
fAA_a_mig = rep(0,H)
fAB_a_mig = rep(0,H)
fBB_a_mig = rep(0,H)
# MALES
gAA_a_mig = rep(0,H)
gAB_a_mig = rep(0,H)
gBB_a_mig = rep(0,H)

###### X-LINKED LOCUS #######
# FEMALES
fAA_x_mig = rep(0,H)
fAB_x_mig = rep(0,H)
fBB_x_mig = rep(0,H)
# MALES
gA_x_mig = rep(0,H)
gB_x_mig = rep(0,H)


################################
######### 2. SELECTION #########
################################

# vectors to hold genotype frequencies AFTER SEX-SPECIFIC SELECTION

###### AUTOSOMAL LOCUS #######
# FEMALES
fAA_a_sel = rep(0,H)
fAB_a_sel = rep(0,H)
fBB_a_sel = rep(0,H)
# MALES
gAA_a_sel = rep(0,H)
gAB_a_sel = rep(0,H)
gBB_a_sel = rep(0,H)

###### X-LINKED LOCUS #######
# FEMALES
fAA_x_sel = rep(0,H)
fAB_x_sel = rep(0,H)
fBB_x_sel = rep(0,H)
# MALES
gA_x_sel = rep(0,H)
gB_x_sel = rep(0,H)

#### FITNESS ####

###### AUTOSOMAL LOCUS #######
# average population fitness per patch in FEMALES 
Wavg_a = rep(0,H)
# average population fitness per patch in MALES 
Vavg_a = rep(0,H)

# vectors for genotypes fitness
# FEMALES
WAA_a = rep(0,H)
WAB_a = rep(0,H)
WBB_a = rep(0,H)
# MALES
VAA_a = rep(0,H)
VAB_a = rep(0,H)
VBB_a = rep(0,H)

# fitness loops for the autosomal locus in the context of DOMINANCE REVERSAL
for (i in 1:(H/2) ) {
  # FEMALES
  WAA_a[i] = 1-sf
  WAB_a[i] = 1-h*sf
  WBB_a[i] = 1
  # MALES  
  VAA_a[i] = 1-sm
  VAB_a[i] = 1-h*sm
  VBB_a[i] = 1
}

for (i in (H/2+1):H){
  # FEMALES
  WAA_a[i] = 1
  WAB_a[i] = 1-h*sf
  WBB_a[i] = 1-sf
  # MALES  
  VAA_a[i] = 1
  VAB_a[i] = 1-h*sm
  VBB_a[i] = 1-sm
}

###### X-LINKED LOCUS #######
# average population fitness per patch at the autosomal locus in FEMALES 
Wavg_x = rep(0,H)
# average population fitness per patch at the autosomal locus in MALES 
Vavg_x = rep(0,H)

# vectors for genotypes fitness
# FEMALES
WAA_x = rep(0,H)
WAB_x = rep(0,H)
WBB_x = rep(0,H)
# MALES
VA_x = rep(0,H)
VB_x = rep(0,H)

# fitness loops for the X-linked locus in the context of DOMINANCE REVERSAL
for (i in 1:(H/2)){
  # FEMALES
  WAA_x[i] = 1-sf
  WAB_x[i] = 1-h*sf
  WBB_x[i] = 1
  # MALES
  VA_x[i] = 1-sm
  VB_x[i] = 1
}

for (i in (H/2+1):H){
  # FEMALES
  WAA_x[i] = 1
  WAB_x[i] = 1-h*sf
  WBB_x[i] = 1-sf
  # MALES
  VA_x[i] = 1
  VB_x[i] = 1-sm
}


##########################################
######### CLINE SLOPE EXTRACTION #########
##########################################
slope_Auto=rep(0,H)
slope_X=rep(0,H)

########################################################################################
####################################### LOOP ###########################################
########################################################################################

### START LOOP ###
for(j in 1:G){
  
  ### 1. SEX-SPECIFIC MIGRATION
  # 1.1 genotype frequencies after migration in patches 2 to H-1
  
  for(i in 2:(H-1)){
    
    ### AUTOSOMAL ###
    
    # FEMALES
    fAA_a_mig[i] = fAA_a[i]*(1-mf) + (mf/2)*(fAA_a[i-1] + fAA_a[i+1])
    fAB_a_mig[i] = fAB_a[i]*(1-mf) + (mf/2)*(fAB_a[i-1] + fAB_a[i+1])
    fBB_a_mig[i] = fBB_a[i]*(1-mf) + (mf/2)*(fBB_a[i-1] + fBB_a[i+1])
    # MALES
    gAA_a_mig[i] = gAA_a[i]*(1-mm) + (gAA_a[i-1] + gAA_a[i+1])*(mm/2)
    gAB_a_mig[i] = gAB_a[i]*(1-mm) + (gAB_a[i-1] + gAB_a[i+1])*(mm/2)
    gBB_a_mig[i] = gBB_a[i]*(1-mm) + (gBB_a[i-1] + gBB_a[i+1])*(mm/2)
    
    ###### X-LINKED LOCUS #######
    
    # FEMALES
    fAA_x_mig[i] = fAA_x[i]*(1-mf) + (fAA_x[i-1] + fAA_x[i+1])*(mf/2)
    fAB_x_mig[i] = fAB_x[i]*(1-mf) + (fAB_x[i-1] + fAB_x[i+1])*(mf/2) 
    fBB_x_mig[i] = fBB_x[i]*(1-mf) + (fBB_x[i-1] + fBB_x[i+1])*(mf/2)
    # MALES
    gA_x_mig[i] = gA_x[i]*(1-mm) + (gA_x[i-1] + gA_x[i+1])*(mm/2)
    gB_x_mig[i] = gB_x[i]*(1-mm) + (gB_x[i-1] + gB_x[i+1])*(mm/2) 
  }
  
  # 1.2 genotype frequencies after migration in patch 1
  
  ### AUTOSOMAL LOCUS ###
  
  # FEMALES
  fAA_a_mig[1] = fAA_a[2]*(mf/2) + fAA_a[1]*(1-(mf/2))
  fAB_a_mig[1] = fAB_a[2]*(mf/2) + fAB_a[1]*(1-(mf/2))
  fBB_a_mig[1] = fBB_a[2]*(mf/2) + fBB_a[1]*(1-(mf/2))
  # MALES
  gAA_a_mig[1] = gAA_a[2]*(mm/2) + gAA_a[1]*(1-(mm/2))
  gAB_a_mig[1] = gAB_a[2]*(mm/2) + gAB_a[1]*(1-(mm/2))
  gBB_a_mig[1] = gBB_a[2]*(mm/2) + gBB_a[1]*(1-(mm/2))
  
  ###### X-LINKED LOCUS #######
  
  # FEMALES
  fAA_x_mig[1] = fAA_x[2]*(mf/2) + fAA_x[1]*(1-(mf/2))
  fAB_x_mig[1] = fAB_x[2]*(mf/2) + fAB_x[1]*(1-(mf/2))
  fBB_x_mig[1] = fBB_x[2]*(mf/2) + fBB_x[1]*(1-(mf/2))
  # MALES
  gA_x_mig[1] = gA_x[2]*(mm/2) + gA_x[1]*(1-(mm/2))
  gB_x_mig[1] = gB_x[2]*(mm/2) + gB_x[1]*(1-(mm/2))
  
  # 1.3 genotype frequencies after migration in patch H
  
  ### AUTOSOMAL LOCUS ###
  
  # FEMALES
  fAA_a_mig[H] = fAA_a[H-1]*(mf/2) + fAA_a[H]*(1-(mf/2))
  fAB_a_mig[H] = fAB_a[H-1]*(mf/2) + fAB_a[H]*(1-(mf/2))
  fBB_a_mig[H] = fBB_a[H-1]*(mf/2) + fBB_a[H]*(1-(mf/2))
  # MALES
  gAA_a_mig[H] = gAA_a[H-1]*(mm/2) + gAA_a[H]*(1-(mm/2))
  gAB_a_mig[H] = gAB_a[H-1]*(mm/2) + gAB_a[H]*(1-(mm/2))
  gBB_a_mig[H] = gBB_a[H-1]*(mm/2) + gBB_a[H]*(1-(mm/2)) 
  
  ###### X-LINKED LOCUS #######
  
  # FEMALES
  fAA_x_mig[H] = fAA_x[H-1]*(mf/2) + fAA_x[H]*(1-(mf/2))
  fAB_x_mig[H] = fAB_x[H-1]*(mf/2) + fAB_x[H]*(1-(mf/2))
  fBB_x_mig[H] = fBB_x[H-1]*(mf/2) + fBB_x[H]*(1-(mf/2)) 
  # MALES
  gA_x_mig[H] = gA_x[H-1]*(mm/2) + gA_x[H]*(1-(mm/2))
  gB_x_mig[H] = gB_x[H-1]*(mm/2) + gB_x[H]*(1-(mm/2))
  
  
  ### 2. SEX-SPECIFIC SELECTION
  # 2.1. average population fitness
  
  for(i in 1:H){
    
    ### AUTOSOMAL LOCUS ###
    # FEMALE POPULATION
    Wavg_a[i] = WAA_a[i]*fAA_a_mig[i] + WAB_a[i]*fAB_a_mig[i] + WBB_a[i]*fBB_a_mig[i]
    # MALE POPULATION
    Vavg_a[i] = VAA_a[i]*gAA_a_mig[i] + VAB_a[i]*gAB_a_mig[i] + VBB_a[i]*gBB_a_mig[i]
    
    ### X-LINKED LOCUS ###
    # FEMALE POPULATION
    Wavg_x[i] = WAA_x[i]*fAA_x_mig[i] + WAB_x[i]*fAB_x_mig[i] + WBB_x[i]*fBB_x_mig[i]
    # MALE POPULATION
    Vavg_x[i] = VA_x[i]*gA_x_mig[i] + VB_x[i]*gB_x_mig[i]
  }
  
  # 2.2. genotype frequencies after selection
  for (i in 1:H){
    
    ### AUTOSOMAL LOCUS ###
    fAA_a_sel[i] = (WAA_a[i]*fAA_a_mig[i]) / Wavg_a[i]
    fAB_a_sel[i] = (WAB_a[i]*fAB_a_mig[i]) / Wavg_a[i]
    fBB_a_sel[i] = (WBB_a[i]*fBB_a_mig[i]) / Wavg_a[i]
    # MALES
    gAA_a_sel[i] = (VAA_a[i]*gAA_a_mig[i]) / Vavg_a[i]
    gAB_a_sel[i] = (VAB_a[i]*gAB_a_mig[i]) / Vavg_a[i]
    gBB_a_sel[i] = (VBB_a[i]*gBB_a_mig[i]) / Vavg_a[i]  
    
    ### X-LINKED LOCUS ###
    # FEMALES
    fAA_x_sel[i] = (WAA_x[i]*fAA_x_mig[i]) / Wavg_x[i]
    fAB_x_sel[i] = (WAB_x[i]*fAB_x_mig[i]) / Wavg_x[i]
    fBB_x_sel[i] = (WBB_x[i]*fBB_x_mig[i]) / Wavg_x[i] 
    # MALES
    gA_x_sel[i] = (VA_x[i]*gA_x_mig[i]) / Vavg_x[i]
    gB_x_sel[i] = (VB_x[i]*gB_x_mig[i]) / Vavg_x[i]
  }
  
  #3. Random mating (no drift)
  
  for(i in 1:H) {
    
    ### AUTOSOMAL LOCUS ###
    # FEMALES
    fAA_a[i] = fAA_a_sel[i]*gAA_a_sel[i] + fAA_a_sel[i]*gAB_a_sel[i]/2 + fAB_a_sel[i]*gAA_a_sel[i]/2 + fAB_a_sel[i]*gAB_a_sel[i]/4
    fBB_a[i] = fBB_a_sel[i]*gBB_a_sel[i] + fBB_a_sel[i]*gAB_a_sel[i]/2 + fAB_a_sel[i]*gBB_a_sel[i]/2 + fAB_a_sel[i]*gAB_a_sel[i]/4
    fAB_a[i] = 1 - fAA_a[i] - fBB_a[i] 
    # MALES
    gAA_a[i] = fAA_a[i]
    gBB_a[i] = fBB_a[i]
    gAB_a[i] = fAB_a[i]
    
    ### X-LINKED LOCUS ###
    # FEMALES
    fAA_x[i] = fAA_x_sel[i]*gA_x_sel[i] + fAB_x_sel[i]*gA_x_sel[i]/2
    fBB_x[i] = fBB_x_sel[i]*gB_x_sel[i] + fAB_x_sel[i]*gB_x_sel[i]/2
    fAB_x[i] = 1 - fAA_x[i] - fBB_x[i]
    # MALES
    gA_x[i] = fAA_x_sel[i] + fAB_x_sel[i]/2
    gB_x[i] = fBB_x_sel[i] + fAB_x_sel[i]/2
  }
  
  for(i in 2:H){
    slope_Auto[i-1] = ((fAA_a[i] + fAB_a[i]/2 + gAA_a[i] +gAB_a[i]/2)/2) - ((fAA_a[i-1] + fAB_a[i-1]/2 + gAA_a[i-1] +gAB_a[i-1]/2)/2)
    slope_X[i-1]= (2*(fAA_x[i] + fAB_x[i]/2)/3 + gA_x[i]/3) - (2*(fAA_x[i-1] + fAB_x[i-1]/2)/3 + gA_x[i-1]/3)
  }
  #### END OF THE LOOP
}

###############################
####### 4. CALCULATIONS #######
###############################

#Frequency of the allele A at the autosomal locus (assumed equal sex ratio)
fA_a = (fAA_a + fAB_a/2 + gAA_a +gAB_a/2) / 2
fA_x = 2*(fAA_x + fAB_x/2)/3 + gA_x/3

####### Reaction diffusion approxiamtion #######
approx_auto = sqrt((sm + sf)*(6*h + 5)/(48*(mf + mm)))
approx_auto
approx.X = sqrt((sf*(6*h + 5) + 8*sm)/(24*(2*mf + mm)))
approx.X
# APPROX X/Auto
approx.X/approx_auto
#######

###### MID-RANGE cline slopes #######
slope_a = fA_a[H/2+1] - fA_a[H/2]
slope_a
slope_x = fA_x[H/2+1] - fA_x[H/2]
slope_x
# MID-RANGE X/Auto
slope_x/slope_a
#######

####### MAXIMUM cline slopes #######
max.slope_Auto = max(slope_Auto)
max.slope_Auto
max.slope_X = max(slope_X)
max.slope_X
# MAXI X/Auto
max.slope_X / max.slope_Auto

##### PLOTS #####
par(mfrow=c(1,1))
# fA AUTOSOMAL (red line)
plot(1:H,(fAA_a+fAB_a/2+gAA_a+gAB_a/2)/2 ,ylim=c(0,1),type="l",col="red")
par(new=TRUE)
#fA X-LINKED (blue line)
plot(1:H, 2*(fAA_x + fAB_x/2)/3 + gA_x/3 ,ylim=c(0,1),type="l",col="blue")

```
