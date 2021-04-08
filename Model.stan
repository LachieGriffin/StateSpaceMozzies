data
{
  // Model assumes the number of timesteps (trap deployment) is the same across the regions
  int<lower = 55> t_0; // Time which the modelling should start - has to allow for rainfall
  int<lower = 1> numTimes;  // number of time steps in model
  int<lower = 1> numPast; // number of days in past mosquito's can be born
  
  int<lower = 1> numTraps_Mourilyan;  // number of traps
  int<lower = 1> numReleases_Mourilyan; // release of sterile male mozzies
  int<lower = 1> integralPoints_Mourilyan; //Combination of strips
  int<lower = 1> trapIndex_Mourilyan[integralPoints_Mourilyan]; //Trap index used in dot product for data model
  real<lower = 0.0> integralData_Mourilyan[integralPoints_Mourilyan,(numTimes - t_0 + 1)]; // Calculation for integral strips
  real<lower = 0.0> startTime_Mourilyan[numTraps_Mourilyan, numTimes]; // The proportion of the day transpired when a trap starts trapping
  real<lower = 0.0> endTime_Mourilyan[numTraps_Mourilyan, numTimes]; // The proportion of the day transpired when a trap stops trapping
  real<lower = 0.0> collectionDay_Mourilyan[numTraps_Mourilyan,numTimes]; // indicator matrix for collection date
  int<lower = 0> malesCaught_Mourilyan[integralPoints_Mourilyan]; // trapped males indicator matrix
  int<lower = 0> femalesCaught_Mourilyan[integralPoints_Mourilyan]; // trapped females indicator matrix
  real<lower = 0> rainDays_Mourilyan[numTimes, numPast];
  real<lower = 0> sterileRelease_Mourilyan[numTimes]; // Release numbers of mozzies for dot product

  int<lower = 1> numTraps_GoondiBend;  // number of traps
  int<lower = 1> numReleases_GoondiBend; // release of sterile male mozzies
  int<lower = 1> integralPoints_GoondiBend; //Combination of strips
  int<lower = 1> trapIndex_GoondiBend[integralPoints_GoondiBend]; //Trap index used in dot product for data model
  real<lower = 0.0> integralData_GoondiBend[integralPoints_GoondiBend,(numTimes - t_0 + 1)]; // Calculation for integral strips
  real<lower = 0.0> startTime_GoondiBend[numTraps_GoondiBend, numTimes]; // The proportion of the day transpired when a trap starts trapping
  real<lower = 0.0> endTime_GoondiBend[numTraps_GoondiBend, numTimes]; // The proportion of the day transpired when a trap stops trapping
  real<lower = 0.0> collectionDay_GoondiBend[numTraps_GoondiBend,numTimes]; // indicator matrix for collection date
  int<lower = 0> malesCaught_GoondiBend[integralPoints_GoondiBend]; // trapped males indicator matrix
  int<lower = 0> femalesCaught_GoondiBend[integralPoints_GoondiBend]; // trapped females indicator matrix
  real<lower = 0> rainDays_GoondiBend[numTimes, numPast];
  real<lower = 0> sterileRelease_GoondiBend[numTimes]; // Release numbers of mozzies for dot product

  int<lower = 1> numTraps_SouthJohnstone;  // number of traps
  int<lower = 1> numReleases_SouthJohnstone; // release of sterile male mozzies
  int<lower = 1> integralPoints_SouthJohnstone; //Combination of strips
  int<lower = 1> trapIndex_SouthJohnstone[integralPoints_SouthJohnstone]; //Trap index used in dot product for data model
  real<lower = 0.0> integralData_SouthJohnstone[integralPoints_SouthJohnstone,(numTimes - t_0 + 1)]; // Calculation for integral strips
  real<lower = 0.0> startTime_SouthJohnstone[numTraps_SouthJohnstone, numTimes]; // The proportion of the day transpired when a trap starts trapping
  real<lower = 0.0> endTime_SouthJohnstone[numTraps_SouthJohnstone, numTimes]; // The proportion of the day transpired when a trap stops trapping
  real<lower = 0.0> collectionDay_SouthJohnstone[numTraps_SouthJohnstone,numTimes]; // indicator matrix for collection date
  int<lower = 0> malesCaught_SouthJohnstone[integralPoints_SouthJohnstone]; // trapped males indicator matrix
  int<lower = 0> femalesCaught_SouthJohnstone[integralPoints_SouthJohnstone]; // trapped females indicator matrix
  real<lower = 0> rainDays_SouthJohnstone[numTimes, numPast];
  real<lower = 0> sterileRelease_SouthJohnstone[numTimes]; // Release numbers of mozzies for dot product

  int<lower = 1> numTraps_Wangan;  // number of traps
  int<lower = 1> integralPoints_Wangan; // Total number of strips
  real<lower = 0.0> integralData_Wangan[integralPoints_Wangan,(numTimes - t_0 + 1)]; // Calculation for integral strips
  real<lower = 0.0> startTime_Wangan[numTraps_Wangan, numTimes]; // The proportion of the day transpired when a trap starts trapping
  real<lower = 0.0> endTime_Wangan[numTraps_Wangan, numTimes]; // The proportion of the day transpired when a trap stops trapping
  real<lower = 0.0> collectionDay_Wangan[numTraps_Wangan,numTimes]; // indicator matrix for collection date
  int<lower = 0> malesCaught_Wangan[integralPoints_Wangan]; // trapped males indicator matrix
  int<lower = 0> femalesCaught_Wangan[integralPoints_Wangan]; // trapped females indicator matrix
  real<lower = 0> rainDays_Wangan[numTimes, numPast];
  int<lower = 1> trapIndex_Wangan[integralPoints_Wangan]; //Trap index used in dot product for data model
  
  int<lower = 1> numTraps_Belvedere;  // number of traps
  int<lower = 1> integralPoints_Belvedere; // Total number of strips
  real<lower = 0.0> integralData_Belvedere[integralPoints_Belvedere,(numTimes - t_0 + 1)]; // Calculation for integral strips
  real<lower = 0.0> startTime_Belvedere[numTraps_Belvedere, numTimes]; // The proportion of the day transpired when a trap starts trapping
  real<lower = 0.0> endTime_Belvedere[numTraps_Belvedere, numTimes]; // The proportion of the day transpired when a trap stops trapping
  real<lower = 0.0> collectionDay_Belvedere[numTraps_Belvedere,numTimes]; // indicator matrix for collection date
  int<lower = 0> malesCaught_Belvedere[integralPoints_Belvedere]; // trapped males indicator matrix
  int<lower = 0> femalesCaught_Belvedere[integralPoints_Belvedere]; // trapped females indicator matrix
  real<lower = 0> rainDays_Belvedere[numTimes, numPast];
  int<lower = 1> trapIndex_Belvedere[integralPoints_Belvedere]; //Trap index used in dot product for data model
  
  int<lower = 1> numTraps_SouthInnisfail;  // number of traps
  int<lower = 1> integralPoints_SouthInnisfail; // Total number of strips
  real<lower = 0.0> integralData_SouthInnisfail[integralPoints_SouthInnisfail,(numTimes - t_0 + 1)]; // Calculation for integral strips
  real<lower = 0.0> startTime_SouthInnisfail[numTraps_SouthInnisfail, numTimes]; // The proportion of the day transpired when a trap starts trapping
  real<lower = 0.0> endTime_SouthInnisfail[numTraps_SouthInnisfail, numTimes]; // The proportion of the day transpired when a trap stops trapping
  real<lower = 0.0> collectionDay_SouthInnisfail[numTraps_SouthInnisfail,numTimes]; // indicator matrix for collection date
  int<lower = 0> malesCaught_SouthInnisfail[integralPoints_SouthInnisfail]; // trapped males indicator matrix
  int<lower = 0> femalesCaught_SouthInnisfail[integralPoints_SouthInnisfail]; // trapped females indicator matrix
  real<lower = 0> rainDays_SouthInnisfail[numTimes, numPast];
  int<lower = 1> trapIndex_SouthInnisfail[integralPoints_SouthInnisfail]; //Trap index used in dot product for data model
}

parameters
{
  real<lower = 0.0, upper = 1.0> propSurvived_M;
  real<lower = 0.0, upper = 1.0> propSurvived_F;
  real<lower = 0.0> beta_R;
  real<lower = 0.0> gamma_R;
  real<lower = 0.0> beta_S;

  // Treated region trap params
  real<lower = 0.0, upper = 1.0> trapCatchability_Mourilyan[numTraps_Mourilyan];
  real<lower = 0.0, upper = 1.0> trapSterileDeliveryEfficiency_Mourilyan[numTraps_Mourilyan];
  real<lower = 0.0, upper = 1.0> trapCatchability_GoondiBend[numTraps_GoondiBend];
  real<lower = 0.0, upper = 1.0> trapSterileDeliveryEfficiency_GoondiBend[numTraps_GoondiBend];
  real<lower = 0.0, upper = 1.0> trapCatchability_SouthJohnstone[numTraps_SouthJohnstone];
  real<lower = 0.0, upper = 1.0> trapSterileDeliveryEfficiency_SouthJohnstone[numTraps_SouthJohnstone];
  
  // Untreated region trap params
  real<lower = 0.0, upper = 1.0> trapCatchability_Wangan[numTraps_Wangan];
  real<lower = 0.0, upper = 1.0> trapCatchability_Belvedere[numTraps_Belvedere];
  real<lower = 0.0, upper = 1.0> trapCatchability_SouthInnisfail[numTraps_SouthInnisfail];
  
  // Initial population
  
  real<lower = 0.0, upper = 1.0> mu_catchAbility;
  real<lower = 0.0> sd_catchAbility;
  
  real<lower = 0.0> phi1;
  real<lower = 0.0> phi2;
  real<lower = 0.0> phi3;
  real<lower = 0.0> phi4;
  
  real<lower = 0.0> sigma_innovations;
  
  # set up innovations RV's
  real<lower = 0.0> Mourilyan_innovations[numTimes - t_0 + 1];
  real<lower = 0.0> GoondiBend_innovations[numTimes - t_0 + 1];
  real<lower = 0.0> SouthJohnstone_innovations[numTimes - t_0 + 1];
  real<lower = 0.0> Wangan_innovations[numTimes - t_0 + 1];
  real<lower = 0.0> Belvedere_innovations[numTimes - t_0 + 1];
  real<lower = 0.0> SouthInnisfail_innovations[numTimes - t_0 + 1];
}

model
{
  // Population params 
  real S_Mourilyan[numTimes];
  real M_Mourilyan[numTimes];
  real F_Mourilyan[numTimes];
  
  real S_GoondiBend[numTimes];
  real M_GoondiBend[numTimes];
  real F_GoondiBend[numTimes];
  
  real S_SouthJohnstone[numTimes];
  real M_SouthJohnstone[numTimes];
  real F_SouthJohnstone[numTimes];

  real M_Wangan[numTimes];
  real F_Wangan[numTimes];
  
  real M_Belvedere[numTimes];
  real F_Belvedere[numTimes];
  
  real M_SouthInnisfail[numTimes];
  real F_SouthInnisfail[numTimes];
  
  // Transformed params
  real O_Mourilyan[numTimes];
  real O_GoondiBend[numTimes];
  real O_SouthJohnstone[numTimes];
  
  real p;
  real eta[numTimes];
  real dp;
  int WindowSize;
  real Damp;
  
  real w_R[numPast];
  real w_S[180];
  
  // Priors
  beta_R ~ normal(0.43, 1);
  gamma_R ~ normal(0.0, 100.0);
  beta_S ~ normal(0.0, 100.0);

  propSurvived_M ~ normal(0.8, 0.1);
  propSurvived_F ~ normal(0.8, 0.1);
  
  mu_catchAbility ~ normal(0.0, 0.01);
  sd_catchAbility ~ normal(0.0, 1.0);
  
  sigma_innovations ~ exponential(1);

  // Catchability priors
  trapCatchability_Mourilyan ~ normal(mu_catchAbility, sd_catchAbility);
  trapSterileDeliveryEfficiency_Mourilyan ~ normal(1.0/numTraps_Mourilyan, 0.5);
  
  trapCatchability_GoondiBend ~ normal(mu_catchAbility, sd_catchAbility);
  trapSterileDeliveryEfficiency_GoondiBend ~ normal(1.0/numTraps_GoondiBend, 0.5);
  
  trapCatchability_SouthJohnstone ~ normal(mu_catchAbility, sd_catchAbility);
  trapSterileDeliveryEfficiency_SouthJohnstone ~ normal(1.0/numTraps_SouthJohnstone, 0.5);
  
  trapCatchability_Wangan ~ normal(mu_catchAbility, sd_catchAbility);
  trapCatchability_Belvedere ~ normal(mu_catchAbility, sd_catchAbility);
  trapCatchability_SouthInnisfail ~ normal(mu_catchAbility, sd_catchAbility);
  
  # phi1 male treated, phi2 female treated, phi3 male untreated, phi4 female untreated
  
  phi1 ~ normal(0,10);
  phi2 ~ normal(0,10);
  phi3 ~ normal(0,10);
  phi4 ~ normal(0,10);
  
  // Initialise population vectors
  for(t in 1:(numTimes - t_0 + 1))
  {
    Mourilyan_innovations[t] ~ lognormal(log(1/(sqrt(sigma_innovations^2 + 1))),sqrt(log(1 + sigma_innovations^2)));
    GoondiBend_innovations[t] ~ lognormal(log(1/(sqrt(sigma_innovations^2 + 1))),sqrt(log(1 + sigma_innovations^2)));
    SouthJohnstone_innovations[t] ~ lognormal(log(1/(sqrt(sigma_innovations^2 + 1))),sqrt(log(1 + sigma_innovations^2)));
    Wangan_innovations[t] ~ lognormal(log(1/(sqrt(sigma_innovations^2 + 1))),sqrt(log(1 + sigma_innovations^2)));
    Belvedere_innovations[t] ~ lognormal(log(1/(sqrt(sigma_innovations^2 + 1))),sqrt(log(1 + sigma_innovations^2)));
    SouthInnisfail_innovations[t] ~ lognormal(log(1/(sqrt(sigma_innovations^2 + 1))),sqrt(log(1 + sigma_innovations^2)));
  }
  
  // Set up inital population as normal
  for(t in 1:numTimes)
  {
    S_Mourilyan[t] = 0.0;
    M_Mourilyan[t] = 0.0;
    F_Mourilyan[t] = 0.0;
    S_GoondiBend[t] = 0.0;
    M_GoondiBend[t] = 0.0;
    F_GoondiBend[t] = 0.0;
    S_SouthJohnstone[t] = 0.0;
    M_SouthJohnstone[t] = 0.0;
    F_SouthJohnstone[t] = 0.0;
    M_Wangan[t] = 0.0;
    F_Wangan[t] = 0.0;
    M_Belvedere[t] = 0.0;
    F_Belvedere[t] = 0.0;
    M_SouthInnisfail[t] = 0.0;
    F_SouthInnisfail[t] = 0.0;
    
    O_Mourilyan[t] = 0.0;
    O_GoondiBend[t] = 0.0;
    O_SouthJohnstone[t] = 0.0;
  }
  
  ///////////////////////////////////////
  ///////// SHARED PARAMETERS ///////////
  ///////////////////////////////////////
  
  for(i in 1:42)
  {
  // Rainfall gamma weight
    w_R[i] = gamma_R*(beta_R)*exp(-beta_R*i);
  }
  for(i in 1:180)
  {
  // sterile population weight
   w_S[i] = (beta_S)*exp(-beta_S*i);
  }
  for(i in 1:180)
  {
    w_S[i] = w_S[i]/sum(w_S);
  }
  
  /////////////////////////////////////////
  ///////// Mourilyan [TREATED] ///////////
  /////////////////////////////////////////
  
  // Calculate process model
  for(t in t_0:numTimes)
  {
    // set size for window
    WindowSize = min(189,t);
    
    dp = dot_product(rainDays_Mourilyan[t - 9, 1:numPast], w_R[1:numPast]);
    S_Mourilyan[t] = S_Mourilyan[t-1]*propSurvived_M + sterileRelease_Mourilyan[t];
    
    if(t > t_0)
    {
      O_Mourilyan[t-1] = (S_Mourilyan[t-1])/(M_Mourilyan[t-1]);
    }
    
    Damp = dot_product(w_S[1:(WindowSize-9)], O_Mourilyan[(t - WindowSize + 1):(t-9)]);

    p = 1.0/(1.0 + Damp);
    
    eta[t] = dp*p;
    
    if(t == t_0)
    {
      M_Mourilyan[t] = (eta[t]/(1-propSurvived_M))*Mourilyan_innovations[t-t_0+1];
      F_Mourilyan[t] = (eta[t]/(1-propSurvived_F))*Mourilyan_innovations[t-t_0+1];
    }
    else
    {
      M_Mourilyan[t] = (M_Mourilyan[t-1]*propSurvived_M + eta[t])*Mourilyan_innovations[t-t_0+1];
      F_Mourilyan[t] = (F_Mourilyan[t-1]*propSurvived_F + eta[t])*Mourilyan_innovations[t-t_0+1];
    }
  }
  
  // Data model
  for(j in 1:integralPoints_Mourilyan)
  {
      malesCaught_Mourilyan[j] ~ neg_binomial_2(trapCatchability_Mourilyan[trapIndex_Mourilyan[j]]*(dot_product(M_Mourilyan[t_0:numTimes],integralData_Mourilyan[j,1:(numTimes - t_0 + 1)]) + dot_product(S_Mourilyan[t_0:numTimes],integralData_Mourilyan[j,1:(numTimes - t_0 + 1)])*trapSterileDeliveryEfficiency_Mourilyan[trapIndex_Mourilyan[j]]),phi1);
      femalesCaught_Mourilyan[j] ~ neg_binomial_2(trapCatchability_Mourilyan[trapIndex_Mourilyan[j]]*dot_product(F_Mourilyan[t_0:numTimes], integralData_Mourilyan[j,1:(numTimes - t_0 + 1)]),phi2);
  }
  
  ///////////////////////////////////////////
  ///////// Goondi Bend [TREATED] ///////////
  ///////////////////////////////////////////
  
  // Calculate process model
  for(t in t_0:numTimes)
  {
    // set size for window
    WindowSize = min(189,t);
    
    dp = dot_product(rainDays_GoondiBend[t - 9, 1:numPast], w_R[1:numPast]);
    S_GoondiBend[t] = S_GoondiBend[t-1]*propSurvived_M + sterileRelease_GoondiBend[t];
    
    if(t > t_0)
    {
      O_GoondiBend[t-1] = (S_GoondiBend[t-1])/(M_GoondiBend[t-1]);
    }
    
    Damp = dot_product(w_S[1:(WindowSize-9)], O_GoondiBend[(t - WindowSize + 1):(t-9)]);

    p = 1.0/(1.0 + Damp);
    
    eta[t] = dp*p;
    
    if(t == t_0)
    {
      M_GoondiBend[t] = (eta[t]/(1-propSurvived_M))*GoondiBend_innovations[t-t_0+1];
      F_GoondiBend[t] = (eta[t]/(1-propSurvived_F))*GoondiBend_innovations[t-t_0+1];
    }
    else
    {
      M_GoondiBend[t] = (M_GoondiBend[t-1]*propSurvived_M + eta[t])*GoondiBend_innovations[t-t_0+1];
      F_GoondiBend[t] = (F_GoondiBend[t-1]*propSurvived_F + eta[t])*GoondiBend_innovations[t-t_0+1];
    }
  }
  
  // Data model
  for(j in 1:integralPoints_GoondiBend)
  {
      malesCaught_GoondiBend[j] ~ neg_binomial_2(trapCatchability_GoondiBend[trapIndex_GoondiBend[j]]*(dot_product(M_GoondiBend[t_0:numTimes],integralData_GoondiBend[j,1:(numTimes - t_0 + 1)]) + dot_product(S_GoondiBend[t_0:numTimes],integralData_GoondiBend[j,1:(numTimes - t_0 + 1)])*trapSterileDeliveryEfficiency_GoondiBend[trapIndex_GoondiBend[j]]),phi1);
      femalesCaught_GoondiBend[j] ~ neg_binomial_2(trapCatchability_GoondiBend[trapIndex_GoondiBend[j]]*dot_product(F_GoondiBend[t_0:numTimes], integralData_GoondiBend[j,1:(numTimes - t_0 + 1)]),phi2);
  }
  
  ///////////////////////////////////////////////
  ///////// South Johnstone [TREATED] ///////////
  ///////////////////////////////////////////////
  
  // Calculate process model
  for(t in t_0:numTimes)
  {
    // set size for window
    WindowSize = min(189,t);
    
    dp = dot_product(rainDays_SouthJohnstone[t - 9, 1:numPast], w_R[1:numPast]);
    S_SouthJohnstone[t] = S_SouthJohnstone[t-1]*propSurvived_M + sterileRelease_SouthJohnstone[t];
    
    if(t > t_0)
    {
      O_SouthJohnstone[t-1] = (S_SouthJohnstone[t-1])/(M_SouthJohnstone[t-1]);
    }
    
    Damp = dot_product(w_S[1:(WindowSize-9)], O_SouthJohnstone[(t - WindowSize + 1):(t-9)]);
    p = 1.0/(1.0 + Damp);
    
    eta[t] = dp*p;
    
    if(t == t_0)
    {
      M_SouthJohnstone[t] = (eta[t]/(1-propSurvived_M))*SouthJohnstone_innovations[t-t_0+1];
      F_SouthJohnstone[t] = (eta[t]/(1-propSurvived_F))*SouthJohnstone_innovations[t-t_0+1];
    }
    else
    {
      M_SouthJohnstone[t] = (M_SouthJohnstone[t-1]*propSurvived_M + eta[t])*SouthJohnstone_innovations[t-t_0+1];
      F_SouthJohnstone[t] = (F_SouthJohnstone[t-1]*propSurvived_F + eta[t])*SouthJohnstone_innovations[t-t_0+1];
    }
  }
  
  // Data model
  for(j in 1:integralPoints_SouthJohnstone)
  {
      malesCaught_SouthJohnstone[j] ~ neg_binomial_2(trapCatchability_SouthJohnstone[trapIndex_SouthJohnstone[j]]*(dot_product(M_SouthJohnstone[t_0:numTimes],integralData_SouthJohnstone[j,1:(numTimes - t_0 + 1)]) + dot_product(S_SouthJohnstone[t_0:numTimes],integralData_SouthJohnstone[j,1:(numTimes - t_0 + 1)])*trapSterileDeliveryEfficiency_SouthJohnstone[trapIndex_SouthJohnstone[j]]),phi1);
      femalesCaught_SouthJohnstone[j] ~ neg_binomial_2(trapCatchability_SouthJohnstone[trapIndex_SouthJohnstone[j]]*dot_product(F_SouthJohnstone[t_0:numTimes], integralData_SouthJohnstone[j,1:(numTimes - t_0 + 1)]),phi2);
  }
  
  ///////////////////////////////////////////
  /////////// Wangan [UNTREATED] ////////////
  ///////////////////////////////////////////
  
  // Calculate process model
  for(t in t_0:numTimes)
  {
    dp = dot_product(rainDays_Wangan[t-9, 1:numPast], w_R[1:numPast]);
    eta[t] = dp;
    
    #Wangan_innovations[t-t_0+1] ~ normal(0, sigma_innovations);
    
    if(t == t_0)
    {
      M_Wangan[t] = (eta[t]/(1-propSurvived_M))*Wangan_innovations[t-t_0+1];
      F_Wangan[t] = (eta[t]/(1-propSurvived_F))*Wangan_innovations[t-t_0+1];
    }
    else
    {
      M_Wangan[t] = (M_Wangan[t-1]*propSurvived_M + eta[t])*Wangan_innovations[t-t_0+1];
      F_Wangan[t] = (F_Wangan[t-1]*propSurvived_F + eta[t])*Wangan_innovations[t-t_0+1];
    }
  }
  // Data model
  for(j in 1:integralPoints_Wangan)
  {
      malesCaught_Wangan[j] ~ neg_binomial_2(trapCatchability_Wangan[trapIndex_Wangan[j]]*dot_product(M_Wangan[t_0:numTimes],integralData_Wangan[j,1:(numTimes - t_0 + 1)]),phi3);
      femalesCaught_Wangan[j] ~ neg_binomial_2(trapCatchability_Wangan[trapIndex_Wangan[j]]*dot_product(F_Wangan[t_0:numTimes],integralData_Wangan[j,1:(numTimes - t_0 + 1)]),phi4);
  }
  
  //////////////////////////////////////////////
  /////////// Belvedere [UNTREATED] ////////////
  //////////////////////////////////////////////
  
  // Calculate process model
  for(t in t_0:numTimes)
  {
    dp = dot_product(rainDays_Belvedere[t-9, 1:numPast], w_R[1:numPast]);
    eta[t] = dp;
    
    if(t == t_0)
    {
      M_Belvedere[t] = (eta[t]/(1-propSurvived_M))*Belvedere_innovations[t-t_0+1];
      F_Belvedere[t] = (eta[t]/(1-propSurvived_F))*Belvedere_innovations[t-t_0+1];
    }
    else
    {
      M_Belvedere[t] = (M_Belvedere[t-1]*propSurvived_M + eta[t])*Belvedere_innovations[t-t_0+1];
      F_Belvedere[t] = (F_Belvedere[t-1]*propSurvived_F + eta[t])*Belvedere_innovations[t-t_0+1];
    }
  }
  // Data model
  for(j in 1:integralPoints_Belvedere)
  {
      malesCaught_Belvedere[j] ~ neg_binomial_2(trapCatchability_Belvedere[trapIndex_Belvedere[j]]*dot_product(M_Belvedere[t_0:numTimes],integralData_Belvedere[j,1:(numTimes - t_0 + 1)]),phi3);
      femalesCaught_Belvedere[j] ~ neg_binomial_2(trapCatchability_Belvedere[trapIndex_Belvedere[j]]*dot_product(F_Belvedere[t_0:numTimes],integralData_Belvedere[j,1:(numTimes - t_0 + 1)]),phi4);
  }
  
  ///////////////////////////////////////////////////
  /////////// SouthInnisfail [UNTREATED] ////////////
  ///////////////////////////////////////////////////
  
  // Calculate process model
  for(t in t_0:numTimes)
  {
    dp = dot_product(rainDays_SouthInnisfail[t-9, 1:numPast], w_R[1:numPast]);
    eta[t] = dp;
    
    if(t == t_0)
    {
      M_SouthInnisfail[t] = (eta[t]/(1-propSurvived_M))*SouthInnisfail_innovations[t-t_0+1];
      F_SouthInnisfail[t] = (eta[t]/(1-propSurvived_F))*SouthInnisfail_innovations[t-t_0+1];
    }
    else
    {
      M_SouthInnisfail[t] = (M_SouthInnisfail[t-1]*propSurvived_M + eta[t])*SouthInnisfail_innovations[t-t_0+1];
      F_SouthInnisfail[t] = (F_SouthInnisfail[t-1]*propSurvived_F + eta[t])*SouthInnisfail_innovations[t-t_0+1];
    }
  }
  // Data model
  for(j in 1:integralPoints_SouthInnisfail)
  {
      malesCaught_SouthInnisfail[j] ~ neg_binomial_2(trapCatchability_SouthInnisfail[trapIndex_SouthInnisfail[j]]*dot_product(M_SouthInnisfail[t_0:numTimes],integralData_SouthInnisfail[j,1:(numTimes - t_0 + 1)]),phi3);
      femalesCaught_SouthInnisfail[j] ~ neg_binomial_2(trapCatchability_SouthInnisfail[trapIndex_SouthInnisfail[j]]*dot_product(F_SouthInnisfail[t_0:numTimes],integralData_SouthInnisfail[j,1:(numTimes - t_0 + 1)]),phi4);
  }
}
generated quantities
{
  // Population params 
  real S_Mourilyan[numTimes];
  real M_Mourilyan[numTimes];
  real F_Mourilyan[numTimes];
  
  real S_GoondiBend[numTimes];
  real M_GoondiBend[numTimes];
  real F_GoondiBend[numTimes];
  
  real S_SouthJohnstone[numTimes];
  real M_SouthJohnstone[numTimes];
  real F_SouthJohnstone[numTimes];

  real M_Wangan[numTimes];
  real F_Wangan[numTimes];
  
  real M_Belvedere[numTimes];
  real F_Belvedere[numTimes];
  
  real M_SouthInnisfail[numTimes];
  real F_SouthInnisfail[numTimes];
  
  // Transformed params
  real O_Mourilyan[numTimes];
  real O_GoondiBend[numTimes];
  real O_SouthJohnstone[numTimes];
  
  real p;
  real eta[numTimes];
  real dp;
  int WindowSize;
  real Damp;
  
  real w_R[numPast];
  real w_S[180];
  
  // Set up inital population as normal
  for(t in 1:numTimes)
  {
    S_Mourilyan[t] = 0.0;
    M_Mourilyan[t] = 0.0;
    F_Mourilyan[t] = 0.0;
    S_GoondiBend[t] = 0.0;
    M_GoondiBend[t] = 0.0;
    F_GoondiBend[t] = 0.0;
    S_SouthJohnstone[t] = 0.0;
    M_SouthJohnstone[t] = 0.0;
    F_SouthJohnstone[t] = 0.0;
    M_Wangan[t] = 0.0;
    F_Wangan[t] = 0.0;
    M_Belvedere[t] = 0.0;
    F_Belvedere[t] = 0.0;
    M_SouthInnisfail[t] = 0.0;
    F_SouthInnisfail[t] = 0.0;
    
    O_Mourilyan[t] = 0.0;
    O_GoondiBend[t] = 0.0;
    O_SouthJohnstone[t] = 0.0;
  }
  
  ///////////////////////////////////////
  ///////// SHARED PARAMETERS ///////////
  ///////////////////////////////////////
  
  for(i in 1:42)
  {
  // Rainfall gamma weight
    w_R[i] = gamma_R*(beta_R)*exp(-beta_R*i);
  }
  for(i in 1:180)
  {
  // sterile population weight
   w_S[i] = (beta_S)*exp(-beta_S*i);
  }
  for(i in 1:180)
  {
    w_S[i] = w_S[i]/sum(w_S);
  }
  
  /////////////////////////////////////////
  ///////// Mourilyan [TREATED] ///////////
  /////////////////////////////////////////
  
  // Calculate process model
  for(t in t_0:numTimes)
  {
    // set size for window
    WindowSize = min(189,t);
    
    dp = dot_product(rainDays_Mourilyan[t - 9, 1:numPast], w_R[1:numPast]);
    S_Mourilyan[t] = S_Mourilyan[t-1]*propSurvived_M + sterileRelease_Mourilyan[t];
    
    if(t > t_0)
    {
      O_Mourilyan[t-1] = (S_Mourilyan[t-1])/(M_Mourilyan[t-1]);
    }
    
    Damp = dot_product(w_S[1:(WindowSize-9)], O_Mourilyan[(t - WindowSize + 1):(t-9)]);

    p = 1.0/(1.0 + Damp);
    
    eta[t] = dp*p;
    
    if(t == t_0)
    {
      M_Mourilyan[t] = (eta[t]/(1-propSurvived_M))*Mourilyan_innovations[t-t_0+1];
      F_Mourilyan[t] = (eta[t]/(1-propSurvived_F))*Mourilyan_innovations[t-t_0+1];
    }
    else
    {
      M_Mourilyan[t] = (M_Mourilyan[t-1]*propSurvived_M + eta[t])*Mourilyan_innovations[t-t_0+1];
      F_Mourilyan[t] = (F_Mourilyan[t-1]*propSurvived_F + eta[t])*Mourilyan_innovations[t-t_0+1];
    }
  }
  
  ///////////////////////////////////////////
  ///////// Goondi Bend [TREATED] ///////////
  ///////////////////////////////////////////
  
  // Calculate process model
  for(t in t_0:numTimes)
  {
    // set size for window
    WindowSize = min(189,t);
    
    dp = dot_product(rainDays_GoondiBend[t - 9, 1:numPast], w_R[1:numPast]);
    S_GoondiBend[t] = S_GoondiBend[t-1]*propSurvived_M + sterileRelease_GoondiBend[t];
    
    if(t > t_0)
    {
      O_GoondiBend[t-1] = (S_GoondiBend[t-1])/(M_GoondiBend[t-1]);
    }
    
    Damp = dot_product(w_S[1:(WindowSize-9)], O_GoondiBend[(t - WindowSize + 1):(t-9)]);

    p = 1.0/(1.0 + Damp);
    
    eta[t] = dp*p;
    
    if(t == t_0)
    {
      M_GoondiBend[t] = (eta[t]/(1-propSurvived_M))*GoondiBend_innovations[t-t_0+1];
      F_GoondiBend[t] = (eta[t]/(1-propSurvived_F))*GoondiBend_innovations[t-t_0+1];
    }
    else
    {
      M_GoondiBend[t] = (M_GoondiBend[t-1]*propSurvived_M + eta[t])*GoondiBend_innovations[t-t_0+1];
      F_GoondiBend[t] = (F_GoondiBend[t-1]*propSurvived_F + eta[t])*GoondiBend_innovations[t-t_0+1];
    }
  }
  
  ///////////////////////////////////////////////
  ///////// South Johnstone [TREATED] ///////////
  ///////////////////////////////////////////////
  
  // Calculate process model
  for(t in t_0:numTimes)
  {
    // set size for window
    WindowSize = min(189,t);
    
    dp = dot_product(rainDays_SouthJohnstone[t - 9, 1:numPast], w_R[1:numPast]);
    S_SouthJohnstone[t] = S_SouthJohnstone[t-1]*propSurvived_M + sterileRelease_SouthJohnstone[t];
    
    if(t > t_0)
    {
      O_SouthJohnstone[t-1] = (S_SouthJohnstone[t-1])/(M_SouthJohnstone[t-1]);
    }
    
    Damp = dot_product(w_S[1:(WindowSize-9)], O_SouthJohnstone[(t - WindowSize + 1):(t-9)]);
    p = 1.0/(1.0 + Damp);
    
    eta[t] = dp*p;
    
    if(t == t_0)
    {
      M_SouthJohnstone[t] = (eta[t]/(1-propSurvived_M))*SouthJohnstone_innovations[t-t_0+1];
      F_SouthJohnstone[t] = (eta[t]/(1-propSurvived_F))*SouthJohnstone_innovations[t-t_0+1];
    }
    else
    {
      M_SouthJohnstone[t] = (M_SouthJohnstone[t-1]*propSurvived_M + eta[t])*SouthJohnstone_innovations[t-t_0+1];
      F_SouthJohnstone[t] = (F_SouthJohnstone[t-1]*propSurvived_F + eta[t])*SouthJohnstone_innovations[t-t_0+1];
    }
  }
  
  ///////////////////////////////////////////
  /////////// Wangan [UNTREATED] ////////////
  ///////////////////////////////////////////
  
  // Calculate process model
  for(t in t_0:numTimes)
  {
    dp = dot_product(rainDays_Wangan[t-9, 1:numPast], w_R[1:numPast]);
    eta[t] = dp;
    
    #Wangan_innovations[t-t_0+1] ~ normal(0, sigma_innovations);
    
    if(t == t_0)
    {
      M_Wangan[t] = (eta[t]/(1-propSurvived_M))*Wangan_innovations[t-t_0+1];
      F_Wangan[t] = (eta[t]/(1-propSurvived_F))*Wangan_innovations[t-t_0+1];
    }
    else
    {
      M_Wangan[t] = (M_Wangan[t-1]*propSurvived_M + eta[t])*Wangan_innovations[t-t_0+1];
      F_Wangan[t] = (F_Wangan[t-1]*propSurvived_F + eta[t])*Wangan_innovations[t-t_0+1];
    }
  }
  
  //////////////////////////////////////////////
  /////////// Belvedere [UNTREATED] ////////////
  //////////////////////////////////////////////
  
  // Calculate process model
  for(t in t_0:numTimes)
  {
    dp = dot_product(rainDays_Belvedere[t-9, 1:numPast], w_R[1:numPast]);
    eta[t] = dp;
    
    if(t == t_0)
    {
      M_Belvedere[t] = (eta[t]/(1-propSurvived_M))*Belvedere_innovations[t-t_0+1];
      F_Belvedere[t] = (eta[t]/(1-propSurvived_F))*Belvedere_innovations[t-t_0+1];
    }
    else
    {
      M_Belvedere[t] = (M_Belvedere[t-1]*propSurvived_M + eta[t])*Belvedere_innovations[t-t_0+1];
      F_Belvedere[t] = (F_Belvedere[t-1]*propSurvived_F + eta[t])*Belvedere_innovations[t-t_0+1];
    }
  }
  
  ///////////////////////////////////////////////////
  /////////// SouthInnisfail [UNTREATED] ////////////
  ///////////////////////////////////////////////////
  
  // Calculate process model
  for(t in t_0:numTimes)
  {
    dp = dot_product(rainDays_SouthInnisfail[t-9, 1:numPast], w_R[1:numPast]);
    eta[t] = dp;
    
    if(t == t_0)
    {
      M_SouthInnisfail[t] = (eta[t]/(1-propSurvived_M))*SouthInnisfail_innovations[t-t_0+1];
      F_SouthInnisfail[t] = (eta[t]/(1-propSurvived_F))*SouthInnisfail_innovations[t-t_0+1];
    }
    else
    {
      M_SouthInnisfail[t] = (M_SouthInnisfail[t-1]*propSurvived_M + eta[t])*SouthInnisfail_innovations[t-t_0+1];
      F_SouthInnisfail[t] = (F_SouthInnisfail[t-1]*propSurvived_F + eta[t])*SouthInnisfail_innovations[t-t_0+1];
    }
  }
}
