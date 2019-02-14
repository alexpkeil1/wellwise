
   // example with horseshoe prior
   data{
   // This programming block will be used for:
   //  reading in data in specific formats/types
   //   e.g. a 1/0 variable would be most efficiently coded if we define it as
   //   an *int* variable with lower and upper bounds. A continuous
   //   variable could be coded as a *real* or as a *vector*, which will 
   //   determine the efficiency of the coding (i.e. we could probably use
   //   these interchangeably for the same statistical model, but the 
   //   model code will vary (we may be able to vectorize some operations)
   //   and one implementation may be more efficient (e.g. vectorized code
   //   will often be slightly faster than looping)
    int<lower=0> N;
    int<lower=0> p;
    int<lower=0> dx;
    matrix[N,dx] X; // functions like Nx5 matrix (but each column is real[])
    int<lower=0, upper=1> y[N];
   }
   transformed data{
    // This programming block will be used for:
    // any necessary transformations to help make code/sampling more efficient
    // e.g. continuous outcomes and exposures will usually benefit from being
    //  centered and divided by their standard deviation, which helps substantially
    //  with the efficiency of MCMC. Note that such transformations will require
    //  back-transformations in the *generated quantitites* programming block, which
    //  may or may not be tricky.

    // here we use it for convenience when reading in data
    vector[N] Arsenic = col(X,1);
    vector[N] Manganese = col(X,2);
    vector[N] Lead = col(X,3);
    vector[N] Cadmium = col(X,4);
    vector[N] Copper = col(X,5);
   }
   parameters{
    // This programming block will be used for:
    // defining the parameters in the statstical model, including hyperparameters
    vector<lower=0>[p] lambda; // local shrinkage
    real<lower=0> tau; // global shrinkage
    vector[p] b_;
    real b0; // given uniform prior
    
   }
   transformed parameters{
  // This programming block will be used for:
  //   using variable transformations to allow for more efficient MCMC sampling
  //  e.g. see the "matt-trick" in stan user group forums
    // one example is using this block to simplify sampling under the horseshoe prior
    // e.g. b_ is easy to sample, N(0,1)
    vector[p] beta = tau * lambda .* b_;
  }
   model{
  // This programming block will be used for:
  //  the actual statistical model which estimates the risk of birth outcomes
  //  conditional on: exposures, confounders, hyperparameters
    {
     vector[N] mu;
     b0 ~ normal(0, 5);
     // horseshoe priors
     tau ~ cauchy(0, 1);
     lambda ~ cauchy(0, 1);
     b_ ~ normal(0, 1);
     mu = b0 + beta[1]*Arsenic + beta[2]*Manganese + beta[3]*Lead + beta[4]*Cadmium + beta[5]*Copper +
          beta[6]*Arsenic .* Arsenic + beta[7]*Arsenic .* Manganese + beta[8]*Arsenic .* Lead + beta[9]*Arsenic .* Cadmium +
          beta[10]*Arsenic .* Copper + beta[11]*Manganese .* Manganese + beta[12]*Lead .* Manganese + beta[13]*Cadmium .* Manganese +
          beta[14]*Copper .* Manganese + beta[15]*Lead .* Lead + beta[16]*Cadmium .* Lead +
          beta[17]*Copper .* Lead + beta[18]*Cadmium .* Cadmium + beta[19]*Cadmium .* Copper + beta[20]*Copper .* Copper;
;
     y ~ bernoulli_logit(mu);
    }
   }
   generated quantities{
  // This programming block will be used for:
  //  risk estimation (e.g. getting population average risk under a specific intervention)
  //  effect estimation (e.g. getting risk difference contrasting two interventions)
  //  intervention modifiers: e.g. water chemistry parameters assumed to be known
  //  cost-effectiveness g-formula estimation (e.g. add in cost per treated individual)
   real rd[5];
   {
     vector[N] r1;
     vector[N] r0;
     
     matrix[5,4] intprop = [[.98, .98, .98, .98], // filter 1 reduces metal 1 (Arsenic) by 98%, metal 2 (Cadmium) by 98%, metal 3 (Mn) by 98%, metal 4 (Pb) by 98% 
                            [0, .985, 0, .965],
                            [0, .918, 0, 0],
                            [.994, .98, .98, .987],
                            [0, 0, 0, .99]]; // ROWS = 5 filters, cols = 4 metals
     
     r0 = inv_logit(b0 + beta[1]*Arsenic + beta[2]*Manganese + beta[3]*Lead + beta[4]*Cadmium + beta[5]*Copper +
     beta[6]*Arsenic .* Arsenic + beta[7]*Arsenic .* Manganese + beta[8]*Arsenic .* Lead + beta[9]*Arsenic .* Cadmium +
     beta[10]*Arsenic .* Copper + beta[11]*Manganese .* Manganese + beta[12]*Lead .* Manganese + beta[13]*Cadmium .* Manganese +
     beta[14]*Copper .* Manganese + beta[15]*Lead .* Lead + beta[16]*Cadmium .* Lead +
     beta[17]*Copper .* Lead + beta[18]*Cadmium .* Cadmium + beta[19]*Cadmium .* Copper + beta[20]*Copper .* Copper);
     
     for(j in 1:5) { //looping over interventions
        // This is intervention to set some expsoures to the level expected under filtration
       r1 = inv_logit(b0 + beta[1]*((1-intprop[j,1])*Arsenic) + beta[2]*((1-intprop[j,3])*Manganese) + beta[3]*((1-intprop[j,4])*Lead) + 
           beta[4]*((1-intprop[j,2])*Cadmium) + beta[5]*Copper + beta[6]*((1-intprop[j,1])*Arsenic) .* ((1-intprop[j,1])*Arsenic) + 
           beta[7]*((1-intprop[j,1])*Arsenic) .* ((1-intprop[j,3])*Manganese) + beta[8]*((1-intprop[j,1])*Arsenic) .* ((1-intprop[j,4])*Lead) +
           beta[9]*((1-intprop[j,1])*Arsenic) .* ((1-intprop[j,2])*Cadmium) +  beta[10]*((1-intprop[j,1])*Arsenic) .* Copper + 
           beta[11]*((1-intprop[j,3])*Manganese) .* ((1-intprop[j,3])*Manganese) + beta[12]*((1-intprop[j,4])*Lead) .* ((1-intprop[j,3])*Manganese) + 
           beta[13]*((1-intprop[j,2])*Cadmium) .* ((1-intprop[j,3])*Manganese) + beta[14]*Copper .* ((1-intprop[j,3])*Manganese) + 
           beta[15]*((1-intprop[j,4])*Lead) .* ((1-intprop[j,4])*Lead) + beta[16]*((1-intprop[j,2])*Cadmium) .* ((1-intprop[j,4])*Lead) +  
           beta[17]*Copper .* ((1-intprop[j,4])*Lead) + beta[18]*((1-intprop[j,2])*Cadmium) .* ((1-intprop[j,2])*Cadmium) + 
           beta[19]*((1-intprop[j,2])*Cadmium) .* Copper + beta[20]*Copper .* Copper);
      rd[j] = mean(r1)-mean(r0);
     } 
    
   }
}


