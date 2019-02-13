
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
    vector[dx] X[N]; // functions like Nx5 matrix (but each column is real[])
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
    //real meany;
    //vector[N] ycen;
    matrix[N,p] D;

    //meany = mean(y);
    //ycen = y-meany;
   // todo: transform y to center!
    for(c in 1:dx){
     D[,c] = to_vector(X[,c]);
     //D[,c] = X[,c];
    }
    for(c in (dx+1):p){
     D[,c] = to_vector(X[,c-5]) .* to_vector(X[,c-5]);
     //D[,c] = X[,c-5] .* X[,c-5];
    }
   }
   parameters{
  // This programming block will be used for:
  // defining the parameters in the statstical model, including hyperparameters
    vector<lower=0>[p] lambda; // local shrinkage
    vector[p] beta;
    real<lower=0> sig_b; 
    real mu_b;
    real b0; // given uniform prior
    
   }
   transformed parameters{
  // This programming block will be used for:
  //   using variable transformations to allow for more efficient MCMC sampling
  //  e.g. see the "matt-trick" in stan user group forums
  }
   model{
  // This programming block will be used for:
  //  the actual statistical model which estimates the risk of birth outcomes
  //  conditional on: exposures, confounders, hyperparameters
    {
     vector[N] mu;
     b0 ~ normal(0, 5);
     sig_b ~ cauchy(0, 1);
     beta ~ normal(mu_b, sig_b);
     mu = b0 + D * beta;
     y ~ bernoulli_logit(mu);
    }
   }
   generated quantities{
  // This programming block will be used for:
  //  risk estimation (e.g. getting population average risk under a specific intervention)
  //  effect estimation (e.g. getting risk difference contrasting two interventions)
  //  intervention modifiers: e.g. water chemistry parameters assumed to be known
  //  cost-effectiveness g-formula estimation (e.g. add in cost per treated individual)
   real rd;
    {
     vector[N] r1;
     vector[N] r0;
     matrix[N,p] D1;
     matrix[N,p] D0;
      D1=D;
      D0=D;
     for(i in 1:N){
        // This is intervention to set exposure 1 to 1.0 versus 0.0 (i.e. both main effect
        // and the self interaction term go to 1.0.
        // Will need to generalize this block to allow for more general statistical
        //    models.
        D1[i,1] = 1.0;
        D1[i,6] = 1.0;
        D0[i,1] = 0.0;
        D0[i,6] = 0.0;
      }
    
      r1 =  inv_logit(b0 + D1 * beta);
      r0 =  inv_logit(b0 + D0 * beta);
      rd = mean(r1)-mean(r0);
    }
   }


