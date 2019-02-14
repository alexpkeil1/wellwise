
   // example with horseshoe prior
   data{
    int<lower=0> N;
    int<lower=0> p;
    int<lower=0> dx;
    vector[dx] X[N]; // functions like Nx5 matrix (but each column is real[])
    int<lower=0, upper=1> y[N];
   }
   transformed data{
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
    }
   }
   parameters{
    vector[p] beta;
    real<lower=0> sig_b; 
    real mu_b;
    real b0; // given uniform prior
    
   }
   transformed parameters{}
   model{
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
   real rd;
    {
     vector[N] r1;
     vector[N] r0;
     matrix[N,p] D1;
     matrix[N,p] D0;
      D1=D;
      D0=D;
     for(i in 1:N){
        // this is intervention to set exposure 1 to 1.0 versus 0.0 (i.e. both main effect
        // and the self interaction term go to 1.0
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


