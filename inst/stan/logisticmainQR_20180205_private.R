
// example with standard logistic model
   data{
    int<lower=0> N;
    int<lower=0> p;
    int<lower=0> dx;
    matrix[N, dx] X; // functions like Nx5 matrix (but each column is real[])
    int<lower=0, upper=1> y[N];
   }
   transformed data{
   }
   parameters{
    vector[p] b_;
    real<lower=0> sig_b; 
    real mu_b;
    real b0; // given uniform prior
   }
   transformed parameters{
     vector[p] beta = b_ * sig_b + mu_b;
   }
   model{
   {
     vector[N] mu;
        mu_b ~ normal(0, 5);
        sig_b ~ cauchy(0, 1);
        b_ ~ normal(0, 1);
     
     mu = b0 + X * beta;
     y ~ bernoulli_logit(mu);
   }
   }
   generated quantities{
    real muhat;
    {
    vector[N] r0;
     r0 = inv_logit(b0 + X * beta);
     muhat = mean(r0);
    }
   }


stanmodelQR <- 
"
// example with logistic model with QR decomposition
   data{
    int<lower=0> N;
    int<lower=0> p;
    int<lower=0> dx;
    matrix[N, dx] X; 
    int<lower=0, upper=1> y[N];
   }
   transformed data{
    // Compute, thin, and then scale QR decomposition
    matrix[N, dx] Q = qr_Q(X)[, 1:dx] * N;
    matrix[dx, dx] R = qr_R(X)[1:dx, ] / N;
    matrix[dx, dx] R_inv = inverse(R);
   }
   parameters{
    vector[dx] b_;
    real b0; // given uniform prior
    real sig_b; // given uniform prior
    real mu_b; // given uniform prior
   }
   transformed parameters{
     vector[dx] beta_tilde = b_ * sig_b + mu_b;
     vector[p] beta = R_inv * beta_tilde;
   }
   model{
   {
     vector[N] mu;
     beta ~ normal(0, 100);
     mu = b0 + Q * beta_tilde;
     y ~ bernoulli_logit(mu);
   }
   }
   generated quantities{
    real muhat;
    {
    vector[N] r0;
     r0 = inv_logit(b0 + X * beta);
     muhat = mean(r0);
    }
   }
"

