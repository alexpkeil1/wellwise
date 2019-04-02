
// example with standard logistic model
   data{
    int<lower=0> N;
    int<lower=0> p;
    int<lower=0> dx;
    matrix[N, dx] X; // functions like Nx5 matrix (but each column is real[])
    int<lower=0, upper=1> y[N];
   }
   transformed data{
    vector[N] Arsenic = col(X,1);
    vector[N] Manganese = col(X,2);
    vector[N] Lead = col(X,3);
    vector[N] Cadmium = col(X,4);
    vector[N] Copper = col(X,5);
   }
   parameters{
    vector[p] b_;
    //vector[p] beta;
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
        //sig_b ~ cauchy(0, 0.9);// divergent with half cauchy(0,1), cauchy(0,0.9), ok with 0.1
        sig_b ~ student_t(1, 0, 1.0);// divergent at df < 25
        //sig_b ~ inv_gamma(10, 1);// non-divergent
        b_ ~ normal(0, 1);
        //beta ~ normal(0, 1); //works
        //beta ~ normal(0, sig_b);// divergent with half cauchy(0,1)
        //beta ~ normal(mu_b, 1);// works
     
     mu = b0 + beta[1]*Arsenic + beta[2]*Manganese + beta[3]*Lead + beta[4]*Cadmium + beta[5]*Copper;
     y ~ bernoulli_logit(mu);
   }
   }
   generated quantities{
   real rd[5];
   real r0m;
   real r1m[5];
   {
     vector[N] r1;
     vector[N] r0;
     
     matrix[5,4] intprop = [[.98, .98, .98, .98], // filter 1 reduces metal 1 (Arsenic) by 98%, metal 2 (Cadmium) by 98%, metal 3 (Mn) by 98%, metal 4 (Pb) by 98% 
                            [0, .985, 0, .965],
                            [0, .918, 0, 0],
                            [.994, .98, .98, .987],
                            [0, 0, 0, .99]]; // ROWS = 5 filters, cols = 4 metals
     
     r0 = inv_logit(b0 + beta[1]*Arsenic + beta[2]*Manganese + beta[3]*Lead + beta[4]*Cadmium + beta[5]*Copper);
     r0m = mean(r0);
     for(j in 1:5) { //looping over interventions
       r1 = inv_logit(b0 + beta[1]*((1-intprop[j,1])*Arsenic) + beta[2]*((1-intprop[j,3])*Manganese) + beta[3]*((1-intprop[j,4])*Lead) + 
           beta[4]*((1-intprop[j,2])*Cadmium) + beta[5]*Copper);
      r1m[j] = mean(r1);
      rd[j] = r1m[j]-r0m;
     } 
   }
   }

