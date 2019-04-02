
// example with standard logistic model
   data{
    int<lower=0> N;
    int<lower=0> p;
    int<lower=0> dx;
    matrix[N, dx] X; // functions like Nx5 matrix (but each column is real[])
    int<lower=0, upper=1> y[N];
   }
   transformed data{
    vector[N] Arsenic = log(col(X,1));
    vector[N] Manganese = log(col(X,2));
    vector[N] Lead = log(col(X,3));
    vector[N] Cadmium = log(col(X,4));
    vector[N] Copper = log(col(X,5));
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
     
     mu = b0 + beta[1]*Arsenic + beta[2]*Manganese + beta[3]*Lead + beta[4]*Cadmium + beta[5]*Copper +
          beta[6]*Arsenic .* Arsenic + beta[7]*Arsenic .* Manganese + beta[8]*Arsenic .* Lead + beta[9]*Arsenic .* Cadmium +
          beta[10]*Arsenic .* Copper + beta[11]*Manganese .* Manganese + beta[12]*Lead .* Manganese + beta[13]*Cadmium .* Manganese +
          beta[14]*Copper .* Manganese + beta[15]*Lead .* Lead + beta[16]*Cadmium .* Lead +
          beta[17]*Copper .* Lead + beta[18]*Cadmium .* Cadmium + beta[19]*Cadmium .* Copper + beta[20]*Copper .* Copper;
     y ~ bernoulli_logit(mu);
   }
   }
   generated quantities{
   real rd[5];
   {
     vector[N] r1;
     vector[N] r0;
     
   matrix[5,5] intprop = [[.98, .98, .98, .98, 0], // filter 1 reduces metal 1 (Arsenic) by 98%, metal 2 (Cadmium) by 98%, metal 3 (Mn) by 98%, metal 4 (Pb) by 98%, metal 5 (Copper) by 0% 
                          [0, .985, 0, .965, 0],
                          [0, .918, 0, 0, 0],
                          [.994, .98, .98, .987, 0],
                          [0, 0, 0, .99, 0]]; // ROWS = 5 filters, cols = 5 metals
     
     r0 = inv_logit(b0 + beta[1]*Arsenic + beta[2]*Manganese + beta[3]*Lead + beta[4]*Cadmium + beta[5]*Copper +
         beta[6]*Arsenic .* Arsenic + beta[7]*Arsenic .* Manganese + beta[8]*Arsenic .* Lead + beta[9]*Arsenic .* Cadmium +
         beta[10]*Arsenic .* Copper + beta[11]*Manganese .* Manganese + beta[12]*Lead .* Manganese + beta[13]*Cadmium .* Manganese +
         beta[14]*Copper .* Manganese + beta[15]*Lead .* Lead + beta[16]*Cadmium .* Lead +
         beta[17]*Copper .* Lead + beta[18]*Cadmium .* Cadmium + beta[19]*Cadmium .* Copper + beta[20]*Copper .* Copper);
     
     for(j in 1:5) { //looping over interventions
       r1 = inv_logit(b0 + beta[1]*(log(1-intprop[j,1])+Arsenic) + beta[2]*(log(1-intprop[j,3])+Manganese) + beta[3]*(log(1-intprop[j,4])+Lead) + 
            beta[4]*(log(1-intprop[j,2])+Cadmium) + beta[5]*(log(1-intprop[j,5])+Copper) + beta[6]*(log(1-intprop[j,1])+Arsenic) .* (log(1-intprop[j,1])+Arsenic) + 
            beta[7]*(log(1-intprop[j,1])+Arsenic) .* (log(1-intprop[j,3])+Manganese) + beta[8]*(log(1-intprop[j,1])+Arsenic) .* (log(1-intprop[j,4])+Lead) + 
            beta[9]*(log(1-intprop[j,1])+Arsenic) .* (log(1-intprop[j,2])+Cadmium) + beta[10]*(log(1-intprop[j,1])+Arsenic) .* (log(1-intprop[j,5])+Copper) + 
            beta[11]*(log(1-intprop[j,3])+Manganese) .* (log(1-intprop[j,3])+Manganese) + beta[12]*(log(1-intprop[j,4])+Lead) .* (log(1-intprop[j,3])+Manganese) + 
            beta[13]*(log(1-intprop[j,2])+Cadmium) .* (log(1-intprop[j,3])+Manganese) + beta[14]*(log(1-intprop[j,5])+Copper) .* (log(1-intprop[j,3])+Manganese) + 
            beta[15]*(log(1-intprop[j,4])+Lead) .* (log(1-intprop[j,4])+Lead) + beta[16]*(log(1-intprop[j,2])+Cadmium) .* (log(1-intprop[j,4])+Lead) +  
            beta[17]*(log(1-intprop[j,5])+Copper) .* (log(1-intprop[j,4])+Lead) + beta[18]*(log(1-intprop[j,2])+Cadmium) .* (log(1-intprop[j,2])+Cadmium) + 
            beta[19]*(log(1-intprop[j,2])+Cadmium) .* (log(1-intprop[j,5])+Copper) + beta[20]*(log(1-intprop[j,5])+Copper) .* (log(1-intprop[j,5])+Copper));
      rd[j] = mean(r1)-mean(r0);
     } 
   }
   }

