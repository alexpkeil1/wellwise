  
// example with horseshoe prior
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
  vector[p] beta_;
  vector<lower=0>[p] lambda_; 
  real<lower=0> tau_; 
  real<lower=0> zeta; 
  vector<lower=0>[p] nu; 
  real b0; // given uniform prior
  //real<lower=0> k; // goes wild
}
transformed parameters{

  vector[p] beta;
  vector<lower=0>[p] lambda; 
  real<lower=0> tau; 
  // tau has t-dist with k d.f. if its a inv-gamma mixture of normals
  lambda = lambda_ .* sqrt(nu);
  tau = tau_ * sqrt(zeta);
  beta = tau * lambda .* beta_;
}
model{
{
  vector[N] mu;
  real k = 1; // d.f. for t distribution of global shrinkage parameter 
                //(can make into a r.v., as well, but that is very slow)
  b0 ~ normal(0, 5);
  beta_ ~ normal(0, 1);
  // horseshoe priors
  nu ~ inv_gamma(1/2., k/2.); // scale mixture prior for student t with 1 d.f., (cauchy(0,k))
  zeta ~ inv_gamma(1/2., k/2.); // scale mixture prior for student t with 1 d.f. (cauchy(0,k))
  tau_ ~ normal(0, 1);
  lambda_ ~ normal(0, 1);

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
  
  matrix[5,5] intprop = [[.98, .98, .98, .98, 0], // filter 1 reduces metal 1 (Arsenic) by 98%, metal 2 (Cadmium) by 98%, metal 3 (Mn) by 98%, metal 4 (Pb) by 98% 
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
       r1 = inv_logit(b0 + beta[1]*((1-intprop[j,1])*Arsenic) + beta[2]*((1-intprop[j,3])*Manganese) + beta[3]*((1-intprop[j,4])*Lead) + 
            beta[4]*((1-intprop[j,2])*Cadmium) + beta[5]*((1-intprop[j,5])*Copper) + beta[6]*((1-intprop[j,1])*Arsenic) .* ((1-intprop[j,1])*Arsenic) +
            beta[7]*((1-intprop[j,1])*Arsenic) .* ((1-intprop[j,3])*Manganese) + beta[8]*((1-intprop[j,1])*Arsenic) .* ((1-intprop[j,4])*Lead) + 
            beta[9]*((1-intprop[j,1])*Arsenic) .* ((1-intprop[j,2])*Cadmium) + beta[10]*((1-intprop[j,1])*Arsenic) .* ((1-intprop[j,5])*Copper) + 
            beta[11]*((1-intprop[j,3])*Manganese) .* ((1-intprop[j,3])*Manganese) + beta[12]*((1-intprop[j,4])*Lead) .* ((1-intprop[j,3])*Manganese) + 
            beta[13]*((1-intprop[j,2])*Cadmium) .* ((1-intprop[j,3])*Manganese) + beta[14]*((1-intprop[j,5])*Copper) .* ((1-intprop[j,3])*Manganese) + 
            beta[15]*((1-intprop[j,4])*Lead) .* ((1-intprop[j,4])*Lead) + beta[16]*((1-intprop[j,2])*Cadmium) .* ((1-intprop[j,4])*Lead) + 
            beta[17]*((1-intprop[j,5])*Copper) .* ((1-intprop[j,4])*Lead) + beta[18]*((1-intprop[j,2])*Cadmium) .* ((1-intprop[j,2])*Cadmium) + 
            beta[19]*((1-intprop[j,2])*Cadmium) .* ((1-intprop[j,5])*Copper) + beta[20]*((1-intprop[j,5])*Copper) .* ((1-intprop[j,5])*Copper));
      rd[j] = mean(r1)-mean(r0);
     } 
}
}

