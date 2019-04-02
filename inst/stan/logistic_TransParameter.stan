  
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
vector[p] beta;
real<lower=0> sig_b; 
real<lower=0> sig_i; 
real mu_b;
real b0; // given uniform prior

}
transformed parameters{}
model{
{
  vector[N] mu;
  b0 ~ normal(0, 5);
  sig_b ~ cauchy(0, 1);
  sig_i ~ cauchy(0, 5);
  
  for (i in 1:dx) {
  beta[i] ~ normal(mu_b, sig_b);
  }
  for (i in (dx+1):20) {
  beta[i] ~ normal(mu_b, sig_i);
  // different priors for interactions
  }
  mu = b0 + beta[1]*Arsenic + beta[2]*Manganese + beta[3]*Lead + beta[4]*Cadmium + beta[5]*Copper +
          beta[6]*Arsenic .* Arsenic + beta[7]*Arsenic .* Manganese + beta[8]*Arsenic .* Lead + beta[9]*Arsenic .* Cadmium +
          beta[10]*Arsenic .* Copper + beta[11]*Manganese .* Manganese + beta[12]*Lead .* Manganese + beta[13]*Cadmium .* Manganese +
          beta[14]*Copper .* Manganese + beta[15]*Lead .* Lead + beta[16]*Cadmium .* Lead +
          beta[17]*Copper .* Lead + beta[18]*Cadmium .* Cadmium + beta[19]*Cadmium .* Copper + beta[20]*Copper .* Copper;

  y ~ bernoulli_logit(mu);
}
}
generated quantities{
real rd;
{
  vector[N] r1;
  vector[N] r0;
  r1 =  inv_logit(b0 + beta[1]*1 + beta[2]*Manganese + beta[3]*Lead + beta[4]*Cadmium + beta[5]*Copper +
          beta[6]*1 + beta[7]*Manganese + beta[8]*Lead + beta[9]*Cadmium +
          beta[10]*Copper + beta[11]*Manganese .* Manganese + beta[12]*Lead .* Manganese + beta[13]*Cadmium .* Manganese +
          beta[14]*Copper .* Manganese + beta[15]*Lead .* Lead + beta[16]*Cadmium .* Lead +
          beta[17]*Copper .* Lead + beta[18]*Cadmium .* Cadmium + beta[19]*Cadmium .* Copper + beta[20]*Copper .* Copper);
  r0 =  inv_logit(b0 + beta[1]*0 + beta[2]*Manganese + beta[3]*Lead + beta[4]*Cadmium + beta[5]*Copper +
          beta[6]*0 + beta[7]*0 + beta[8]*0 + beta[9]*0 + beta[10]*0 + 
          beta[11]*Manganese .* Manganese + beta[12]*Lead .* Manganese + beta[13]*Cadmium .* Manganese +
          beta[14]*Copper .* Manganese + beta[15]*Lead .* Lead + beta[16]*Cadmium .* Lead +
          beta[17]*Copper .* Lead + beta[18]*Cadmium .* Cadmium + beta[19]*Cadmium .* Copper + beta[20]*Copper .* Copper);
  rd = mean(r1)-mean(r0);
}
}

