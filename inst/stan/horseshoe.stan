
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
  vector[N] cen_Arsenic;
  vector[N] cen_Manganese;
  vector[N] cen_Lead;
  vector[N] cen_Cadmium;
  vector[N] cen_Copper;
  real mArsenic = mean(Arsenic);
  real mManganese = mean(Manganese);
  real mLead = mean(Lead);
  real mCadmium = mean(Cadmium);
  real mCopper = mean(Copper);
  real sArsenic = sd(Arsenic);
  real sManganese = sd(Manganese);
  real sLead = sd(Lead);
  real sCadmium = sd(Cadmium);
  real sCopper = sd(Copper);
  cen_Arsenic = (Arsenic-mArsenic) ./ sArsenic;
  cen_Manganese = (Manganese-mManganese) ./ sManganese;
  cen_Lead = (Lead-mLead) ./ sLead;
  cen_Cadmium = (Cadmium-mCadmium) ./ sCadmium;
  cen_Copper = (Copper-mCopper) ./ sCopper;
}
parameters{
  vector[p] beta_;
  vector<lower=0>[p] lambda_; 
  real<lower=0> tau_; 
  real<lower=0> zeta; 
  vector<lower=0>[p] nu; 
  real beta0; // given uniform prior
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
  beta0 ~ normal(0, 5);
  beta_ ~ normal(0, 1);
  // horseshoe priors
  nu ~ inv_gamma(1/2., k/2.); // scale mixture prior for student t with 1 d.f., (cauchy(0,k))
  zeta ~ inv_gamma(1/2., k/2.); // scale mixture prior for student t with 1 d.f. (cauchy(0,k))
  tau_ ~ normal(0, 1);
  lambda_ ~ normal(0, 1);

  mu = beta0 + beta[1]*cen_Arsenic + beta[2]*cen_Manganese + beta[3]*cen_Lead +
           beta[4]*cen_Cadmium + beta[5]*cen_Copper + beta[6]*cen_Arsenic .*cen_Arsenic
           + beta[7]*cen_Arsenic .*cen_Manganese + beta[8]*cen_Arsenic .*cen_Lead +
           beta[9]*cen_Arsenic .*cen_Cadmium + beta[10]*cen_Arsenic .*cen_Copper +
           beta[11]*cen_Manganese .*cen_Manganese + beta[12]*cen_Lead .*cen_Manganese
           + beta[13]*cen_Cadmium .*cen_Manganese + beta[14]*cen_Copper .*cen_Manganese
           + beta[15]*cen_Lead .*cen_Lead + beta[16]*cen_Cadmium .*cen_Lead +
           beta[17]*cen_Copper .*cen_Lead + beta[18]*cen_Cadmium .*cen_Cadmium +
           beta[19]*cen_Cadmium .*cen_Copper + beta[20]*cen_Copper .*cen_Copper;
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
  
  r0 = inv_logit(beta0 + beta[1]*((Arsenic-mArsenic) ./ sArsenic) + beta[2]*((Manganese-
           mManganese) ./ sManganese) + beta[3]*((Lead-mLead) ./ sLead) +
           beta[4]*((Cadmium-mCadmium) ./ sCadmium) + beta[5]*((Copper-mCopper) ./ sCopper)
           + beta[6]*((Arsenic-mArsenic) ./ sArsenic) .*((Arsenic-mArsenic) ./ sArsenic)
           + beta[7]*((Arsenic-mArsenic) ./ sArsenic) .*((Manganese-mManganese) ./
           sManganese) + beta[8]*((Arsenic-mArsenic) ./ sArsenic) .*((Lead-mLead) ./
           sLead) + beta[9]*((Arsenic-mArsenic) ./ sArsenic) .*((Cadmium-mCadmium) ./
           sCadmium) + beta[10]*((Arsenic-mArsenic) ./ sArsenic) .*((Copper-mCopper) ./
           sCopper) + beta[11]*((Manganese-mManganese) ./ sManganese) .*((Manganese-
           mManganese) ./ sManganese) + beta[12]*((Lead-mLead) ./ sLead) .*((Manganese-
           mManganese) ./ sManganese) + beta[13]*((Cadmium-mCadmium) ./
           sCadmium) .*((Manganese-mManganese) ./ sManganese) + beta[14]*((Copper-
           mCopper) ./ sCopper) .*((Manganese-mManganese) ./ sManganese) + beta[15]*((Lead-
           mLead) ./ sLead) .*((Lead-mLead) ./ sLead) + beta[16]*((Cadmium-mCadmium) ./
           sCadmium) .*((Lead-mLead) ./ sLead) + beta[17]*((Copper-mCopper) ./
           sCopper) .*((Lead-mLead) ./ sLead) + beta[18]*((Cadmium-mCadmium) ./
           sCadmium) .*((Cadmium-mCadmium) ./ sCadmium) + beta[19]*((Cadmium-mCadmium) ./
           sCadmium) .*((Copper-mCopper) ./ sCopper) + beta[20]*((Copper-mCopper) ./
           sCopper) .*((Copper-mCopper) ./ sCopper));
  
  for(j in 1:5) { //looping over interventions
    r1 = inv_logit(beta0 + beta[1]*((((1-intprop[j,1])*Arsenic)-mArsenic) ./ sArsenic)
           + beta[2]*((((1-intprop[j,2])*Manganese)-mManganese) ./ sManganese) +
           beta[3]*((((1-intprop[j,3])*Lead)-mLead) ./ sLead) + beta[4]*((((1-intprop[j,
           4])*Cadmium)-mCadmium) ./ sCadmium) + beta[5]*((((1-intprop[j,5])*Copper)-
           mCopper) ./ sCopper) + beta[6]*((((1-intprop[j,1])*Arsenic)-mArsenic) ./
           sArsenic) .*((((1-intprop[j,1])*Arsenic)-mArsenic) ./ sArsenic) + beta[7]*((((1-
           intprop[j,1])*Arsenic)-mArsenic) ./ sArsenic) .*((((1-intprop[j,2])*Manganese)-
           mManganese) ./ sManganese) + beta[8]*((((1-intprop[j,1])*Arsenic)-mArsenic) ./
           sArsenic) .*((((1-intprop[j,3])*Lead)-mLead) ./ sLead) + beta[9]*((((1-
           intprop[j,1])*Arsenic)-mArsenic) ./ sArsenic) .*((((1-intprop[j,4])*Cadmium)-
           mCadmium) ./ sCadmium) + beta[10]*((((1-intprop[j,1])*Arsenic)-mArsenic) ./
           sArsenic) .*((((1-intprop[j,5])*Copper)-mCopper) ./ sCopper) + beta[11]*((((1-
           intprop[j,2])*Manganese)-mManganese) ./ sManganese) .*((((1-intprop[j,
           2])*Manganese)-mManganese) ./ sManganese) + beta[12]*((((1-intprop[j,3])*Lead)-
           mLead) ./ sLead) .*((((1-intprop[j,2])*Manganese)-mManganese) ./ sManganese)
           + beta[13]*((((1-intprop[j,4])*Cadmium)-mCadmium) ./ sCadmium) .*((((1-
           intprop[j,2])*Manganese)-mManganese) ./ sManganese) + beta[14]*((((1-
           intprop[j,5])*Copper)-mCopper) ./ sCopper) .*((((1-intprop[j,2])*Manganese)-
           mManganese) ./ sManganese) + beta[15]*((((1-intprop[j,3])*Lead)-mLead) ./
           sLead) .*((((1-intprop[j,3])*Lead)-mLead) ./ sLead) + beta[16]*((((1-intprop[j,
           4])*Cadmium)-mCadmium) ./ sCadmium) .*((((1-intprop[j,3])*Lead)-mLead) ./
           sLead) + beta[17]*((((1-intprop[j,5])*Copper)-mCopper) ./ sCopper) .*((((1-
           intprop[j,3])*Lead)-mLead) ./ sLead) + beta[18]*((((1-intprop[j,4])*Cadmium)-
           mCadmium) ./ sCadmium) .*((((1-intprop[j,4])*Cadmium)-mCadmium) ./ sCadmium) +
           beta[19]*((((1-intprop[j,4])*Cadmium)-mCadmium) ./ sCadmium) .*((((1-intprop[j,
           5])*Copper)-mCopper) ./ sCopper) + beta[20]*((((1-intprop[j,5])*Copper)-
           mCopper) ./ sCopper) .*((((1-intprop[j,5])*Copper)-mCopper) ./ sCopper));
   rd[j] = mean(r1)-mean(r0);
  } 
 
}
}


