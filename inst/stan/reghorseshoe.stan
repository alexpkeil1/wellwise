  
// example with regularized horseshoe prior
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
  vector[p] b_;
  vector<lower=0>[p] lambda; 
  real<lower=0> tau; 
  real<lower=0> sig_b; 
  real<lower=0> sig_i; 
  real beta0; // given uniform prior
  real<lower=0> Csq;
}
transformed parameters{
  vector[p] beta;
  {
   vector[p] lambda_t;
   lambda_t = sqrt( Csq * lambda ) ./ sqrt( Csq + tau * lambda .* lambda );
   beta = tau * lambda_t .* b_;
  }
}
model{
{
  vector[N] mu;
  beta0 ~ normal(0, 5);
  sig_b ~ cauchy(0, 1);
  sig_i ~ cauchy(0, .5);
  b_ ~ normal(0, 1);
  // horseshoe priors (not pure if using nu <> 1)
  {
    real nu = 1.;
    real k = 1;
    real s2 = 1;
    tau ~ student_t(1., 0, 1);
    lambda ~ student_t(1., 0, 1);
    // regulazation factor prior
    Csq ~ inv_gamma(k/2, k*s2/2); // translates to student_t with k d.f., scale=sqrt(s2) for coefs far from zero
  }
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
  
  matrix[5,5] intprop = [[.98, .98, .98, .98, 0], // filter 1 reduces metal 1 (Arsenic) by 98%, metal 2 (Cadmium) by 98%, metal 3 (Mn) by 98%, metal 4 (Pb) by 98%, metal 5 (Cu) by 0% 
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

