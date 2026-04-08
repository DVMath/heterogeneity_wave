functions {
  array[] real SIR(real t, array[] real y, array[] real theta,
  data array[] real x_r, data array[] int x_i) {
    real S     = y[1];
    real I     = y[2];
    real R     = y[3];
    real D     = y[4];
    real beta  = theta[1];
    real gamma = theta[2];
    real omega = theta[3];
    real n1 = x_r[1];
    real nt1 = S+I+R;
    int inc = x_i[1];
    array[5] real dydt;
    dydt[1] = -beta * S * I / (nt1);
    dydt[2] =  beta * S * I / (nt1) - gamma * I - omega * I;
    dydt[3] =  gamma * I;
    dydt[4] =  omega * I;
    dydt[5] =  beta * S * I/(nt1);
    return dydt;
  }
}
data {
  /* version with cumulative viral titer */
  int N_ind;
  int N_obs;
  int N_contact;
  array[N_obs] real y_obs;
  array[N_obs] int ind;
  array[N_obs] int t;
  array[N_ind] int ysev;
  array[N_ind] real peak;
  array[N_ind] real ratevir;
  array[N_ind] real area;
  array[N_ind] int male;
  array[N_ind] int over15;
  array[N_ind] int serotype;
  array[N_ind] int tevent;
  array[N_contact] int inf_contact;
  array[N_contact] real day_VL;
  int nserotype;
 // real mean_t;
//  real sd_t;
  real eps0;
  real eps1;
  int myinf;
  real d_fixed;
  real sigma_meas;
  int nresp;
  /* population data */
    /* population data */
  int<lower=1> n_t;   // number of observation times
  real t0;              // initial time (e.g. 0)
  array[n_t] real ts; // times at which y is evaluated
  int n_tcobs;
  array[n_tcobs] int ncases;   // observed infected counts
  array[n_tcobs] int ndeaths;   // observed infected counts
  array[n_tcobs] int tc_obs; // times at which y is evaluated
  int<lower=1> npop;   // total population size
  real popnorm;
  real Y0;
  real CFRrough;
  real betarough;
}
transformed data {
  array[1] real x_r = {popnorm};  // Empty array of reals
  array[1] int x_i ={1};   // Empty array of integers
  real gamma_fixed = 1/6.5;  // From biology
  //  CFR = mu /(gamma+mu)
  //  CFR*gamma + CFR*mu = mu
  //  mu = CFR*gamma/(1-CFR)
  real omega_fixed = 0.087; //CFRrough*gamma_fixed/(1-CFRrough) ;
//  real phi_cases_fixed = 1.5;  // From previous fit
//  real phi_deaths_fixed = 2.5;  // From previous fit
  real phi_cases_fixed = 20;  // From previous fit
  real phi_deaths_fixed = 30;  // From previous fit
  real  Y0_fixed = .001;
}
parameters {
    real beta_raw;
  real<upper=0> omega_raw;  // CFR ~ 0.5-2%
  real<upper=0> gamma_raw;
  real R0_raw;
  real<lower=0.001, upper=.05> popfactor;
   real<lower=0> phi;  // Overdispersion parameter
   real<lower=0> phid;  // Overdispersion parameter
 // real Y0p;
  real alpha;
  real beta;
  // real tau;
  real c;
  //real d;
  real psisex;
  real psiage;
  real s2;
  real s3;
  real s4;
  array[N_ind] real betai;
  array[N_ind] real alphai;
  array[N_ind] real lf;
 // real sigma;
  real gamma0;
  real gammaa;
  real gammab;
  real gammap;
//  real gammar;
  real gammar_l;
  real betat;
  real trans0;
  real trans1l;
   real c_mu;
  real<lower=0> c_sigma;
  vector[N_ind] c_raw;          // non-centered
}
transformed parameters {
 // array[N_ind] real<lower =0> tau;
  array[N_obs] real yt;
  array[N_ind] real betas;
  array[N_ind] real alphas;
  array[N_ind] real f;
  array[N_ind] real EIP;
  array[N_ind] real maxi;
  array[N_ind] real psev;
  array[N_ind] real csev;
  array[N_contact] real transv;
  array[nserotype] real fserotype;
  array[N_ind] real ci;
 // real alpha;
  real a;
  real b;
  real alphabeta;
  real e_alphabeta;
  array[N_ind, nresp] real vt;           // Cumulative viral load (or similar)
  array[N_ind, nresp] real sumvt;           // Cumulative viral load (or similar)
  array[N_ind, nresp] real sumlvt;           // Cumulative viral load (or similar)
  array[N_ind, nresp] real lpsurv;          // Log probability of survival
  array[N_ind, nresp] real lptrans;         // Log probability of transmission
  array[N_ind, nresp] real rv;              // Combined risk
  vector[N_ind] sumrv;                   // Sum of rv across time
  vector[N_ind] sumtrans;                // Sum of transmission probabilities
  real sumtrans_avg;
  real sumrv_avg;
  real constant_beta;
  real constant_R0;
  array[5] real y_init;             // [S(0), I(0), R(0)]
  array[n_t,5] real y_hat; // solution at each ts
  array[n_tcobs] real lambda;     // Poisson rates
  array[n_tcobs] real lambdad;     // Poisson rates
  array[3] real<lower=0> theta; // {beta, gamma}
  real<lower=0,upper=1> gamma_trans;
  real<lower=0,upper=1> omega_trans;
  real<lower=0> beta_trans;
 // real constant_beta;
  real<lower=0, upper=5> R0_trans;
  real popnorm_factor;
  real Y0p;
  real Y0f;
  popnorm_factor = npop/popnorm;
    omega_trans = exp(omega_raw);
    gamma_trans = exp(gamma_raw);
  // beta_trans = betarough+ beta_raw;
   R0_trans = 1 + exp(R0_raw);
   beta_trans = R0_trans*(gamma_trans + omega_trans);
//   R0_trans = beta_trans/(gamma_trans+omega_trans);
  Y0f = Y0;
  Y0p = Y0f;
  y_init[1] = popfactor*popnorm-Y0p;
  y_init[2] = Y0p;
  y_init[3] = 0;
  y_init[4] = 0;
  y_init[5] = Y0p;
  for (i in 1:N_ind) {
    ci[i] = c_mu + 
           c_raw[i] * c_sigma;
  //  ci[i] = c_raw[i];
  }
  // y_hat = ode_rk45(SIR, y_init, t0, ts, theta);
 for (lj in 1:N_contact) {
  transv[lj] = trans0 + trans1l * day_VL[lj]; 
 }
 fserotype[1] = 0;
 if (nserotype>1) {
  fserotype[2] = s2;
  fserotype[3] = s3;
  fserotype[4] = s4;
 }
  a= exp(alpha);
  b = exp(beta);
  alphabeta = alpha+beta;
  e_alphabeta = exp(alphabeta);
  for (i in 1:N_ind) {
  // betai[i] = prodi[i] - alphai[i];
  f[i] = exp(lf[i])/(1+exp(lf[i]));
  betas[i] = exp(beta+betai[i] + psisex*male[i] + psiage*over15[i] + fserotype[serotype[i]]);
  alphas[i] = f[i]*betas[i];
  EIP[i] = (log(alphas[i]/(betas[i]-alphas[i])) + ci[i])/betas[i];
  maxi[i] = d_fixed + alphas[i]*EIP[i] - log1p_exp(betas[i]*EIP[i] - ci[i]);
 // opcoes ratevir (indireta) ou alphas (transmissibilidade) 
 // psev[i] = 1.0/(1.0+exp(-(gamma0+ (gammar + gammaa*male[i] +gammab*over15[i])*log(area[i]))));
// csev[i] = gamma0 + gammaa*male[i] +gammab*over15[i] + gammar*log(area[i]);
 // psev[i] = inv_logit(csev[i]);
//  psev[i] = 1/(1+exp(-(gamma0+ gammar*alphas[i] +gammaa*male[i] +gammab*over15[i] )));
}
  /* using the cummulative viral titer */
  for (ii in 1:N_obs) {
 yt[ii] = d_fixed + alphas[ind[ii]]*(t[ii]) - log1p_exp(betas[ind[ii]]*(t[ii]) - ci[ind[ii]] );
}
 for (jj in 1:N_ind) {
   real logS=0.0;
    real gamma_linear = gamma0 + gammaa*male[jj] + gammab*over15[jj];
    // Time 1
    vt[jj, 1] = d_fixed + alphas[jj] - log1p_exp(betas[jj] - ci[jj]);
    sumlvt[jj,1] = fmax(vt[jj,1], 0.0);
    sumvt[jj, 1] = exp(vt[jj, 1]); // d_fixed + alphas[jj] - log1p_exp(betas[jj] - ci[jj]);
    if (1==1) {
      lpsurv[jj, 1] = log1m_inv_logit(gamma_linear + betat* 1 + gammar_l * maxi[jj]);
    lptrans[jj, 1] = log_inv_logit(trans0 + trans1l * vt[jj, 1]);
    rv[jj, 1] = exp(lptrans[jj, 1] + logS);
    logS += log1m_exp(lpsurv[jj,1]);
    // Times 2:30
    for (ll in 2:nresp) {
      vt[jj, ll] = d_fixed + alphas[jj]*ll - log1p_exp(betas[jj]*ll - ci[jj]);
      sumlvt[jj, ll] = sumlvt[jj,ll-1] + fmax(vt[jj, ll], 0.0);
      sumvt[jj, ll] = sumvt[jj, ll-1] + exp(vt[jj, ll]);
      if (is_nan(sumvt[jj,ll]) || is_inf(sumvt[jj,ll])) reject("bad sumvt");
     // print("vt ", jj, " ll ", ll, " value ", vt[jj,ll], " sumvt ", sumvt[jj,ll]);
      lpsurv[jj, ll] = log1m_inv_logit(gamma_linear + betat * ll + gammar_l * maxi[jj]);
      lptrans[jj, ll] = log_inv_logit(trans0 + trans1l * vt[jj, ll]);
      rv[jj, ll] = exp(lptrans[jj, ll] + logS);
      logS += log1m_exp(lpsurv[jj,ll]);
    }
    sumrv[jj] = sum(rv[jj, ]);
    sumtrans[jj] = sum(exp(lptrans[jj, ]));
    }
  }
  sumtrans_avg = sum(sumtrans)/N_ind;
  sumrv_avg = sum(sumrv)/N_ind;
  for (i in 1:N_ind) {
// csev[i] = gamma0 + gammaa*male[i] +gammab*over15[i] + gammar*log(area[i]);
// csev[i] = gamma0 + gammaa*male[i] +gammab*over15[i] + betat * nresp + gammar*log(sumvt[i,nresp]);
// csev[i] = gamma0 + gammaa*male[i] +gammab*over15[i] + betat * nresp + gammar*vt[i,nresp];
 csev[i] = gamma0 + gammaa*male[i] +gammab*over15[i] + betat * nresp + gammar_l*maxi[i];
 psev[i] = inv_logit(csev[i]);
//  psev[i] = 1/(1+exp(-(gamma0+ gammar*alphas[i] +gammaa*male[i] +gammab*over15[i] )));
}
 theta[1] = beta_trans;
//  gamma_trans = beta_trans/R0_trans - omega_trans;
// theta[1] = beta_trans;
theta[2] = gamma_trans;
theta[3] = omega_trans;
constant_beta = beta_trans/sumtrans_avg;
constant_R0 = R0_trans/sumrv_avg;
  // print("y_init: ", y_init);  // Debug print 
//  print("theta: ", theta);  // Debug print  // print("y_init: ", y_init);  // Debug print 
 // print("theta: ", theta);  // Debug print  // print("y_init: ", y_init);  // Debug print 
 // print("theta: ", theta);  // Debug print
 y_hat = integrate_ode_rk45(SIR, y_init, t0, ts, theta, x_r, x_i);
 // Find when peak occurs
  real max_I = 0;
  real peak_day = 0;
  for (i in 1:n_t) {
    if (y_hat[i, 2] > max_I) {
      max_I = y_hat[i, 2];
      peak_day = ts[i];  // Actual time (e.g., day 25)
    }
  }
   //   print("y_hat 5: ", y_hat[,5]);  // Debug print
  if (peak_day > 50) {
    reject("Peak too late");
  }
  if (peak_day <= 1) {
    reject("no peak"); 
  }
  lambda[1] = log(fmax((y_hat[tc_obs[1], 5]-y_hat[1,5])*popnorm_factor, 1e-6));
  lambdad[1] = log(fmax((y_hat[tc_obs[1], 4]-y_hat[1,4])*popnorm_factor, 1e-6));
  for (i in 2:n_tcobs) {
    lambda[i] = log(fmax((y_hat[tc_obs[i], 5] - y_hat[tc_obs[i-1],5])*popnorm_factor, 1e-6));
    lambdad[i] = log(fmax((y_hat[tc_obs[i], 4] - y_hat[tc_obs[i-1], 4])*popnorm_factor, 1e-6));
   // sigmac[i] = phi + phiscale*(lambda[i])^0.5;
  //  sigmad[i] = phid + phiscaled*(lambdad[i])^0.5;
  }
 //  print("lambda ",  exp(lambda));
}
/* model */
model {
/*priors */
popfactor ~ normal(0,1);
R0_raw ~ normal(0,1);
 beta_raw ~ normal(0, 1);
  omega_raw ~ normal(log(0.01), 1);  // CFR ~ 0.5-2%
  gamma_raw ~ normal(log(1/6.5), 1);
 phi ~ exponential(.1);  // Weakly informative
 phid ~ exponential(.1);  // Weakly informative
// p ~ normal(0, 1);
// q ~ normal(0,1);
// pi ~ normal(0, 1);
// qi ~ normal(0,1);
/* alpha, beta */
alpha ~ normal(0,1);
//alphai0 ~ normal(0,1);
//alphai1 ~ normal(0,1);
//za ~ normal(0,1);
 phi ~ exponential(0.1);  // Weakly informative
 phid ~ exponential(0.1);  // Weakly informative
trans0 ~normal(0,10);
trans1l ~normal(0,10);
beta ~ normal(0,1);
// Y0f ~ normal(0, 0.1);
// tau ~ normal(0,10);
c ~ normal(0,1);
 c_mu ~ normal(30, 1);  // rough scale
  c_sigma ~ exponential(1);
  c_raw ~ std_normal();
// d ~ normal(0,1);
/* tau */
for (j in 1:N_ind) {
  //ltau[j] ~ normal(p, q);
  betai[j] ~ normal(0, 1);
//   prodi[j] ~ normal(0, 1);
   alphai[j] ~ normal(0, 1);
   lf[j] ~ normal(0,1);
 //  ysev[j] ~ bernoulli_logit(csev[j]);
}
inf_contact ~ bernoulli_logit(transv);
/* sigma */
// sigma ~  normal(0,1);
gamma0 ~ normal(0,10);
gammaa ~ normal(0,10);
gammab ~ normal(0,10);
gammap ~ normal(0,10);
gammar_l ~ normal(0,10);
betat ~ normal(0,1);
//l_beta_trans ~ normal(0,1);
//l_omega_trans ~ normal(0,10);
//l_R0_trans_add ~ normal(0,1);
// gamma_trans ~ normal(0,10);
psiage ~ normal(0,1);
psisex ~ normal(0,1);
s2 ~ normal(0,1);
s3 ~ normal(0,1);
s4 ~ normal(0,1);
/* observations */
for (ii in 1:N_obs) {
  y_obs[ii] ~ normal(yt[ii], sigma_meas);
}
/* instead of psev for now */
//  csev[i] = gamma0 + gammaa*male[i] +gammab*over15[i] + gammar*log(sumvt[i,nresp]);
  for (i in 1:N_ind)
    for (tt in 1:tevent[i])
      (ysev[i] == 1 && tt == tevent[i] ? 1 : 0) 
        ~ bernoulli_logit(gamma0 + gammaa*male[i] + gammab*over15[i]
                          + betat*tt + gammar_l*maxi[i]);
/* observations */
for (ii in 1:n_tcobs) {
//  ncases[ii] ~ poisson(lambda[ii]);
//  ndeaths[ii] ~ poisson(lambdad[ii]);
  ncases[ii] ~ neg_binomial_2(exp(lambda[ii]), phi);
  ndeaths[ii] ~ neg_binomial_2(exp(lambdad[ii]), phid);
}
}
generated quantities {
  array[n_t,5] real y_pred; // solution at each ts
  vector[n_tcobs] lambda_pred;  // Predicted case counts
  vector[n_tcobs] lambdad_pred;  // Predicted case counts
//  array[n_t] vector[4] y_rep; // solution at each ts
  array[n_tcobs] int y_rep;  // Posterior predictive samples
  array[n_tcobs] int yd_rep;  // Posterior predictive samples
  array[N_ind] real R0_ind;
  array[N_ind] real beta_ind;
  array[N_ind] real gamma_ind;
  array[N_ind] real omega_ind;
  array[N_ind] real pop_ind;
  array[N_ind] real R0_estim;
  array[3] real theta_pred;
  real trans1 = log(10)*trans1l;
  real gammar = log(10)*gammar_l;
  /* within-host quantities */
/* population level evaluation */
  theta_pred[1] = beta_trans;
  theta_pred[2] = gamma_trans;
  theta_pred[3] = omega_trans;
  // Solve ODE with sampled parameters
//  y_pred = ode_rk45(SIR, y_init, t0, ts, theta_pred, 1e-8, 1e-6, 1e6);
//  y_pred = integrate_ode_rk45(SIR, y_init, t0, ts, theta_pred, x_r, x_i, 1e-6, 1e-6, 1e5);
  y_pred = integrate_ode_rk45(SIR, y_init, t0, ts, theta_pred, x_r, x_i);
  lambda_pred[1] = log(fmax((y_pred[tc_obs[1], 5]-y_pred[1,5])*popnorm_factor, 1e-6));
  lambdad_pred[1] = log(fmax((y_pred[tc_obs[1], 4]-y_pred[1,4])*popnorm_factor, 1e-6));
  y_rep[1] = neg_binomial_2_rng(exp(lambda_pred[1]), phi);
  yd_rep[1] = neg_binomial_2_rng(exp(lambdad_pred[1]), phid);
  // Generate predictions
  for (i in 2:n_tcobs) {
    lambda_pred[i] = log(fmax((y_pred[tc_obs[i], 5] - y_pred[tc_obs[i-1],5])*popnorm_factor, 1e-6));
    lambdad_pred[i] = log(fmax((y_pred[tc_obs[i], 4] - y_pred[tc_obs[i-1], 4])*popnorm_factor, 1e-6));
    y_rep[i] = neg_binomial_2_rng(exp(lambda_pred[i]), phi);  // Simulate new data
    yd_rep[i] = neg_binomial_2_rng(exp(lambdad_pred[i]), phid);  // Simulate new data
  }
  for (ii in 1:N_ind) {
/* changed to output*/
//  beta_ind[ii] = constant_R0 * sumtrans[ii]; // beta_trans;
  gamma_ind[ii] = gamma_trans;
  /* psev   */
  omega_ind[ii] = psev[ii]*gamma_trans/(1-psev[ii]);
  pop_ind[ii] =  popnorm*popfactor;
  R0_estim[ii] = constant_R0*sumrv[ii];
  R0_ind[ii] = constant_beta * sumtrans[ii]/(gamma_trans + psev[ii]*gamma_trans/(1-psev[ii]));
  beta_ind[ii] = R0_estim[ii]*(gamma_ind[ii]+omega_ind[ii]); // beta_trans;
}
}
