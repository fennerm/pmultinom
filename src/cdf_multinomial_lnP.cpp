#include <Rcpp.h>
using namespace Rcpp;
#include <gsl/gsl_cdf.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_randist.h>
#include <math.h>

/* Computes the (log) CDF of a multinomial distribution with parameters p and N
 * i.e. P(X1<=n1,X2<=n2,...,XK<=nK)
 * Uses a approximation from
 * "A Representation for Multinomial Cumulative Distribution Functions",
 * Bruce Levin, The Annals of Statistics, v.9, n.5, pp.1123-1126, 1981
 */
// [[Rcpp::export]]
double cdf_multinomial_lnP(const int K, const int N, Rcpp::NumericVector p,
                           Rcpp::NumericVector n) {
  int k = 0;
  double gamma1 = 0.0, gamma2 = 0.0, sum_s2 = 0.0, sum_mu = 0.0, sp = 0.0,
         x = 0.0, x2 = 0.0, PWN = 0.0;
  double s = static_cast<double>(N);
  double log_cdf = -log(gsl_ran_poisson_pdf(N, s));

  /* This is the P(W=N) bit */
  for (k = 0; k < K; k++) {
    double mf2, mf3, mf4, mr, mu, mu3, mu4, pcdf, s2;

    sp = s * p[k];

    pcdf = gsl_cdf_poisson_P(n[k], sp);
    log_cdf += log(pcdf);

    mu = sp * (1 - gsl_ran_poisson_pdf(n[k], sp) / pcdf);
    s2 = mu - (n[k] - mu) * (sp - mu);

    /* factorial moments? */
    mr = n[k];
    mf2 = sp * mu - mr * (sp - mu);

    mr *= n[k] - 1;
    mf3 = sp * mf2 - mr * (sp - mu);

    mr *= n[k] - 2;
    mf4 = sp * mf3 - mr * (sp - mu);

    /* Central Moments */
    // mu2 = mf2 + mu * (1 - mu);
    mu3 = mf3 + mf2 * (3 - 3 * mu) + mu * (1 + mu * (-3 + 2 * mu));
    mu4 = mf4 + mf3 * (6 - 4 * mu) + mf2 * (7 + mu * (-12 + 6 * mu)) +
          mu * (1 + mu * (-4 + mu * (6 - 3 * mu)));

    /* accumulate coef skewness and excess */
    gamma1 += mu3;
    gamma2 += mu4 - 3 * s2 * s2;
    sum_mu += mu;
    sum_s2 += s2;
  }

  sp = sqrt(sum_s2);
  gamma1 /= sum_s2 * sp;
  gamma2 /= sum_s2 * sum_s2;

  x = (N - sum_mu) / sp;
  x2 = x * x;
  PWN =
      -x2 / 2 +
      log(1 + gamma1 / 6 * x * (x2 - 3) + gamma2 / 24 * (x2 * x2 - 6 * x2 + 3) +
          gamma1 * gamma1 / 72 * (((x2 - 15) * x2 + 45) * x2 - 15)) -
      log(2 * M_PI) / 2 - log(sp);

  log_cdf += PWN;

  return log_cdf;
}
