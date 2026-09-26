// The background under a reflection: Parkhurst et al. (2016), a robust
// generalized linear model with a Poisson link.
//
// Why not something simpler. On this data the background is 0.3 counts a
// pixel, so most pixels are 0 or 1 and a normal approximation to the Poisson
// is not one. Measured on 1482 reflections of insulin against DIALS' own GLM:
// a five per cent truncated mean is 22.5 per cent low, three sigma clipping
// 14.2 per cent low, and the median is exactly zero because most pixels are
// zero. Every rejection scheme that assumes normality throws away the upper
// tail of a distribution that has nothing else, and a background biased low is
// an intensity biased high for every reflection in the dataset.
//
// THE MODEL
// ---------
// Pixels y_i are Poisson with a common mean mu, and theta = log(mu) is what is
// estimated. With r_i = (y_i - mu)/sqrt(mu) the Pearson residual and the Huber
// weight
//
//     psi_c(r) = r                 |r| <= c
//              = c * sgn(r)        otherwise
//
// the quasi-likelihood score and its Fisher information are, for the constant
// model of Parkhurst Appendix B,
//
//     U = [sum_i psi_c(r_i) - n C1(mu)] sqrt(mu)
//     I = n C2(mu) sqrt(mu)
//     theta <- theta + [sum_i psi_c(r_i) - n C1(mu)] / [n C2(mu)]
//
// with C1(mu) = E[psi_c(r)] and C2(mu) = E[psi_c(r) (y - mu)], the correction
// that makes the estimator Fisher consistent -- without it, clipping one tail
// of an asymmetric distribution biases the answer, which is the whole problem
// being solved.
//
// THE EXPECTATIONS ARE SUMMED, NOT SOLVED
// ---------------------------------------
// Appendix B gives closed forms for C1 and C2 in terms of regularized gamma
// functions, to avoid summing over the distribution. At a background of 0.3
// counts a pixel the distribution is over by j = 20, so the sum is a dozen
// terms and is exact; the closed forms are an optimisation for a regime this
// is not in. They are also the part of the paper whose exponents did not
// survive the PDF's text layer, and transcribing them from a guess would put
// an error somewhere no test would find it.
//
// The sum is taken to `mu + 10 sqrt(mu) + 30`, which holds better than 1e-12
// of a Poisson for any mean these pixels have.

#pragma once

#include <cstddef>
#include <vector>

namespace mxi {

struct BackgroundResult {
  bool valid = false;
  //: The fitted mean, in counts per pixel.
  double mean = 0.0;
  //: Iterations taken; `converged` is false if it ran out.
  int iterations = 0;
  bool converged = false;
  //: How many pixels went in.
  std::size_t n = 0;
};

struct BackgroundOptions {
  //: The Huber tuning constant. 1.345 is 95 per cent efficient on a normal
  //: distribution and is what Parkhurst uses; infinity would be the plain
  //: maximum likelihood estimate, which is the mean.
  double tuning = 1.345;
  int max_iterations = 100;
  //: Convergence on theta, which is log(mu), so this is a relative tolerance
  //: on the background itself.
  double tolerance = 1e-10;
};

//: Fit a constant background to these pixel values.
//:
//: Values are counts and may be zero; a negative value is not a count and the
//: caller should have removed it, along with anything the detector marked bad.
BackgroundResult glm_background(const std::vector<double> &values,
                                const BackgroundOptions &options = {});

//: E[psi_c(r)] and E[psi_c(r) (y - mu)] for a Poisson of mean `mu`, by direct
//: summation. Exposed because they are what the algorithm turns on and a test
//: that cannot see them can only check the answer, not the reason.
double glm_c1(double mu, double tuning);
double glm_c2(double mu, double tuning);

}  // namespace mxi
