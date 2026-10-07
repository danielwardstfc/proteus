// constants.hpp

#include "mfem.hpp"
#include "lmmhd/constants.hpp"

using namespace std;
using namespace mfem;

// Define constants.
real_t Re(1.0);     // Reynolds number.
real_t Ha(1.0);    // Hartmann number.

real_t reciprocal_Re(1 / Re);
real_t alpha(1.0); // alpha = 1 (default).
real_t alpha1(alpha + reciprocal_Re);
real_t neg_alpha1(-alpha1);
real_t kappa_val(Ha * Ha / Re);     // Stuart number.
real_t neg_kappa_val(-kappa_val);

real_t Bx(0.0);
real_t By(0.0);
real_t Bz(1.0);