#include "mfem.hpp"

using namespace std;
using namespace mfem;

inline real_t velocity_nbc(const Vector & x)
{
   return 0.0;
}

inline real_t pressure_nbc(const Vector & x)
{
   return 0.0;
}


inline real_t pressure_dbc(const Vector & x)
{
   if (x(0) > 0.0)
   { // Value at outlet.
      return 1.0;
   }
   else
   { // Value at inlet.
      return 1.0;
   }
}

inline real_t zero_dbc(const Vector & x)
{
   return 0.0;
}

inline real_t outflow_term_func(const Vector & x)
{
   // Half because of eq. (2.2a) and definition of O.
   return 0.5;
}


inline void velocity_dbc_vec_func(const Vector & x, Vector & f)
{
   f = 0.0;

   real_t y_mid = 1.0;
   real_t z_mid = 1.0;
   real_t u_avg = 1.0;
   
   if (x(0) < 1e-6)
   { // Parabolic condition in x-direction at inlet.
      //f(0) = (9. / 4.) * u_avg * (1. - ((x(1) - y_mid)*(x(1) - y_mid)) / (y_mid * y_mid)) * 
      //               (1. - ((x(2) - z_mid)*(x(2) - z_mid)) / (z_mid * z_mid));
      f(0) = 0.0;
      f(1) = 0.0;
   }
   else
   { // Zero on top and bottom boundaries.
      f(0) = 0.0;
      f(1) = 0.0;
   }

   // Lid-driven cavity.
   if (x(2) > 0.999)
   {
      f(0) = 1.0;
   }

   //if (x.Size() == 3)
   //{
   //   f(2) = 0.0;
   //}
}

inline void currentD_dbc_vec_func(const Vector & x, Vector & f)
{
   f = 0.0;
}

inline real_t electPot_dbc(const Vector & x)
{
   return 0.0;
}

inline void zero_func(const Vector & x, Vector & f)
{
   f = 0.0;
}


inline void u_exact(const mfem::Vector & x, mfem::Vector & f)
{

   real_t r_max = 0.1;
   real_t u_avg = 1.0;

   f(0) = u_avg * (1. - ((x(1)-r_max)*(x(1)-r_max)) / (r_max*r_max));
   f(1) = 0.0;
   if (x.Size() == 3)
   {
      f(2) = 0.0;
   }
   
}
