#ifndef VECTOR_CONVECTION_INTEGRATOR_HPP
#define VECTOR_CONVECTION_INTEGRATOR_HPP

#include "mfem.hpp"

using namespace mfem;

// $\alpha (Q \cdot \nabla u, v)$
class VectorConvectionIntegrator : public BilinearFormIntegrator
{
private:
   VectorCoefficient &velocity_coeff;  // Rename to vel_coeff or similar to be clear it's a coefficient.
   real_t alpha;

public:
   VectorConvectionIntegrator(VectorCoefficient &v, real_t a = 1.0);

   virtual void AssembleElementMatrix(const FiniteElement &el,
                                      ElementTransformation &Trans,
                                      mfem::DenseMatrix &elmat) override;
};


// $-\alpha (v, Q \cdot \nabla v)$, negative transpose of VectorConvectionIntegrator
class ConservativeVectorConvectionIntegrator : public TransposeIntegrator
{
public:
   ConservativeVectorConvectionIntegrator(VectorCoefficient &q, real_t a = 1.0)
      : TransposeIntegrator(new VectorConvectionIntegrator(q, -a)) { }
};

#endif // VECTOR_CONVECTION_INTEGRATOR_HPP