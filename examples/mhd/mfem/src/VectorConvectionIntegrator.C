#include "lmmhd/VectorConvectionIntegrator.hpp"

VectorConvectionIntegrator::VectorConvectionIntegrator(VectorCoefficient &v, real_t a)
   : velocity_coeff(v), alpha(a) { }

void VectorConvectionIntegrator::AssembleElementMatrix(
   const FiniteElement &el, ElementTransformation &Trans, DenseMatrix &elmat)
{
    const int nd   = el.GetDof();
    const int dim  = el.GetDim();
    const int vdim = Trans.GetSpaceDim();

    elmat.SetSize(vdim * nd, vdim * nd);
    elmat = 0.0;

    Vector shape(nd);
    Vector vel(vdim);

    DenseMatrix dshape(nd, dim);

    // Pick quadrature rule.
    const IntegrationRule *ir = &IntRules.Get(el.GetGeomType(), 2 * el.GetOrder());

    // Loop through quadrature points and determine element matrix contributions.
    for (int i = 0; i < ir->GetNPoints(); i++)
    {
        const IntegrationPoint &ip = ir->IntPoint(i);

        Trans.SetIntPoint(&ip);

        // Calculate Jacobian weight to scale contribution.
        const real_t w = ip.weight * Trans.Weight();
        // Evaluate convection field (w in paper).
        velocity_coeff.Eval(vel, Trans, ip);
        // Evaluate shape function values.
        el.CalcShape(ip, shape);
        // Compute derivates of shape functions in physical space.
        el.CalcPhysDShape(Trans, dshape);

        // Loop through DOFs (for test (phi_j)).
        for (int j = 0; j < nd; j++)
        {
            // Loop through DOFs (for trial (phi_k)).
            for (int k = 0; k < nd; k++)
            {
                // Calculate dot product w \cdot \nabla phi_k (basis function).
                real_t dot = 0.0;
                for (int d = 0; d < dim; d++) dot += vel(d) * dshape(k,d);

                // Loop through dimensions.
                for (int vd = 0; vd < vdim; vd++)
                {
                    // Compute integrand for each component (vd - x,y,z).
                    int row = vd * nd + j;
                    int col = vd * nd + k;
                    elmat(row, col) += alpha * dot * shape(j) * w;
                    // Double check which way around the matrix entries need to be.
                }
            }
        }
    }
}