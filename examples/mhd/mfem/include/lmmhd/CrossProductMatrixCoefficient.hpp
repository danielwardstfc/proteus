class CrossProductMatrixCoefficient : public MatrixCoefficient
{
private:
   Vector B;

public:
   CrossProductMatrixCoefficient(const Vector &Bvec)
      : MatrixCoefficient(Bvec.Size()), B(Bvec) { }

   virtual void Eval(mfem::DenseMatrix &M, ElementTransformation &T,
                     const IntegrationPoint &ip) override
   {
      M.SetSize(3);
      real_t bx = B(0), by = B(1), bz = B(2);

      M(0,0) = 0.0;   M(0,1) = -bz;   M(0,2) =  by;
      M(1,0) =  bz;   M(1,1) = 0.0;   M(1,2) = -bx;
      M(2,0) = -by;   M(2,1) =  bx;   M(2,2) = 0.0;
   }
};