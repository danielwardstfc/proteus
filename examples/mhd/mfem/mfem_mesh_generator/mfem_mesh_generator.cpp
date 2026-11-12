#include "mfem.hpp"
#include <memory>
#include <iostream>
#include <fstream>
#include <cmath>

using namespace std;
using namespace mfem;

mfem::real_t cluster_symmetric(mfem::real_t xi, mfem::real_t factor)
{
    // xi in [0,1], cluster away from 0.5
    mfem::real_t x_shifted = xi - 0.5;
    return 0.5 + 0.5 * tanh(factor * x_shifted) / tanh(factor / 2.0);
}

int main(int argc, char *argv[])
{
   // Mesh parameters
   int nx = 40;
   int ny = 30;
   int nz = 30;

   real_t Lx = 4.0;
   real_t Ly = 1.0;
   real_t Lz = 1.0;

   real_t clusterX = 1.0;
   real_t clusterY = 1.8;
   real_t clusterZ = 1.8;

   // Generate mesh.
   Mesh mesh = Mesh::MakeCartesian3D(nx, ny, nz, mfem::Element::Type::HEXAHEDRON, Lx, Ly, Lz);
   int dim = mesh.Dimension();

   // Cluster vertices towards the walls in x, y and z.
   int numVertices = mesh.GetNV();
   for (int i = 0; i < numVertices; i++)
   {
      real_t *v = mesh.GetVertex(i);
      real_t xi = v[0] / Lx; // Normalize x
      real_t yi = v[1] / Ly; // Normalize y
      real_t zi = v[2] / Lz; // Normalize z

      v[0] = cluster_symmetric(xi, clusterX) * Lx;
      v[1] = cluster_symmetric(yi, clusterY) * Ly;
      v[2] = cluster_symmetric(zi, clusterZ) * Lz;
   }

   mesh.Save("box_clustered.mesh"); // Save mesh to file.

   mesh.Clear();

   return 0;
}
