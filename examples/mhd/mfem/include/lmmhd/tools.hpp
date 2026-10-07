#pragma once
#include "mfem.hpp"
#include <iostream>
#include <fstream>

using namespace mfem;
using namespace std;

// Print all parameter read in from input file.
void PrintParams(int nx, int ny, int nz, 
		real_t Lx, real_t Ly, real_t Lz,
                real_t clusterX, real_t clusterY, real_t clusterZ,
                real_t t_final, real_t dt, int vis_steps, 
		real_t Re, real_t Ha);

// Print finite element space sizes.
void PrintFESpaces(int j_space_size, int phi_space_size, 
        int v_space_size, int p_space_size);

// Clustering function (symmetric around center).
real_t cluster_symmetric(real_t xi, real_t factor);

// Checkpoint function for debugging.
void checkpoint(int num);

// Calculate relative L2-norm of grid functions.
//real_t rel_L2_norm(Vector X, Vector Xn_1, ParFiniteElementSpace *fespace);
real_t rel_L2_norm(const ParGridFunction &gf, const ParGridFunction &gfN_1, ParFiniteElementSpace *fespace);

class RemoveMeanProjector
{
private:

   ParLinearForm *mass_lf = nullptr;
   ConstantCoefficient onecoeff;
   Vector one_vec, mass_vec;
   real_t volume = 0.0;

public:

   RemoveMeanProjector() = default;

   RemoveMeanProjector(ParFiniteElementSpace &fes);

   void RemoveMean(Vector &v) const;
};

// Helper class for outputting to log file.
class Logger
{
private:
    std::ofstream &file;

public:
    Logger(std::ofstream &file_) : file(file_) {}

    template <typename T>
    Logger &operator<<(const T &value)
    {
        if (mfem::Mpi::Root())
        {
            file << value;
        }
        return *this;
    }

    // Support std::endl
    Logger &operator<<(std::ostream& (*manip)(std::ostream&))
    {
        if (mfem::Mpi::Root())
        {
            manip(file);
        }
        return *this;
    }
};

// Helper class for accessing FGMRES residuals and outputting to log file.
class FGMRESLogMonitor : public mfem::IterativeSolverMonitor
{
private:
    Logger logger;

public:
    FGMRESLogMonitor(Logger &logger_)
        : logger(logger_)
    {}

    virtual void MonitorResidual(int it,
                                 mfem::real_t norm,
                                 const mfem::Vector &r,
                                 bool final) override
    {
        if (mfem::Mpi::Root())
        {
            logger << "FGMRES iteration : "
                     << it
                     << " || r || = "
                     << norm
                     << std::endl;
        }
    }
};