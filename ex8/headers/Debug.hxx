#pragma once

#include <cmath>
#include <limits>
#include <iostream>
#include <string>

#include "MFEMMGIS/NonLinearEvolutionProblem.hxx"
#include "MFEMMGIS/Material.hxx"

#ifdef MFEM_USE_MPI
#include <mpi.h>
#endif

inline void debug_print_physics_stats(
    mfem_mgis::Context& ctx, 
    mfem_mgis::NonLinearEvolutionProblem& heat_transfer,
    mfem_mgis::NonLinearEvolutionProblem& mechanics,
    bool parallel,
    const SetupPropertiesResult& setup) 
{
    // Helper lambda for nodal fields
    auto print_nodal_stats = [&](const mfem::GridFunction& gf, const std::string& name, bool is_vector) {
        const mfem::FiniteElementSpace* fes = gf.FESpace();
        int ndofs = fes->GetNDofs();
        int vdim = fes->GetVDim();

        double local_min = std::numeric_limits<double>::max();
        double local_max = -std::numeric_limits<double>::max();
        double local_sum = 0.0;
        double local_sum_sq = 0.0;
        long long local_count = ndofs;

        for (int i = 0; i < ndofs; ++i) {
            double val = 0.0;
            if (is_vector) {
                double mag_sq = 0.0;
                for (int d = 0; d < vdim; ++d) {
                    double comp = gf(fes->DofToVDof(i, d));
                    mag_sq += comp * comp;
                }
                val = std::sqrt(mag_sq);
            } else {
                val = gf(fes->DofToVDof(i, 0));
            }

            if (val < local_min) local_min = val;
            if (val > local_max) local_max = val;
            local_sum += val;
            local_sum_sq += val * val;
        }

        double g_min = local_min, g_max = local_max, g_sum = local_sum, g_sum_sq = local_sum_sq;
        long long g_count = local_count;

        #ifdef MFEM_USE_MPI
        if (parallel) {
            MPI_Comm comm = MPI_COMM_WORLD;
            if (auto pfes = dynamic_cast<const mfem::ParFiniteElementSpace*>(fes)) {
                comm = pfes->GetComm();
            }
            
            MPI_Allreduce(&local_min, &g_min, 1, MPI_DOUBLE, MPI_MIN, comm);
            MPI_Allreduce(&local_max, &g_max, 1, MPI_DOUBLE, MPI_MAX, comm);
            MPI_Allreduce(&local_sum, &g_sum, 1, MPI_DOUBLE, MPI_SUM, comm);
            MPI_Allreduce(&local_sum_sq, &g_sum_sq, 1, MPI_DOUBLE, MPI_SUM, comm);
            MPI_Allreduce(&local_count, &g_count, 1, MPI_LONG_LONG, MPI_SUM, comm);
        }
        #endif

        if (mfem_mgis::getMPIrank() == 0 && g_count > 0) {
            double mean = g_sum / g_count;
            double var = (g_sum_sq / g_count) - (mean * mean);
            double std_dev = (var > 0.0) ? std::sqrt(var) : 0.0;

            std::cout << " DEBUG STATS : " << name << std::endl;
            std::cout << "   -> Global MIN : " << g_min << std::endl;
            std::cout << "   -> Global MAX : " << g_max << std::endl;
            std::cout << "   -> MEAN       : " << mean << std::endl;
            std::cout << "   -> STD DEV    : " << std_dev << std::endl;
            std::cout << "   -> (Nodes)    : " << g_count << std::endl;
            std::cout << "------------------------------------------------" << std::endl;
        }
    };

    // Temperature
    auto thermo_fed = heat_transfer.getFiniteElementDiscretizationPointer();
    #ifdef MFEM_USE_MPI
    if (parallel) {
        auto& thermo_fes = thermo_fed->getFiniteElementSpace<true>();
        mfem::ParGridFunction T_pgf(&thermo_fes);
        // Sync true DOFs with ghost nodes
        T_pgf.SetFromTrueDofs(heat_transfer.getUnknowns(mfem_mgis::bts));
        print_nodal_stats(T_pgf, "Temperature", false);
    } else {
        auto& thermo_fes = thermo_fed->getFiniteElementSpace<false>();
        mfem::GridFunction T_gf(&thermo_fes, heat_transfer.getUnknowns(mfem_mgis::bts).GetData());
        print_nodal_stats(T_gf, "Temperature", false);
    }
    #else
    auto& thermo_fes = thermo_fed->getFiniteElementSpace<false>();
    mfem::GridFunction T_gf(&thermo_fes, heat_transfer.getUnknowns(mfem_mgis::bts).GetData());
    print_nodal_stats(T_gf, "Temperature", false);
    #endif

    // Displacement
    auto mech_fed = mechanics.getFiniteElementDiscretizationPointer();
    #ifdef MFEM_USE_MPI
    if (parallel) {
        auto& mech_fes = mech_fed->getFiniteElementSpace<true>();
        mfem::ParGridFunction U_pgf(&mech_fes);
        U_pgf.SetFromTrueDofs(mechanics.getUnknowns(mfem_mgis::bts));
        print_nodal_stats(U_pgf, "Displacement Magnitude (||U||)", true);
    } else {
        auto& mech_fes = mech_fed->getFiniteElementSpace<false>();
        mfem::GridFunction U_gf(&mech_fes, mechanics.getUnknowns(mfem_mgis::bts).GetData());
        print_nodal_stats(U_gf, "Displacement Magnitude (||U||)", true);
    }
    #else
    auto& mech_fes = mech_fed->getFiniteElementSpace<false>();
    mfem::GridFunction U_gf(&mech_fes, mechanics.getUnknowns(mfem_mgis::bts).GetData());
    print_nodal_stats(U_gf, "Displacement Magnitude (||U||)", true);
    #endif

    // Swelling
    if (setup.swelling_model) {
        auto& m_sw = setup.swelling_model->getMaterial();
        
        const auto& vals = m_sw.s1.internal_state_variables; 
        
        double local_min = std::numeric_limits<double>::max();
        double local_max = -std::numeric_limits<double>::max();
        double local_sum = 0.0;
        double local_sum_sq = 0.0;
        
        const long long num_ips = vals.size();
        const long long local_count = num_ips;

        for (long long i = 0; i < num_ips; ++i) {
            const double val = vals[i]; 
            
            if (val < local_min) local_min = val;
            if (val > local_max) local_max = val;
            local_sum += val;
            local_sum_sq += val * val;
        }

        double g_min = local_min, g_max = local_max, g_sum = local_sum, g_sum_sq = local_sum_sq;
        long long g_count = local_count;

        #ifdef MFEM_USE_MPI
        if (parallel) {
            MPI_Comm comm = MPI_COMM_WORLD;
            auto mech_fed_local = mechanics.getFiniteElementDiscretizationPointer();
            auto& m_fes_local = mech_fed_local->getFiniteElementSpace<true>();
            if (auto pfes = dynamic_cast<const mfem::ParFiniteElementSpace*>(&m_fes_local)) {
                comm = pfes->GetComm();
            }
            MPI_Allreduce(&local_min, &g_min, 1, MPI_DOUBLE, MPI_MIN, comm);
            MPI_Allreduce(&local_max, &g_max, 1, MPI_DOUBLE, MPI_MAX, comm);
            MPI_Allreduce(&local_sum, &g_sum, 1, MPI_DOUBLE, MPI_SUM, comm);
            MPI_Allreduce(&local_sum_sq, &g_sum_sq, 1, MPI_DOUBLE, MPI_SUM, comm);
            MPI_Allreduce(&local_count, &g_count, 1, MPI_LONG_LONG, MPI_SUM, comm);
        }
        #endif

        if (mfem_mgis::getMPIrank() == 0 && g_count > 0) {
            double mean = g_sum / g_count;
            double var = (g_sum_sq / g_count) - (mean * mean);
            double std_dev = (var > 0.0) ? std::sqrt(var) : 0.0;

            std::cout << " DEBUG STATS : SolidSwelling (Fuel)" << std::endl;
            std::cout << "   -> Global MIN : " << g_min << std::endl;
            std::cout << "   -> Global MAX : " << g_max << std::endl;
            std::cout << "   -> MEAN       : " << mean << std::endl;
            std::cout << "   -> STD DEV    : " << std_dev << std::endl;
            std::cout << "   -> (Int pts)  : " << g_count << std::endl;
            std::cout << "------------------------------------------------" << std::endl;
        }
    }

    // Power 
    if (!setup.fields.empty() && setup.fields[0].Pow_s1_sw) {
        const auto& pow_vals = *setup.fields[0].Pow_s1_sw;
        
        double local_min = std::numeric_limits<double>::max();
        double local_max = -std::numeric_limits<double>::max();
        double local_sum = 0.0;
        double local_sum_sq = 0.0;
        
        const long long num_ips = pow_vals.size(); 
        const long long local_count = num_ips;

        for (long long i = 0; i < num_ips; ++i) {
            const double val = pow_vals[i]; 
            if (val < local_min) local_min = val;
            if (val > local_max) local_max = val;
            local_sum += val;
            local_sum_sq += val * val;
        }

        double g_min = local_min, g_max = local_max, g_sum = local_sum, g_sum_sq = local_sum_sq;
        long long g_count = local_count;

        #ifdef MFEM_USE_MPI
        if (parallel) {
            MPI_Comm comm = MPI_COMM_WORLD;
            auto mech_fed_local = mechanics.getFiniteElementDiscretizationPointer();
            auto& m_fes_local = mech_fed_local->getFiniteElementSpace<true>();
            if (auto pfes = dynamic_cast<const mfem::ParFiniteElementSpace*>(&m_fes_local)) {
                comm = pfes->GetComm();
            }
            MPI_Allreduce(&local_min, &g_min, 1, MPI_DOUBLE, MPI_MIN, comm);
            MPI_Allreduce(&local_max, &g_max, 1, MPI_DOUBLE, MPI_MAX, comm);
            MPI_Allreduce(&local_sum, &g_sum, 1, MPI_DOUBLE, MPI_SUM, comm);
            MPI_Allreduce(&local_sum_sq, &g_sum_sq, 1, MPI_DOUBLE, MPI_SUM, comm);
            MPI_Allreduce(&local_count, &g_count, 1, MPI_LONG_LONG, MPI_SUM, comm);
        }
        #endif

        if (mfem_mgis::getMPIrank() == 0 && g_count > 0) {
            double mean = g_sum / g_count;
            double var = (g_sum_sq / g_count) - (mean * mean);
            double std_dev = (var > 0.0) ? std::sqrt(var) : 0.0;

            std::cout << " DEBUG STATS : PowerDensity" << std::endl;
            std::cout << "   -> Global MIN : " << g_min << std::endl;
            std::cout << "   -> Global MAX : " << g_max << std::endl;
            std::cout << "   -> MEAN       : " << mean << std::endl;
            std::cout << "   -> STD DEV    : " << std_dev << std::endl;
            std::cout << "   -> (Int pts)  : " << g_count << std::endl;
            std::cout << "------------------------------------------------" << std::endl;
        }
    } else {
        if (mfem_mgis::getMPIrank() == 0) std::cout << "[WARNING] PowerDensity not found." << std::endl;
    }
}