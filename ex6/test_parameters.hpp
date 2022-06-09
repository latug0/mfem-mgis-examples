#pragma once
#include "timer.hpp"
#include <memory>
#include <cstdlib>
#include <iostream>
#include "mfem/general/optparser.hpp"
#include "mfem/linalg/solvers.hpp"
#include "mfem/fem/datacollection.hpp"
#include "MFEMMGIS/MFEMForward.hxx"
#include "MFEMMGIS/Material.hxx"
#include "MFEMMGIS/Profiler.hxx"
#include "MFEMMGIS/AnalyticalTests.hxx"
#include "MFEMMGIS/NonLinearEvolutionProblemImplementation.hxx"
#include "MFEMMGIS/PeriodicNonLinearEvolutionProblem.hxx"

#ifdef MFEM_USE_PETSC
#include "mfem/linalg/petsc.hpp"
#endif /* MFEM_USE_PETSC */

#ifdef MFEM_USE_PETSC
#include "mfem/linalg/mumps.hpp"
#endif /* MFEM_USE_MUMPS */

// We need this class for test case sources
constexpr double xmax = 1.;
struct TestParameters {
	const char* mesh_file = "cas-cible-1.mesh";
	const char* behaviour = "Elasticity";
	const char* library = "src/libBehaviour.so";
	const char* reference_file = "Elasticity.ref";
	int order = 1;
	int tcase = 1;
	double xmax = 1.;
	double ymax = 1.;
	double zmax = 1.;
	int linearsolver = 1;
	int preconditioner = 1;
	bool parallel = true;
};

TestParameters parseCommandLineOptions(int& argc, char* argv[]);
