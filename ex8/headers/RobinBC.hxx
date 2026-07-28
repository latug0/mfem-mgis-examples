/*!
 * \file   include/MFEMMGIS/RobinBC.hxx
 * \brief  Definition of the Robin boundary condition and its nonlinear form integrator.
 * \author Julien Rigal
 * \date   22/07/2026
 */

#ifndef MFEM_MGIS_ROBIN_BC_HXX
#define MFEM_MGIS_ROBIN_BC_HXX

#include <memory>
#include "MFEMMGIS/AbstractBoundaryCondition.hxx"
#include "MFEMMGIS/FiniteElementDiscretization.hxx"
#include "mfem.hpp"

namespace mfem_mgis {

  /*!
   * \brief Integrator for Robin boundary condition (Jacobian + residual)
   */
  struct RobinNonlinearFormIntegrator : public mfem::NonlinearFormIntegrator {
    //! \brief base heat transfer coefficient
    double h_base;
    //! \brief ambient temperature
    double T_inf;
    //! \brief optional displacement field to account for geometry changes
    mfem::GridFunction* u_disp;

    /*!
     * \brief constructor
     * \param[in] h_: base heat transfer coefficient
     * \param[in] T_inf_: ambient temperature
     * \param[in] u_disp_: optional displacement field
     */
    RobinNonlinearFormIntegrator(double h_, double T_inf_, mfem::GridFunction* u_disp_);
    
    /*!
     * \brief Assemble the element vector (residual)
     * \param[in] e: finite element
     * \param[in] tr: element transformation
     * \param[in] T_el: element temperatures
     * \param[out] R: assembled residual vector
     */
    void AssembleElementVector(const mfem::FiniteElement &e,
                               mfem::ElementTransformation &tr,
                               const mfem::Vector &T_el,
                               mfem::Vector &R) override;

    /*!
     * \brief Assemble the element gradient (Jacobian)
     * \param[in] e: finite element
     * \param[in] tr: element transformation
     * \param[in] elfun: element function values
     * \param[out] K: assembled local Jacobian matrix
     */
    void AssembleElementGrad(const mfem::FiniteElement &e,
                             mfem::ElementTransformation &tr,
                             const mfem::Vector &elfun,
                             mfem::DenseMatrix &K) override;
  }; // end of struct RobinNonlinearFormIntegrator

  /*!
   * \brief class used to simplify the definition of Robin boundary conditions.
   */
  class RobinBC : public AbstractBoundaryCondition {
   public:
    /*!
     * \brief constructor
     * \param[in] fed: finite element discretization
     * \param[in] tag: boundary tag
     * \param[in] h: base heat transfer coefficient
     * \param[in] T_inf: ambient temperature
     * \param[in] u_disp: optional displacement field
     */
    RobinBC(std::shared_ptr<mfem_mgis::FiniteElementDiscretization> fed, int tag, double h, double T_inf, mfem::GridFunction* u_disp);
    
    //! \brief destructor
    ~RobinBC() override;

    void setup(const real, const real) override {}

    bool addNonlinearFormIntegrator(Context &, NonlinearForm<false> &, const mfem::Vector &) noexcept override;
    bool addNonlinearFormIntegrator(Context &, NonlinearForm<true> &, const mfem::Vector &) noexcept override;

    bool addLinearFormIntegrators(Context &, BilinearForm<false> &, LinearForm<false> &,
                                  const mfem::Vector &, const real, const real) noexcept override { return true; }
    bool addLinearFormIntegrators(Context &, BilinearForm<true> &, LinearForm<true> &,
                                  const mfem::Vector &, const real, const real) noexcept override { return true; }

   private:
    mfem::Array<int> bdr_marker;
    RobinNonlinearFormIntegrator* nfi;
  }; // end of class RobinBC

} // end of namespace mfem_mgis

#endif /* MFEM_MGIS_ROBIN_BC_HXX */