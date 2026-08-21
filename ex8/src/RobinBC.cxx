/*!
 * \file   RobinBC.cxx
 * \brief  Implementation of the Robin boundary condition and its nonlinear form integrator.
 * \author Julien Rigal
 * \date   22/07/2026
 */

#include "../headers/RobinBC.hxx"
#include <algorithm> 
#include <cmath>

namespace mfem_mgis {

  RobinNonlinearFormIntegrator::RobinNonlinearFormIntegrator(double h_, double T_inf_, mfem::GridFunction* u_disp_)
      : h_base(h_), T_inf(T_inf_), u_disp(u_disp_) {
  } // end of RobinNonlinearFormIntegrator

  void RobinNonlinearFormIntegrator::AssembleElementVector(const mfem::FiniteElement &e,
                                                          mfem::ElementTransformation &tr,
                                                          const mfem::Vector &T_el,
                                                          mfem::Vector &R) {
      const int nnodes = e.GetDof();
      mfem::Vector shape(nnodes);
      R.SetSize(nnodes);
      R = 0.0;
      const auto *ir = this->IntRule;
      if (!ir) ir = &mfem::IntRules.Get(e.GetGeomType(), e.GetOrder());

      mfem::DenseMatrix grad_u;
      mfem::FaceElementTransformations *face_tr = dynamic_cast<mfem::FaceElementTransformations *>(&tr);

      for (int i = 0; i < ir->GetNPoints(); ++i) {
          const auto &ip = ir->IntPoint(i);
          tr.SetIntPoint(&ip);
          e.CalcShape(ip, shape);
          const double T_pt = shape * T_el;
          const double w = ip.weight * tr.Weight();
          
          double h_eff = h_base;
          if (u_disp != nullptr && face_tr != nullptr) {
              mfem::IntegrationPoint eip;
              face_tr->Loc1.Transform(ip, eip);
              face_tr->Elem1->SetIntPoint(&eip);

              u_disp->GetVectorGradient(*(face_tr->Elem1), grad_u);
              double dilatation_volumique = grad_u(0,0) + grad_u(1,1) + grad_u(2,2);
              h_eff = h_base * (1.0 + dilatation_volumique * (2.0 / 3.0));
          }

          for (int ni = 0; ni < nnodes; ++ni)
              R[ni] -= w * h_eff * (T_pt - T_inf) * shape[ni]; 
      }
  } // end of AssembleElementVector

  void RobinNonlinearFormIntegrator::AssembleElementGrad(const mfem::FiniteElement &e,
                                                        mfem::ElementTransformation &tr,
                                                        const mfem::Vector &,
                                                        mfem::DenseMatrix &K) {
      const int nnodes = e.GetDof();
      mfem::Vector shape(nnodes);
      K.SetSize(nnodes, nnodes);
      K = 0.0;
      const auto *ir = this->IntRule;
      if (!ir) ir = &mfem::IntRules.Get(e.GetGeomType(), e.GetOrder());

      mfem::DenseMatrix grad_u;
      mfem::FaceElementTransformations *face_tr = dynamic_cast<mfem::FaceElementTransformations *>(&tr);
      
      for (int i = 0; i < ir->GetNPoints(); ++i) {
          const auto &ip = ir->IntPoint(i);
          tr.SetIntPoint(&ip);
          e.CalcShape(ip, shape);
          const double w = ip.weight * tr.Weight();
          
          double h_eff = h_base;
          if (u_disp != nullptr && face_tr != nullptr) {
              mfem::IntegrationPoint eip;
              face_tr->Loc1.Transform(ip, eip);
              face_tr->Elem1->SetIntPoint(&eip);
              
              u_disp->GetVectorGradient(*(face_tr->Elem1), grad_u);
              double dilatation_volumique = grad_u(0,0) + grad_u(1,1) + grad_u(2,2);
              h_eff = h_base * (1.0 + dilatation_volumique * (2.0 / 3.0));
          }

          for (int ni = 0; ni < nnodes; ++ni)
              for (int nj = 0; nj < nnodes; ++nj)
                  K(ni, nj) -= w * h_eff * shape[ni] * shape[nj]; 
      }
  } // end of AssembleElementGrad

  RobinBC::RobinBC(std::shared_ptr<mfem_mgis::FiniteElementDiscretization> fed, int tag, double h, double T_inf, mfem::GridFunction* u_disp)
      : nfi(new RobinNonlinearFormIntegrator(h, T_inf, u_disp)) {
      
      int local_max = 0;
      
      if (fed->describesAParallelComputation()) {
  #ifdef MFEM_USE_MPI
          auto& fes = fed->getFiniteElementSpace<true>();
          auto* mesh = fes.GetMesh();
          local_max = (mesh->bdr_attributes.Size() > 0) ? mesh->bdr_attributes.Max() : 0;
  #else
          mfem_mgis::raise("RobinBC: unsupported parallel computations");
  #endif
      } else {
            auto& fes = fed->getFiniteElementSpace<false>();
            auto* mesh = fes.GetMesh();
            local_max = (mesh->bdr_attributes.Size() > 0) ? mesh->bdr_attributes.Max() : 0;
      }

      bdr_marker.SetSize(std::max(local_max, tag));
      bdr_marker = 0;
      bdr_marker[tag - 1] = 1;
  } // end of RobinBC

  RobinBC::~RobinBC() {
  } // end of ~RobinBC

  bool RobinBC::addNonlinearFormIntegrator(Context &, NonlinearForm<false> &f, const mfem::Vector &) noexcept {
      f.AddBoundaryIntegrator(nfi, bdr_marker);
      return true;
  } // end of addNonlinearFormIntegrator

  bool RobinBC::addNonlinearFormIntegrator(Context &, NonlinearForm<true> &f, const mfem::Vector &) noexcept {
      f.AddBoundaryIntegrator(nfi, bdr_marker);
      return true;
  } // end of addNonlinearFormIntegrator

} // end of namespace mfem_mgis