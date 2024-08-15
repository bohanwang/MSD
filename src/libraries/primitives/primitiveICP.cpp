#include "primitiveICP.h"
#include "templatePrimitive.h"

#include "pgoLogging.h"
#include "cgalInterface.h"
#include "libiglInterface.h"
#include "EigenSupport.h"
#include "EigenMKLPardisoSupport.h"
#include "geometryQuery.h"
#include "createTriMesh.h"
#include "triMeshNeighbor.h"
#include "predicates.h"

#include <tbb/parallel_for.h>
#include <tbb/enumerable_thread_specific.h>
#include <tbb/spin_mutex.h>

#include <fmt/format.h>

#include <chrono>
#include <unordered_set>
#include <filesystem>
#include <numeric>

namespace MedialAxisRepresentation
{
namespace ES = pgo::EigenSupport;

struct PrimitiveICPSolverData
{
  PrimitiveICPSolverData(int nv)
  {
    numVtx = nv;
  }

  void init()
  {
    primitiveDCurr.resize(numVtx);
    primitiveDCurr.setZero();

    primitiveDDCurr.resize(numVtx);
    primitiveDDCurr.setZero();

    primitiveDDProj.resize(numVtx);
    primitiveDDProj.setZero();

    primitiveDDelta.resize(numVtx);
    primitiveDDelta.setZero();

    primitiveDDMax.resize(numVtx);
    bSmoothness.resize(numVtx);
    bExpansion.resize(numVtx);
    bContact.resize(numVtx);
    bProj.resize(numVtx);
    rhs.resize(numVtx);

    AExpansion.resize(numVtx, numVtx);
    AExpansion.setIdentity();

    AContact.resize(numVtx, numVtx);
    ES::VXd IVec(numVtx);
    IVec.setZero();

    pgo::libiglInterface::diag(IVec, AContact);

    AProj.resize(numVtx, numVtx);
    AProj.setIdentity();

    // ES::Mm(graphL, graphL, graphLTL, 1);
    graphLTL = graphL * (-1);

    sys = graphLTL + AExpansion + AContact;
    rhs = bExpansion + bContact;

    // hessianMatrixMapping
    std::vector<int> dofs(numVtx);
    std::iota(dofs.begin(), dofs.end(), 0);

    ES::SpMatI mappingL;
    ES::small2Big(graphLTL, sys, dofs, mappingL);

    ES::SpMatI mappingExpansion;
    ES::small2Big(AExpansion, sys, dofs, mappingExpansion);

    ES::SpMatI mappingContact;
    ES::small2Big(AContact, sys, dofs, mappingContact);

    ES::SpMatI mappingProj;
    ES::small2Big(AProj, sys, dofs, mappingProj);
    hessianMatrixMappings = { mappingL, mappingExpansion, mappingContact, mappingProj };

    // solver
    // solver = std::make_shared<ES::EigenPardisoSupport>(sys, ES::EigenPardisoSupport::MatrixType::REAL_SYM_INDEFINITE,
    //   ES::EigenPardisoSupport::ReorderingType::NESTED_DISSECTION, 0, 0, 0, 0, 0, 0);
    // solver->analyze(sys);

    // debug eigen solver
    eigenSolver = std::make_shared<ES::SymSolver>();
    eigenSolver->analyzePattern(sys);
  }

  // update according to primiviteDDelta
  void update(const std::shared_ptr<TemplatePrimitive> &primitiveTemp)
  {
    // tbb::parallel_for(0, numVtx, [&](int i) {
    for (int i = 0; i < numVtx; i++) {
      ES::V3d dir_i = primitiveTemp->rayDir.row(i);
      ES::V3d init_i = primitiveCurrPos.row(i);
      ES::V3d v_i = primitiveDDelta(i) * dir_i + init_i;

      primitiveCurrPos.row(i) = v_i;
      primitiveMeshCurr.pos(i) = v_i;
      primitiveDDCurr(i) += primitiveDDelta(i);
      primitiveDCurr(i) += primitiveDDelta(i);

      // PGO_ALOG(std::abs(primitiveDCurr(i) - (primitiveCurrPos.row(i) - primitiveRestPos.row(i)).norm() - primitiveTemp->radius) < 1e-6);
    }
    // });
    primitiveMeshPseudoNormal.updateVertexPositions(primitiveMeshCurr);
  }

  pgo::Mesh::TriMeshGeo primitiveMeshCurr, primitiveMeshRest;
  pgo::Mesh::TriMeshPseudoNormal primitiveMeshPseudoNormal;
  ES::MXd primitiveRestPos;
  ES::MXi primitiveRestTri;
  ES::MXd primitiveCurrPos;

  ES::VXd primitiveDCurr;
  ES::VXd primitiveDDProj;
  ES::VXd primitiveDDCurr;
  ES::VXd primitiveDDelta;
  ES::VXd primitiveDDMax;

  // Smoothness
  ES::SpMatD graphL;
  ES::SpMatD graphLTL;
  ES::VXd bSmoothness;

  // Expansion
  ES::SpMatD AExpansion;  // is identy matrix and is fixed, same as ATA
  ES::VXd bExpansion;     // same as ATb

  // Contact
  ES::SpMatD AContact;  // is a diagonal matrix
  ES::VXd bContact;

  // Projective Constraints
  ES::SpMatD AProj;  // is a diagonal matrix
  ES::VXd bProj;

  // All Energies
  ES::SpMatD sys;
  ES::VXd rhs;

  // Store three index mapping for hessian assembly
  std::vector<ES::SpMatI, Eigen::aligned_allocator<ES::SpMatI>> hessianMatrixMappings;

  // Solver
  std::shared_ptr<ES::EigenMKLPardisoSupport> solver;
  std::shared_ptr<ES::SymSolver> eigenSolver;

  int numVtx;
};

PrimitiveICPSolver::PrimitiveICPSolver(const ES::MXd &centers, const ES::VXd &centerRadii, const PrimitiveICPSolverParameters &param_, int primitiveType_):
  params(param_), primitiveType(primitiveType_)
{
  primitiveTemp = std::make_shared<TemplatePrimitive>(param_.sphereRadius, param_.numCylinderSlices, primitiveType_);
  primitiveTemp->init(centers, centerRadii);  // transform to match the input centers

  int nv = primitiveTemp->primitiveTemplateMesh.numVertices();
  pd = new PrimitiveICPSolverData(nv);

  // build primitive normals
  pd->primitiveMeshRest = primitiveTemp->primitiveTemplateMesh;
  pd->primitiveMeshCurr = primitiveTemp->primitiveTemplateMesh;
  pd->primitiveMeshPseudoNormal.buildPseudoNormals(pd->primitiveMeshCurr);

  // Smoothness
  pgo::Mesh::triMeshGeoToMatrices(primitiveTemp->primitiveTemplateMesh, pd->primitiveRestPos, pd->primitiveRestTri);
  pd->primitiveCurrPos = pd->primitiveRestPos;

  pgo::libiglInterface::graphLaplacianMatrix(pd->primitiveRestTri, pd->graphL);

  // init
  pd->init();

  // TODO base hessian, probably used by barycentric coordinate
}

PrimitiveICPSolver::~PrimitiveICPSolver()
{
  delete pd;
}

int PrimitiveICPSolver::sim(const ES::MXd &centers, const pgo::Mesh::TriMeshGeo &targetMesh, const pgo::Mesh::TriMeshBVTree &targetMeshBVTree, const pgo::Mesh::TriMeshPseudoNormal &targetMeshNormals,
  const std::string &prefix, pgo::Mesh::TriMeshGeo &meshOut)
{
  PGO_ALOG(centers.rows() == primitiveType);
  PGO_ALOG(centers.cols() == 3);

  // Update primitiveMeshRest, primitiveMeshCurr, primitiveMeshPseudoNormal, primitiveCurrPos
  pd->primitiveMeshRest = primitiveTemp->primitiveTemplateMesh;
  pd->primitiveMeshCurr = primitiveTemp->primitiveTemplateMesh;
  pd->primitiveMeshPseudoNormal.updateVertexPositions(pd->primitiveMeshCurr);

  pgo::Mesh::BoundingBox bb(targetMesh.positions());
  double goodDelta = bb.diameter() * 0.1;

  int numVtx = pd->primitiveMeshRest.numVertices();

  // Update primitiveDDMax, primitiveDCurr
  for (int i = 0; i < numVtx; i++) {
    ES::V3d start = pd->primitiveMeshRest.pos(i);
    ES::V3d dir = primitiveTemp->rayDir.row(i);
    start -= dir / dir.norm() * primitiveTemp->radius;
    ES::V3d end = start + dir / dir.norm() * 10;

    double segWeight[2], triWeights[3];
    int triID = targetMeshBVTree.lineSegmentFirstIntersectionPointSemiExact(targetMesh, start, end, segWeight, triWeights);

    bool foundDMax = false;
    if (triID >= 0) {
      ES::V3d n = targetMeshNormals.triNormal(triID);
      if (n.dot(dir) > 0) {
        pd->primitiveDDMax(i) = segWeight[1] * (start - end).norm() - primitiveTemp->radius;
        foundDMax = true;
      }
    }

    // if (!foundDMax) {
    //   end = start - ES::toVec3d(dir / dir.norm() * 10);
    //   triID = targetMeshBVTree.lineSegmentFirstIntersectionPointSemiExact(targetMesh, start, end, segWeight, triWeights);
    //   if (triID >= 0) {
    //     ES::V3d n = targetMeshNormals.triNormal(triID);
    //     if (dot(n, ES::V3d(dir.data())) > 0) {
    //       pd->primitiveDDMax(i) = -segWeight[1] * len(start - end) - primitiveTemp->radius;
    //       foundDMax = true;
    //     }
    //   }
    // }

    if (foundDMax != true) {
      return 1;
    }

    pd->primitiveDCurr(i) = primitiveTemp->radius;
  }

  // initialize primitiveDDCurr with min(DMax)
  // double minDMax = pd->primitiveDDMax.minCoeff();
  pd->primitiveDDelta.setConstant(0.0);
  pd->update(primitiveTemp);

  // params.expansionRadius = minDMax;
  // fmt::print("expansion radius: {}\n", params.expansionRadius);

  pd->primitiveMeshCurr.save(fmt::format("{}/init.obj", prefix));
  int niter = params.maxNumIter;
  if (centers.rows() == 1)
    niter = 3;

  for (int iter = 0; iter < niter; iter++) {
    // double smoothnessCoeff = 1e-1, expansionCoeff = 1.0, contactCoeff = 1.0; // sphere
    // double smoothnessCoeff = 1.0, expansionCoeff = 1.0, contactCoeff = 1.0; // cylinder
    double smoothnessCoeff = params.smoothnessCoeff * std::pow(params.smoothnessDecay, iter),
           expansionCoeff = params.expansionCoeff,
           contactCoeff = params.contactCoeff;

    if (params.verbose) {
      fmt::print("Iter: {}\n", iter);
    }
    // fmt::print("Alpha: {}\n", w_cur);

    if (params.solType == SolverType::LINEARSOLVER) {  // WRONG, deprecated
      // For sys, only need to update AContact
      // For rhs, need to update bContact and bExpansion
      for (int i = 0; i < numVtx; i++) {
        ES::V3d q_i = pd->primitiveMeshCurr.pos(i);
        auto ret = targetMeshBVTree.closestTriangleQuery(targetMesh, q_i);
        ES::V3d dir = ret.closestPosition - q_i;
        ES::V3d target_n = targetMeshNormals.getPseudoNormal(targetMesh.triangles().data(), ret.triID, ret.feature);
        double dist = dir.norm();
        dir /= dist;
        double dotValue = dir.dot(target_n);

        pd->bExpansion(i) = std::min(pd->primitiveDDMax(i) - pd->primitiveDDCurr(i), 0.2 * params.expansionRadius) * (-1);

        double AVal, bVal;
        if (dotValue <= 0 || dist < 1e-3) {
          ES::M3d nnT;
          ES::V3d n = target_n;
          ES::tensorProduct(nnT, n, n);
          ES::V3d dir_i = primitiveTemp->rayDir.row(i);
          ES::V3d closPos = ret.closestPosition;
          ES::V3d init_i = pd->primitiveMeshCurr.pos(i);

          AVal = dir_i.transpose() * nnT * dir_i;
          bVal = (init_i - closPos).transpose() * nnT * dir_i;
        }
        else {
          AVal = 0;
          bVal = 0;
        }
        double prevACoeff = pd->AContact.coeff(i, i);
        pd->AContact.coeffRef(i, i) += AVal - prevACoeff;
        pd->bContact(i) = bVal;
      }

      // Solve
      std::memset(pd->sys.valuePtr(), 0, sizeof(double) * pd->sys.nonZeros());
      ES::addSmallToBig(smoothnessCoeff, pd->graphLTL, pd->sys, 1.0, pd->hessianMatrixMappings[0]);
      ES::addSmallToBig(expansionCoeff, pd->AExpansion, pd->sys, 1.0, pd->hessianMatrixMappings[1]);
      ES::addSmallToBig(contactCoeff, pd->AContact, pd->sys, 1.0, pd->hessianMatrixMappings[2]);

      pd->rhs = contactCoeff * pd->bContact;
      pd->rhs += expansionCoeff * pd->bExpansion;
      pd->rhs *= -1;

      // ES::VXd x(pd->primitiveDDelta.size());
      // x.setZero();
      // ES::SpMatD sysDebug = pd->sys;
      // ES::VXd rhsDebug = pd->rhs;
      // auto solver = std::make_shared<ES::EigenPardisoSupport>(sysDebug, ES::EigenPardisoSupport::MatrixType::REAL_SYM_INDEFINITE,
      //   ES::EigenPardisoSupport::ReorderingType::NESTED_DISSECTION, 0, 1, 1, 0, 0, 0);
      // solver->analyze(sysDebug);
      // solver->factorize(sysDebug);
      // solver->solve(x.data(), rhsDebug.data(), 1);
      // ES::VXd solveRes(x.size());
      // ES::Mv(sysDebug, x, solveRes);
      // solveRes -= rhsDebug;
      // double solveResNorm = solveRes.norm();
      // fmt::print("Solve Res Norm: {}\n", solveResNorm);

      pd->eigenSolver->factorize(pd->sys);
      pd->primitiveDDelta.setZero();
      pd->primitiveDDelta = pd->eigenSolver->solve(pd->rhs);

      // update
      pd->update(primitiveTemp);
    }
    else if (params.solType == SolverType::PDSOLVER) {
      // set expansion
      for (int i = 0; i < numVtx; i++) {
        double dmax = std::min(primitiveTemp->rayInitialLength(i) * 0.1, goodDelta);

        double rmin = 1e100;
        for (int vi : primitiveTemp->primitiveTemplateMeshVertexNeighboringVertices[i]) {
          rmin = std::min(pd->primitiveDCurr(vi), rmin);
        }
        double diff = pd->primitiveDCurr(i) - rmin;

        dmax = std::min(dmax, dmax * 2 - diff);
        dmax = std::max(dmax, 0.0);

        // std::cout << diff << ',' << dmax << std::endl;

        pd->bExpansion(i) = std::min(pd->primitiveDDMax(i) - pd->primitiveDDCurr(i), dmax) * (-1);
      }

      pd->bSmoothness.resize(numVtx);
      ES::mv(pd->graphLTL, pd->primitiveDCurr, pd->bSmoothness);

      pd->primitiveDDelta.setZero();

      ES::VXd dlast = pd->primitiveDDelta;
      // projective dynamics

      for (int PDIter = 0; PDIter < params.maxPDNumIter; PDIter++) {
        if (params.verbose) {
          fmt::print("PDIter: {}\n", PDIter);
        }
        for (int i = 0; i < numVtx; i++) {
          // project to the constraint set (contact)
          // directly project to the constraint set
          // update AProj, bProj
          // fmt::print("pd->AProj: {}, numVtx: {}\n", pd->AProj.nonZeros(), numVtx);
          if (pd->primitiveDDCurr(i) + pd->primitiveDDelta(i) > pd->primitiveDDMax(i)) {
            pd->primitiveDDProj(i) = pd->primitiveDDMax(i);
            pd->AProj.valuePtr()[i] = 1;
            pd->bProj(i) = pd->primitiveDDCurr(i) - pd->primitiveDDProj(i);
          }
          else {
            pd->primitiveDDProj(i) = pd->primitiveDDCurr(i) + pd->primitiveDDelta(i);
            pd->AProj.valuePtr()[i] = 0;
            pd->bProj(i) = 0;
          }
        }
        // Solve
        std::memset(pd->sys.valuePtr(), 0, sizeof(double) * pd->sys.nonZeros());
        ES::addSmallToBig(smoothnessCoeff, pd->graphLTL, pd->sys, 1.0, pd->hessianMatrixMappings[0]);
        ES::addSmallToBig(expansionCoeff, pd->AExpansion, pd->sys, 1.0, pd->hessianMatrixMappings[1]);

        pd->rhs = expansionCoeff * pd->bExpansion;
        ES::addSmallToBig(contactCoeff, pd->AProj, pd->sys, 1.0, pd->hessianMatrixMappings[2]);
        pd->rhs += contactCoeff * pd->bProj;
        pd->rhs += smoothnessCoeff * pd->bSmoothness;

        pd->rhs *= -1;
        pd->eigenSolver->factorize(pd->sys);
        pd->primitiveDDelta.setZero();
        pd->primitiveDDelta = pd->eigenSolver->solve(pd->rhs);

        double maxDiff = (pd->primitiveDDelta - dlast).cwiseAbs().maxCoeff();
        if (params.verbose) {
          fmt::print("  |dx|_inf={}\n", maxDiff);
          std::cout << std::flush;
        }

        if (maxDiff < 1e-5)
          break;

        dlast.noalias() = pd->primitiveDDelta;
      }

      // project to constraint set at the final iteration
      if (iter == params.maxNumIter - 1) {
        for (int i = 0; i < numVtx; i++) {
          if (pd->primitiveDDCurr(i) + pd->primitiveDDelta(i) > pd->primitiveDDMax(i)) {
            pd->primitiveDDelta(i) = pd->primitiveDDMax(i) - pd->primitiveDDCurr(i);
          }
        }
      }
      // update
      pd->update(primitiveTemp);
    }

    if (params.verbose) {
      // save
      pd->primitiveMeshCurr.save(fmt::format("{}/iter{}.obj", prefix, iter));
    }
  }
  // return
  meshOut = pd->primitiveMeshCurr;
  primitiveTemp->update(centers, pd->primitiveDCurr);
  return 0;
}
}  // namespace MedialAxisRepresentation