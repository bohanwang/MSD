#include "sphereICP.h"

#include "pgoLogging.h"
#include "CGALUtilities.h"
#include "basicUtilities.h"
#include "libiglInterface.h"

#include "EigenSupport.h"
#include "geometryQuery.h"
#include "createTriMesh.h"
#include "finiteDifference.h"
#include "automaticDifferentiation_autodiff.h"
#include "triMeshNeighbor.h"
#include "minimizeEnergy.h"
#include "predicates.h"

#include <tbb/parallel_for.h>
#include <tbb/enumerable_thread_specific.h>
#include <tbb/spin_mutex.h>

#include <fmt/format.h>

#include <unordered_set>
#include <filesystem>

using namespace MedialAxisRepresentation;

namespace MedialAxisRepresentation
{
namespace ES = pgo::EigenSupport;

struct ICPData
{
  ES::VXd coeffs;
  ES::VXd targetPositions;
  ES::VXd targetNormals;
  std::vector<int> vertexIndices;

  ICPData(int n):
    coeffs(n), targetPositions(n * 3), targetNormals(n * 3), vertexIndices(n)
  {
    std::iota(vertexIndices.begin(), vertexIndices.end(), 0);
  }
};

struct SphereICPSolverData
{
  SphereICPSolverData(int nv, int ntri, int numPositionDOFs, bool movableCenter, bool hasConstraints):
    icpData(nv), expansionData(nv), vtxAreas(nv), distAll(nv), icpSampleData(nv), expansionSampleData(nv)
  {
    int numExtraDOFs = 0;
    if (movableCenter) {
      numExtraDOFs = 12;
    }

    qlow.resize(nv * numPositionDOFs + numExtraDOFs);
    qhi.resize(nv * numPositionDOFs + numExtraDOFs);
    q.resize(nv * numPositionDOFs + numExtraDOFs);

    curVel.resize(nv * 3);
    curPosition.resize(nv * 3);
    curPosition1.resize(nv * 3);
    restPosition.resize(nv * 3);
    temp.resize(nv * 3);

    if (hasConstraints) {
      ALOG(numPositionDOFs == 3);

      g.resize(ntri);
      lambda.resize(ntri);
      clow.resize(ntri);
      chi.resize(ntri);
    }

    edgeVertexIDs << 1, 2,
      2, 0,
      0, 1;
  }

  ~SphereICPSolverData()
  {
    ES::DestroySymbolicMmData(mulL3);
    ES::DestroySymbolicMmData(mulL12);
  }

  TriMeshGeo dumpConstraints(bool useSamples = false) const
  {
    if (useSamples) {
      // dump constraints
      TriMeshGeo abc;
      for (int si = 0; si < (int)sampleRestCoeffs.size(); si++) {
        Vec3d p0(0.0);
        for (int vi = 0; vi < 3; vi++) {
          p0 += sphereMeshCur.pos(sampleVertexIndices[si * 3 + vi]) * sampleVertexWeights[si * 3 + vi];
        }

        if (icpSampleData.coeffs[si] > 0) {
          Vec3d p1(icpSampleData.targetPositions.data() + si * 3);
          abc.addMesh(createSingleTriangleMesh(p0, p1, p0 + Vec3d(1e-5)));
        }
        else {
          Vec3d p1(expansionSampleData.targetPositions.data() + si * 3);
          abc.addMesh(createSingleTriangleMesh(p0, p1, p0 + Vec3d(1e-5)));
        }
      }
      return abc;
    }
    else {
      // dump constraints
      TriMeshGeo abc;
      for (int i = 0; i < sphereMeshCur.numVertices(); i++) {
        Vec3d p0 = sphereMeshCur.pos(i);

        if (icpData.coeffs[i] > 0) {
          Vec3d p1(icpData.targetPositions.data() + i * 3);
          abc.addMesh(createSingleTriangleMesh(p0, p1, p0 + Vec3d(1e-5)));
        }
        else {
          Vec3d p1(expansionData.targetPositions.data() + i * 3);
          abc.addMesh(createSingleTriangleMesh(p0, p1, p0 + Vec3d(1e-5)));
        }
      }
      return abc;
    }
  }

  TriMeshGeo sphereMeshCur, sphereMeshRest;
  TriMeshPseudoNormal sphereMeshNormals;
  ES::MXd sphereRestPostions;
  ES::MXi sphereTriangles;

  std::vector<int> sphereSamples;

  ES::VXd g, lambda, clow, chi;

  ES::VXd qlow, qhi;
  ES::VXd q, qvel, q1;
  ES::VXd restPosition, curPosition, curPosition1, curVel;
  ES::VXd temp;

  // vertex constraints
  ICPData icpData, expansionData;
  ES::VXd vtxAreas;
  ES::VXd distAll;
  ES::VXd lastExpansionTargetPosition;

  // sample constraints
  std::vector<int> sampleVertexIndices;
  std::vector<double> sampleVertexWeights;
  ES::VXd sampleRestCoeffs;
  ICPData icpSampleData, expansionSampleData;
  ES::VXd lastExpansionTargetSamplePosition;

  std::unordered_set<int> invalidConstraintsSet;

  ES::SpMatD biL3, M3, MInvL3;
  ES::SpMatD biL12, M12, MInvL12;
  ES::VXd MVec;
  ES::SpMatD baseHessian;

  // inertial
  ES::SpMatD A_inertial;
  ES::VXd b_inertial;
  std::shared_ptr<VegaFEM::PredefinedPotentialEnergies::QuadraticPotentialEnergy> inertiaEnergy;

  // icp
  std::shared_ptr<VegaFEM::NonlinearOptimization::PotentialEnergy> icpEnergy;
  std::shared_ptr<VegaFEM::PredefinedPotentialEnergies::MultipleVertexSliding> icpEnergyQuadratic;

  std::shared_ptr<VegaFEM::PredefinedPotentialEnergies::BarycentricCoordinateSliding> icpEnergyQuadraticDenseSampled;
  std::shared_ptr<VegaFEM::PredefinedPotentialEnergies::MultipleVertexPullingSoftConstraintsPOrder> icpEnergyPOrder;

  // smoothness
  ES::VXd b_smooth;
  std::shared_ptr<VegaFEM::PredefinedPotentialEnergies::QuadraticPotentialEnergy> smoothnessEnergyXDOF;
  std::shared_ptr<VegaFEM::PredefinedPotentialEnergies::SurfaceSmoothnessAbsoluteMeanCurvature> smoothnessEnergyAbsHXDOF;
  std::shared_ptr<VegaFEM::PredefinedPotentialEnergies::SurfaceTriangleDeformation> smoothnessEnergyTriangleStrainXDOF;
  std::shared_ptr<VegaFEM::PredefinedPotentialEnergies::QuadraticPotentialEnergy> smoothnessEnergyAffineDOF;

  // expansion
  ES::VXd expansionForce;
  std::shared_ptr<VegaFEM::PredefinedPotentialEnergies::MultipleVertexSliding> sphereExpansionEnergyQuad;
  std::shared_ptr<VegaFEM::PredefinedPotentialEnergies::BarycentricCoordinateSliding> sphereExpansionEnergyQuadDenseSampled;
  std::shared_ptr<VegaFEM::NonlinearOptimization::PotentialEnergy> expansionEnergy;
  std::shared_ptr<VegaFEM::PredefinedPotentialEnergies::LinearPotentialEnergy> sphereExpansionEnergyLinear;

  // translation
  std::shared_ptr<VegaFEM::PredefinedPotentialEnergies::MultipleVertexConstrainedRigidMotion> sphereTranslationEnergy;

  ES::SpMatD A_centerReg;
  std::shared_ptr<VegaFEM::PredefinedPotentialEnergies::QuadraticPotentialEnergy> sphereCenterRegEnergy;

  // barrier
  ES::VXd triangleRestVols;
  ES::VXd triangleRestAreas;
  std::shared_ptr<SphereInversionBarrierEnergy> inversionBarrierEnergy;

  std::shared_ptr<VegaFEM::PredefinedPotentialEnergies::VertexAffineToPositionEnergy> affineSpaceEnergies;

  // all energy
  std::shared_ptr<VegaFEM::NonlinearOptimization::PotentialEnergies> energiesXDOFs, energiesAffineDOFs;
  std::shared_ptr<SphereTriangleFlippingConstraints> constraints;

  ICPSolver icpSolver;

  // matrix computation
  void computeNewLaplacianAndMass(const TriMeshGeo &mesh);
  void computeCotWeights(const TriMeshGeo &mesh);
  void updateLaplacianAndMass();

  ES::MX3d edgeCotangents, massQuads;
  Eigen::Matrix<int, 3, 2> edgeVertexIDs;
  std::vector<std::array<std::array<std::tuple<std::ptrdiff_t, double>, 4>, 3>> edgeWeightDestinations;

  ES::SpMatD L, L3, L12;
  ES::VXi LToL3Mapping, LToL12Mapping;
  ES::SymbolicMmData *mulL3 = nullptr, *mulL12 = nullptr;
};

void SphereICPSolverData::computeCotWeights(const TriMeshGeo &mesh)
{
  int m = mesh.numTriangles();
  tbb::parallel_for(
    0, m, [&](int i) {
      ES::V3d edgeLengths2;
      edgeLengths2(0) = len2(mesh.pos(i, 1) - mesh.pos(i, 2));
      edgeLengths2(1) = len2(mesh.pos(i, 2) - mesh.pos(i, 0));
      edgeLengths2(2) = len2(mesh.pos(i, 0) - mesh.pos(i, 1));

      ES::V3d p0 = ES::toESV3d(mesh.pos(i, 0));
      ES::V3d p1 = ES::toESV3d(mesh.pos(i, 1));
      ES::V3d p2 = ES::toESV3d(mesh.pos(i, 2));

      ES::V3d e0 = p1 - p0;
      ES::V3d e1 = p2 - p0;

      double dblA = e0.cross(e1).norm();

      // cotangents and diagonal entries for element matrices
      // correctly divided by 4 (alec 2010)
      // Alec: I'm doubtful that using l2 here is actually improving numerics.
      edgeCotangents(i, 0) = (edgeLengths2(1) + edgeLengths2(2) - edgeLengths2(0)) / dblA / 4.0;
      edgeCotangents(i, 1) = (edgeLengths2(2) + edgeLengths2(0) - edgeLengths2(1)) / dblA / 4.0;
      edgeCotangents(i, 2) = (edgeLengths2(0) + edgeLengths2(1) - edgeLengths2(2)) / dblA / 4.0;

      // now mass matrix
      ES::V3d cosines;
      // cosines(0) = (l.col(2).array().pow(2) + l.col(1).array().pow(2) - l.col(0).array().pow(2)) / (l.col(1).array() * l.col(2).array() * 2.0);
      cosines(0) = (edgeLengths2(2) + edgeLengths2(1) - edgeLengths2(0)) / (std::sqrt(edgeLengths2(1) * edgeLengths2(2)) * 2.0);
      // cosines(1) = (l.col(0).array().pow(2) + l.col(2).array().pow(2) - l.col(1).array().pow(2)) / (l.col(2).array() * l.col(0).array() * 2.0);
      cosines(1) = (edgeLengths2(0) + edgeLengths2(2) - edgeLengths2(1)) / (std::sqrt(edgeLengths2(2) * edgeLengths2(0)) * 2.0);
      // cosines(2) = (l.col(1).array().pow(2) + l.col(0).array().pow(2) - l.col(2).array().pow(2)) / (l.col(0).array() * l.col(1).array() * 2.0);
      cosines(2) = (edgeLengths2(1) + edgeLengths2(0) - edgeLengths2(2)) / (std::sqrt(edgeLengths2(0) * edgeLengths2(1)) * 2.0);

      // Matrix<Scalar, Dynamic, 3> barycentric = cosines.array() * l.array();
      // normalize_row_sums(barycentric, barycentric);
      ES::V3d weightedEdgeLength = cosines.array() * edgeLengths2.array().sqrt();
      weightedEdgeLength /= weightedEdgeLength.sum();

      // Matrix<Scalar, Dynamic, 3> partial = barycentric;
      // partial.col(0).array() *= dblA.array() * 0.5;
      // partial.col(1).array() *= dblA.array() * 0.5;
      // partial.col(2).array() *= dblA.array() * 0.5;
      ES::V3d partial = weightedEdgeLength * dblA * 0.5;

      // Matrix<Scalar, Dynamic, 3> quads(partial.rows(), partial.cols());
      ES::V3d quad;
      // quads.col(0) = (partial.col(1) + partial.col(2)) * 0.5;
      quad(0) = (partial(1) + partial(2)) * 0.5;
      // quads.col(1) = (partial.col(2) + partial.col(0)) * 0.5;
      quad(1) = (partial(2) + partial(0)) * 0.5;
      // quads.col(2) = (partial.col(0) + partial.col(1)) * 0.5;
      quad(2) = (partial(0) + partial(1)) * 0.5;

      // quads.col(0) = (cosines.col(0).array() < 0).select(0.25 * dblA, quads.col(0));
      // quads.col(1) = (cosines.col(0).array() < 0).select(0.125 * dblA, quads.col(1));
      // quads.col(2) = (cosines.col(0).array() < 0).select(0.125 * dblA, quads.col(2));
      if (cosines(0) < 0) {
        quad(0) = 0.25 * dblA;
        quad(1) = 0.125 * dblA;
        quad(2) = 0.125 * dblA;
      }

      // quads.col(0) = (cosines.col(1).array() < 0).select(0.125 * dblA, quads.col(0));
      // quads.col(1) = (cosines.col(1).array() < 0).select(0.25 * dblA, quads.col(1));
      // quads.col(2) = (cosines.col(1).array() < 0).select(0.125 * dblA, quads.col(2));
      if (cosines(1) < 0) {
        quad(0) = 0.125 * dblA;
        quad(1) = 0.25 * dblA;
        quad(2) = 0.125 * dblA;
      }

      // quads.col(0) = (cosines.col(2).array() < 0).select(0.125 * dblA, quads.col(0));
      // quads.col(1) = (cosines.col(2).array() < 0).select(0.125 * dblA, quads.col(1));
      // quads.col(2) = (cosines.col(2).array() < 0).select(0.25 * dblA, quads.col(2));
      if (cosines(2) < 0) {
        quad(0) = 0.125 * dblA;
        quad(1) = 0.125 * dblA;
        quad(2) = 0.25 * dblA;
      }

      massQuads.row(i) = quad;
    },
    tbb::static_partitioner());
}

void SphereICPSolverData::computeNewLaplacianAndMass(const TriMeshGeo &mesh)
{
  // compute cotangent
  edgeCotangents.resize(mesh.numTriangles(), 3);
  massQuads.resize(mesh.numTriangles(), 3);

  computeCotWeights(mesh);

  L.setZero();
  L.resize(mesh.numVertices(), mesh.numVertices());
  L.reserve(30 * mesh.numVertices());

  std::vector<Eigen::Triplet<double>> IJV;
  IJV.reserve(mesh.numTriangles() * 3 * 4);

  edgeWeightDestinations.resize(mesh.numTriangles());

  // Loop over triangles
  for (int i = 0; i < mesh.numTriangles(); i++) {
    // loop over edges of element
    for (Eigen::Index e = 0; e < edgeVertexIDs.rows(); e++) {
      int source = mesh.triVtxID(i, edgeVertexIDs(e, 0));
      int dest = mesh.triVtxID(i, edgeVertexIDs(e, 1));

      IJV.emplace_back(source, dest, edgeCotangents(i, e));
      IJV.emplace_back(dest, source, edgeCotangents(i, e));
      IJV.emplace_back(source, source, -edgeCotangents(i, e));
      IJV.emplace_back(dest, dest, -edgeCotangents(i, e));
    }
  }

  L.setFromTriplets(IJV.begin(), IJV.end());

  ES::EntryMap LEntryMap;
  ES::BuildEntryMap(L, LEntryMap);

  for (int i = 0; i < mesh.numTriangles(); i++) {
    // loop over edges of element
    for (Eigen::Index e = 0; e < edgeVertexIDs.rows(); e++) {
      int source = mesh.triVtxID(i, edgeVertexIDs(e, 0));
      int dest = mesh.triVtxID(i, edgeVertexIDs(e, 1));

      int row_cols[4][2] = {
        { source, dest },
        { dest, source },
        { source, source },
        { dest, dest }
      };

      double val[4] = { 1, 1, -1, -1 };

      for (int j = 0; j < 4; j++) {
        auto iter = LEntryMap.find(std::make_pair(row_cols[j][0], row_cols[j][1]));
        if (iter != LEntryMap.end()) {
          edgeWeightDestinations[i][e][j] = std::make_tuple(iter->second, val[j]);
        }
        else {
          SPDLOG_LOGGER_ERROR(Logger::lgr(), "Cannot find entry");
          abort();
        }
      }
    }
  }

  ES::Expand3(L, L3);
  LToL3Mapping.resize(L3.nonZeros());
  for (int rowi = 0; rowi < L3.rows(); rowi++) {
    for (const int *colIter = L3.innerIndexPtr() + L3.outerIndexPtr()[rowi]; colIter != L3.innerIndexPtr() + L3.outerIndexPtr()[rowi + 1]; ++colIter) {
      int col = *colIter;
      std::ptrdiff_t offset = colIter - L3.innerIndexPtr();

      int Lrow = rowi / 3;
      int Lcol = col / 3;

      auto iter = LEntryMap.find(std::make_pair(Lrow, Lcol));
      if (iter != LEntryMap.end()) {
        LToL3Mapping(offset) = iter->second;
      }
      else {
        SPDLOG_LOGGER_ERROR(Logger::lgr(), "Cannot find entry");
        abort();
      }
    }
  }

  ES::ExpandN(L, L12, 12);
  LToL12Mapping.resize(L12.nonZeros());
  for (int rowi = 0; rowi < L12.rows(); rowi++) {
    for (const int *colIter = L12.innerIndexPtr() + L12.outerIndexPtr()[rowi]; colIter != L12.innerIndexPtr() + L12.outerIndexPtr()[rowi + 1]; ++colIter) {
      int col = *colIter;
      std::ptrdiff_t offset = colIter - L12.innerIndexPtr();

      int Lrow = rowi / 12;
      int Lcol = col / 12;

      auto iter = LEntryMap.find(std::make_pair(Lrow, Lcol));
      if (iter != LEntryMap.end()) {
        LToL12Mapping(offset) = iter->second;
      }
      else {
        SPDLOG_LOGGER_ERROR(Logger::lgr(), "Cannot find entry");
        abort();
      }
    }
  }

  // now mass
  MVec.setZero(mesh.numVertices());
  IJV.clear();
  for (int i = 0; i < mesh.numTriangles(); i++) {
    for (int v = 0; v < 3; v++) {
      int vid = mesh.triVtxID(i, v);
      for (int dof = 0; dof < 3; dof++) {
        IJV.emplace_back(vid * 3 + dof, vid * 3 + dof, massQuads(i, v));
      }
      MVec[vid] += massQuads(i, v);
    }
  }

  M3.resize(mesh.numVertices() * 3, mesh.numVertices() * 3);
  M3.reserve(mesh.numVertices() * 3);
  M3.setFromTriplets(IJV.begin(), IJV.end());

  IJV.clear();
  for (int i = 0; i < mesh.numTriangles(); i++) {
    for (int v = 0; v < 3; v++) {
      int vid = mesh.triVtxID(i, v);
      for (int dof = 0; dof < 12; dof++) {
        IJV.emplace_back(vid * 12 + dof, vid * 12 + dof, massQuads(i, v));
      }
    }
  }

  M12.resize(mesh.numVertices() * 12, mesh.numVertices() * 12);
  M12.reserve(mesh.numVertices() * 12);
  M12.setFromTriplets(IJV.begin(), IJV.end());

  MInvL3 = L3;
  tbb::parallel_for(0, (int)MInvL3.rows(), [&](int rowi) {
    for (ES::SpMatD::InnerIterator it(MInvL3, rowi); it; ++it) {
      it.valueRef() /= MVec(rowi / 3);
    }
  });

  MInvL12 = L12;
  tbb::parallel_for(0, (int)MInvL12.rows(), [&](int rowi) {
    for (ES::SpMatD::InnerIterator it(MInvL12, rowi); it; ++it) {
      it.valueRef() /= MVec(rowi / 12);
    }
  });

  ES::SymbolicMm(L3, MInvL3, biL3, &mulL3, 1);
  ES::SymbolicMm(L12, MInvL12, biL12, &mulL12, 1);
}

void SphereICPSolverData::updateLaplacianAndMass()
{
  // update L, M, Mvec
  std::memset(L.valuePtr(), 0, sizeof(double) * L.nonZeros());
  std::memset(M3.valuePtr(), 0, sizeof(double) * M3.nonZeros());
  std::memset(M12.valuePtr(), 0, sizeof(double) * M12.nonZeros());

  MVec.setZero();

  computeCotWeights(sphereMeshRest);

  // Loop over triangles
  for (int i = 0; i < sphereMeshRest.numTriangles(); i++) {
    // loop over edges of element
    for (Eigen::Index e = 0; e < edgeVertexIDs.rows(); e++) {
      for (int j = 0; j < 4; j++) {
        double w = edgeCotangents(i, e) * std::get<1>(edgeWeightDestinations[i][e][j]);
        L.valuePtr()[std::get<0>(edgeWeightDestinations[i][e][j])] += w;
      }
    }

    for (int v = 0; v < 3; v++) {
      for (int dof = 0; dof < 3; dof++) {
        M3.valuePtr()[sphereMeshRest.triVtxID(i, v) * 3 + dof] += massQuads(i, v);
      }

      for (int dof = 0; dof < 12; dof++) {
        M12.valuePtr()[sphereMeshRest.triVtxID(i, v) * 12 + dof] += massQuads(i, v);
      }

      MVec[sphereMeshRest.triVtxID(i, v)] += massQuads(i, v);
    }
  }

  // update L3
  tbb::parallel_for(0, (int)LToL3Mapping.size(), [&](int ei) {
    L3.valuePtr()[ei] = L.valuePtr()[LToL3Mapping[ei]];
  });

  // update L12
  tbb::parallel_for(0, (int)LToL12Mapping.size(), [&](int ei) {
    L12.valuePtr()[ei] = L.valuePtr()[LToL12Mapping[ei]];
  });

  // update MInvL
  std::memcpy(MInvL3.valuePtr(), L3.valuePtr(), MInvL3.nonZeros() * sizeof(double));
  tbb::parallel_for(0, (int)MInvL3.rows(), [&](int rowi) {
    for (ES::SpMatD::InnerIterator it(MInvL3, rowi); it; ++it) {
      it.valueRef() /= MVec(rowi / 3);
    }
  });

  std::memcpy(MInvL12.valuePtr(), L12.valuePtr(), MInvL12.nonZeros() * sizeof(double));
  tbb::parallel_for(0, (int)MInvL12.rows(), [&](int rowi) {
    for (ES::SpMatD::InnerIterator it(MInvL12, rowi); it; ++it) {
      it.valueRef() /= MVec(rowi / 12);
    }
  });

  // update biL
  ES::Mm(L3, MInvL3, mulL3, biL3, 1);
  ES::Mm(L12, MInvL12, mulL12, biL12, 1);
}

enum class OptEnergyType : int
{
  ICP_ENERGY = 0,
  SMOOTHNESS_ENERGY,
  EXPANSION_ENERGY,
  TRANSLATION_ENERGY,
  INVERSION_ENERGY,
  INERTIA_ENERGY,
  NUM_ENERGIES,
};

}  // namespace MedialAxisRepresentation

SphereICPSolver::SphereICPSolver(const TriMeshGeo &sm, const SphereICPSolverParameters &pm)
{
  sphereTemplateMesh = sm;
  params = pm;

  ES::V3d sphereCenter(0.0, 0.0, 0.0);
  for (int i = 0; i < sphereTemplateMesh.numVertices(); i++) {
    sphereCenter += ES::toESV3d(sphereTemplateMesh.pos(i));
  }
  sphereCenter /= sphereTemplateMesh.numVertices();

  // compute sphere radius
  double radius = 0.0;
  for (int i = 0; i < sphereTemplateMesh.numVertices(); i++) {
    radius += (ES::toESV3d(sphereTemplateMesh.pos(i)) - sphereCenter).norm();
    sphereTemplateMesh.pos(i) -= Vec3d(sphereCenter.data());
  }
  radius /= sphereTemplateMesh.numVertices();
  sphereCenter.setZero();

  int nv = sphereTemplateMesh.numVertices();
  int nv3 = nv * 3;

  if (params.centerFlag == CenterOptimizationFlag::OPTIM) {
    sd = new SphereICPSolverData(nv, sphereTemplateMesh.numTriangles(), int(pm.numPositionDOFs), true, pm.addConstraints);
  }
  else {
    sd = new SphereICPSolverData(nv, sphereTemplateMesh.numTriangles(), int(pm.numPositionDOFs), false, pm.addConstraints);
  }

  // create sphere samples
  while (sd->sphereSamples.size() < 100) {
    int furthestOne = -1;
    double maxDist = 0;
    for (int i = 0; i < nv; i++) {
      double closestSampleDist = 1e100;
      for (int vid : sd->sphereSamples) {
        Vec3d diff = sphereTemplateMesh.pos(i) - sphereTemplateMesh.pos(vid);
        closestSampleDist = std::min(closestSampleDist, len2(diff));
      }
      if (closestSampleDist > maxDist) {
        maxDist = closestSampleDist;
        furthestOne = i;
      }
    }

    sd->sphereSamples.emplace_back(furthestOne);
  }

  TriMeshGeo kk;
  for (int vid : sd->sphereSamples) {
    kk.addPos(sphereTemplateMesh.pos(vid));
  }
  kk.save("zz.obj");

  // build sphere normals
  sd->sphereMeshRest = sphereTemplateMesh;
  sd->sphereMeshCur = sphereTemplateMesh;
  sd->sphereMeshNormals.buildPseudoNormals(sd->sphereMeshCur);

  // get matrix from sphere
  BasicUtilities::triMeshGeoToMatrices(sphereTemplateMesh, sd->sphereRestPostions, sd->sphereTriangles);

  // optimization variables
  sd->qlow.setConstant(-1e20);
  sd->qhi.setConstant(1e20);

  if (pm.addConstraints) {
    sd->g.setZero();
    sd->lambda.setZero();
    sd->clow.setConstant(0);
    sd->chi.setConstant(1e20);
  }

  // vtx area
  std::vector<tbb::spin_mutex> lock(nv);

  sd->vtxAreas.setZero();
  tbb::parallel_for(0, sphereTemplateMesh.numTriangles(), [&](int t) {
    double area_t = getTriangleArea(sphereTemplateMesh.pos(t, 0), sphereTemplateMesh.pos(t, 1), sphereTemplateMesh.pos(t, 2));

    for (int ti = 0; ti < 3; ti++) {
      int v_ti = sphereTemplateMesh.triVtxID(t, ti);

      lock.at(v_ti).lock();
      sd->vtxAreas(v_ti) += area_t / 3.0;
      lock.at(v_ti).unlock();
    }
  });

  sd->vtxAreas /= sd->vtxAreas.maxCoeff();

  // base hessian
  sd->baseHessian.resize(nv * 3, nv * 3);
  std::vector<ES::TripletD> baseHessianEntries;
  baseHessianEntries.reserve(81 * sphereTemplateMesh.numTriangles());
  for (int t = 0; t < sphereTemplateMesh.numTriangles(); t++) {
    for (int vi = 0; vi < 3; vi++) {
      for (int vj = 0; vj < 3; vj++) {
        int vid_i = sphereTemplateMesh.triVtxID(t, vi);
        int vid_j = sphereTemplateMesh.triVtxID(t, vj);
        for (int i = 0; i < 3; i++) {
          for (int j = 0; j < 3; j++) {
            baseHessianEntries.emplace_back(vid_i * 3 + i, vid_j * 3 + j, 0.0);
          }
        }
      }
    }
  }
  sd->baseHessian.setFromTriplets(baseHessianEntries.begin(), baseHessianEntries.end());

  //// bi laplacian matrix
  // biLaplacian(sd->biL, sd->sphereRestPostions, sd->sphereTriangles);

  //// mass matrix
  // massMatrix(sd->M, sd->sphereRestPostions, sd->sphereTriangles);
  //  ES::SpMatD L1 = sd->biL, M1 = sd->M;
  sd->computeNewLaplacianAndMass(sphereTemplateMesh);
  // std::cout << (L1 - sd->biL).norm() << ',' << (M1 - sd->M).norm() << std::endl;

  // x
  for (int vi = 0; vi < nv; vi++) {
    sd->curPosition.segment<3>(vi * 3) = sd->sphereRestPostions.row(vi);
    sd->curPosition1.segment<3>(vi * 3) = sd->sphereRestPostions.row(vi);
    sd->restPosition.segment<3>(vi * 3) = sd->sphereRestPostions.row(vi);
    sd->curVel.segment<3>(vi * 3).setZero();
  }

  // inertial
  params.inertial_w = 1.0;
  if (params.timestep > 0) {
    params.inertial_w = 1.0;
  }
  else {
    params.timestep = 1.0;
    params.inertial_w = 0.0;
  }

  // inertia
  // M 1/h(1/h(x - x0) -v0)
  // M (1/h2 x - 1/h2 x0 -1/h v0)
  if (params.inertial_w > 0) {
    sd->A_inertial = sd->M3 / (params.timestep * params.timestep);
    sd->b_inertial.setZero(nv3);
    sd->inertiaEnergy = std::make_shared<VegaFEM::PredefinedPotentialEnergies::QuadraticPotentialEnergy>(sd->A_inertial, sd->b_inertial);
  }

  // icp energy
  if (params.pnorm == 2) {
    if (params.sampleCount > 1) {
      std::vector<int> sampleTriangleID;
      std::vector<std::array<double, 3>> sampleBaryWeights;
      std::vector<double> sampleWeights;
      sampleTriangle(sphereTemplateMesh, 3, sampleTriangleID, sampleWeights, sampleBaryWeights);

      sd->sampleVertexIndices.resize(sampleTriangleID.size() * 3);
      sd->sampleVertexWeights.resize(sampleTriangleID.size() * 3);
      sd->sampleRestCoeffs.resize(sampleTriangleID.size());

      sd->icpSampleData.coeffs.resize(sampleTriangleID.size());
      sd->icpSampleData.targetPositions.resize(sampleTriangleID.size() * 3);
      sd->icpSampleData.targetNormals.resize(sampleTriangleID.size() * 3);
      for (int si = 0; si < (int)sampleTriangleID.size(); si++) {
        sd->sampleRestCoeffs[si] = sampleWeights[si];

        for (int vi = 0; vi < 3; vi++) {
          sd->sampleVertexIndices[si * 3 + vi] = sphereTemplateMesh.triVtxID(sampleTriangleID[si], vi);
          sd->sampleVertexWeights[si * 3 + vi] = sampleBaryWeights[si][vi];
        }
      }

      TriMeshGeo msample;
      for (int si = 0; si < (int)sampleTriangleID.size(); si++) {
        Vec3d p(0.0);
        for (int vi = 0; vi < 3; vi++) {
          p += sphereTemplateMesh.pos(sampleTriangleID[si], vi) * sampleBaryWeights[si][vi];
        }
        msample.addPos(p);
      }
      msample.save("zzs.obj");

      sd->sampleRestCoeffs /= sd->sampleRestCoeffs.maxCoeff();

      sd->icpEnergyQuadraticDenseSampled = std::make_shared<VegaFEM::PredefinedPotentialEnergies::BarycentricCoordinateSliding>(
        sd->baseHessian, sd->restPosition, (int)sampleTriangleID.size(), sd->icpSampleData.coeffs.data(), sd->icpSampleData.targetPositions.data(), sd->icpSampleData.targetNormals.data(),
        sd->sampleVertexIndices.data(), sd->sampleVertexWeights.data(), 3, 1.0, 1, 1, 0);

      sd->icpEnergy = sd->icpEnergyQuadraticDenseSampled;
    }
    else {
      sd->icpEnergyQuadratic = std::make_shared<VegaFEM::PredefinedPotentialEnergies::MultipleVertexSliding>(
        sd->baseHessian, sd->restPosition, nv, sd->icpData.coeffs.data(), sd->icpData.targetPositions.data(), sd->icpData.targetNormals.data(),
        sd->icpData.vertexIndices.data(), 1.0, 1, 0, 0);
      sd->icpEnergy = sd->icpEnergyQuadratic;
    }
  }
  else {
    // p-order icp energy
    sd->icpEnergyPOrder = std::make_shared<VegaFEM::PredefinedPotentialEnergies::MultipleVertexPullingSoftConstraintsPOrder>(
      sd->baseHessian, sd->restPosition, nv, sd->icpData.coeffs.data(), sd->icpData.targetPositions.data(), sd->icpData.targetNormals.data(),
      sd->icpData.vertexIndices.data(), 1.0, 0, 0, params.pnorm);
    sd->icpEnergy = sd->icpEnergyPOrder;
  }

  // smoothness energy
  // 1/2 (x - x0)^T biL (x - x0)
  // 1/2 xT biL x - x0^T biL x + x0^T BiL x0
  sd->b_smooth.resize(nv3);
  ES::Mv(sd->biL3, sd->restPosition, sd->b_smooth);
  sd->b_smooth *= -1.0;
  sd->smoothnessEnergyXDOF = std::make_shared<VegaFEM::PredefinedPotentialEnergies::QuadraticPotentialEnergy>(sd->biL3, sd->b_smooth);

  sd->smoothnessEnergyTriangleStrainXDOF = std::make_shared<VegaFEM::PredefinedPotentialEnergies::SurfaceTriangleDeformation>(sd->restPosition, sd->sphereMeshCur);
  // bending
  // sd->smoothnessEnergyAbsH = std::make_shared<VegaFEM::PredefinedPotentialEnergies::SurfaceSmoothnessAbsoluteMeanCurvature>(sd->xrest, sd->sphereMeshCur, 1e-8);

#if 0
  VegaFEM::NonlinearOptimization::FiniteDifference fd(VegaFEM::NonlinearOptimization::FiniteDifference::M_FIVE_POINT, 1e-7);
  ES::VXd xtest = sd->restPosition;

  ES::MXd aa;
  ES::MXi bb;
  igl::readOBJ("v000-fitting.obj", aa, bb);
  for (int i = 0; i < sphereTemplateMesh.numVertices(); i++) {
    xtest.segment<3>(i * 3) = aa.row(i);
  }
  std::cout << sd->smoothnessEnergyTriangleStrainXDOF->func(sd->restPosition) << std::endl;

  double err[2];
  fd.testEnergy(sd->smoothnessEnergyTriangleStrainXDOF, true, true, -1.0, xtest.data(), -1, err, err + 1);
  std::cout << err[0] << ',' << err[1] << std::endl;
  std::cout << sd->smoothnessEnergyXDOF->func(xtest) << std::endl;
  std::cout << sd->smoothnessEnergyTriangleStrainXDOF->func(xtest) << std::endl;
  exit(1);
#endif

  // expansion
  if (params.expansionConstraintType == ExpansionConstraintType::QUADRATIC) {
    if (params.sampleCount > 1) {
      sd->expansionSampleData.coeffs = sd->sampleRestCoeffs;
      sd->expansionSampleData.targetPositions.resize(sd->sampleRestCoeffs.size() * 3);
      sd->expansionSampleData.targetNormals.resize(sd->sampleRestCoeffs.size() * 3);

      sd->sphereExpansionEnergyQuadDenseSampled = std::make_shared<VegaFEM::PredefinedPotentialEnergies::BarycentricCoordinateSliding>(
        sd->baseHessian, sd->restPosition, (int)sd->expansionSampleData.coeffs.size(), sd->expansionSampleData.coeffs.data(),
        sd->expansionSampleData.targetPositions.data(), sd->expansionSampleData.targetNormals.data(),
        sd->sampleVertexIndices.data(), sd->sampleVertexWeights.data(), 3, 1.0, 0, 0, 0);

      sd->expansionEnergy = sd->sphereExpansionEnergyQuadDenseSampled;
    }
    else {
      sd->sphereExpansionEnergyQuad = std::make_shared<VegaFEM::PredefinedPotentialEnergies::MultipleVertexSliding>(
        sd->baseHessian, sd->restPosition, nv, sd->expansionData.coeffs.data(), sd->expansionData.targetPositions.data(), sd->expansionData.targetNormals.data(),
        sd->expansionData.vertexIndices.data(), 1.0, 0, 0, 0);

      sd->expansionEnergy = sd->sphereExpansionEnergyQuad;
    }
  }
  else if (params.expansionConstraintType == ExpansionConstraintType::LINEAR) {
    sd->expansionForce.setZero(nv3);
    sd->sphereExpansionEnergyLinear = std::make_shared<VegaFEM::PredefinedPotentialEnergies::LinearPotentialEnergy>(sd->expansionForce);
    sd->expansionEnergy = sd->sphereExpansionEnergyLinear;
  }
  else {
    throw std::domain_error("unknown expansion type");
  }

  // translation
  if (params.centerFlag == CenterOptimizationFlag::OPTIM) {
    sd->sphereTranslationEnergy = std::make_shared<VegaFEM::PredefinedPotentialEnergies::MultipleVertexConstrainedRigidMotion>(
      sd->restPosition, sd->sphereSamples, nullptr, 0.0, 1.0, true, false);

    sd->A_centerReg.resize(3, 3);
    sd->A_centerReg.setIdentity();
    sd->sphereCenterRegEnergy = std::make_shared<VegaFEM::PredefinedPotentialEnergies::QuadraticPotentialEnergy>(sd->A_centerReg);

    std::vector<int> regDOFs = { nv3, nv3 + 1, nv3 + 2 };
    sd->sphereCenterRegEnergy->setDOFs(regDOFs);
  }
  else {
    sd->sphereTranslationEnergy = std::make_shared<VegaFEM::PredefinedPotentialEnergies::MultipleVertexConstrainedRigidMotion>(
      sd->restPosition, sd->sphereSamples, nullptr, 0.0, 1.0, false, false);
  }

  // VegaFEM::NonlinearOptimization::FiniteDifference fd(VegaFEM::NonlinearOptimization::FiniteDifference::M_FIVE_POINT, 1e-7);
  // ES::VXd xtest(sd->x.size() + 12);
  // xtest.head(sd->x.size()) = sd->x;
  // xtest.segment<3>(sd->x.size()).setConstant(0.1);
  // ES::Mp<ES::M3d>(xtest.data() + sd->x.size() + 3).setIdentity();

  // ES::M3d R;
  // R = Eigen::AngleAxisd(M_PI / 6.0, ES::V3d(1, 1, 1).normalized());
  // for (int i = 0; i < sd->sphereMeshCur.numVertices(); i++) {
  //   xtest.segment<3>(i * 3) = R * (sd->x.segment<3>(i * 3) - sphereCenter) * 1.5 + sphereCenter;
  // }

  // double err[2];
  // fd.testEnergy(sd->sphereTranslationEnergy, true, true, -1.0, xtest.data(), -1, err, err + 1);
  // std::cout << err[0] << ',' << err[1] << std::endl;
  // exit(1);

  // barrier
  sd->triangleRestVols.resize(sphereTemplateMesh.numTriangles());
  sd->triangleRestAreas.resize(sphereTemplateMesh.numTriangles());
  tbb::parallel_for(0, sphereTemplateMesh.numTriangles(), [&](int t) {
    ES::M3d M;
    for (int v = 0; v < 3; v++) {
      M.col(v) = ES::toESV3d(sphereTemplateMesh.pos(t, v)) - sphereCenter;
    }
    sd->triangleRestVols[t] = VegaFEM::NonlinearOptimization::Determinant::Dim3::det(M.data());
    sd->triangleRestAreas[t] = getTriangleArea(sphereTemplateMesh.pos(t, 0), sphereTemplateMesh.pos(t, 1), sphereTemplateMesh.pos(t, 2));
  });

  double thres = 5e-2;
  if (params.centerFlag == CenterOptimizationFlag::OPTIM) {
    sd->inversionBarrierEnergy = std::make_shared<SphereInversionBarrierEnergy>(sphereCenter, sphereTemplateMesh,
      sd->triangleRestVols, sd->triangleRestAreas, radius, true, thres, thres, thres, params.triangleScale);
  }
  else {
    sd->inversionBarrierEnergy = std::make_shared<SphereInversionBarrierEnergy>(sphereCenter, sphereTemplateMesh,
      sd->triangleRestVols, sd->triangleRestAreas, radius, false, thres, thres, thres, params.triangleScale);
  }

  // VegaFEM::NonlinearOptimization::FiniteDifference fd(VegaFEM::NonlinearOptimization::FiniteDifference::M_FIVE_POINT, 1e-7);
  // ES::VXd xtest(sd->x.size() + 3);
  // xtest.head(sd->x.size()) = sd->x;
  // xtest.tail<3>() = sphereCenter + ES::V3d(1e-3, 1e-3, 1e-3);

  // for (int i = 0; i < sd->sphereMeshCur.numVertices(); i++) {
  //   xtest.segment<3>(i * 3) = (sd->x.segment<3>(i * 3) - sphereCenter) * 1.3 + sphereCenter;
  // }
  // double err[2];
  // fd.testEnergy(sd->inversionBarrierEnergy, true, true, -1.0, xtest.data(), -1, err, err + 1);
  // std::cout << err[0] << ',' << err[1] << std::endl;
  // std::cout << sd->inversionBarrierEnergy->func(xtest) << std::endl;
  // exit(1);

  if (params.centerFlag == CenterOptimizationFlag::OPTIM) {
    sd->energiesXDOFs = std::make_shared<VegaFEM::NonlinearOptimization::PotentialEnergies>(nv3 + 12);
  }
  else {
    sd->energiesXDOFs = std::make_shared<VegaFEM::NonlinearOptimization::PotentialEnergies>(nv3);
  }

  sd->energiesXDOFs->addPotentialEnergy(sd->icpEnergy, params.icp_w);
  sd->energiesXDOFs->addPotentialEnergy(sd->smoothnessEnergyXDOF, params.smooth_w_init);
  sd->energiesXDOFs->addPotentialEnergy(sd->expansionEnergy, params.expansion_w);
  sd->energiesXDOFs->addPotentialEnergy(sd->sphereTranslationEnergy, params.fixTranslation_w);
  sd->energiesXDOFs->addPotentialEnergy(sd->inversionBarrierEnergy, params.barrier_w);

  if (params.inertial_w > 0) {
    sd->energiesXDOFs->addPotentialEnergy(sd->inertiaEnergy, params.inertial_w);
  }

  if (params.centerFlag == CenterOptimizationFlag::OPTIM)
    sd->energiesXDOFs->addPotentialEnergy(sd->sphereCenterRegEnergy, 0.1);

  // triangle quality improvement
  sd->energiesXDOFs->addPotentialEnergy(sd->smoothnessEnergyTriangleStrainXDOF, 0.1);

  sd->energiesXDOFs->init();

  if (params.addConstraints) {
    sd->constraints = std::make_shared<SphereTriangleFlippingConstraints>(sphereCenter, sd->sphereTriangles, nv);
  }

  if (params.numPositionDOFs == NumPositionDOFs::TWELVE) {
    if (params.centerFlag == CenterOptimizationFlag::OPTIM) {
      sd->affineSpaceEnergies = std::make_shared<VegaFEM::PredefinedPotentialEnergies::VertexAffineToPositionEnergy>(nv * 12 + 12, sd->restPosition, sd->energiesXDOFs);
      sd->energiesAffineDOFs = std::make_shared<VegaFEM::NonlinearOptimization::PotentialEnergies>(nv * 12 + 12);
    }
    else {
      sd->affineSpaceEnergies = std::make_shared<VegaFEM::PredefinedPotentialEnergies::VertexAffineToPositionEnergy>(nv * 12, sd->restPosition, sd->energiesXDOFs);
      sd->energiesAffineDOFs = std::make_shared<VegaFEM::NonlinearOptimization::PotentialEnergies>(nv * 12);
    }

    // smoothness energy
    sd->smoothnessEnergyAffineDOF = std::make_shared<VegaFEM::PredefinedPotentialEnergies::QuadraticPotentialEnergy>(sd->biL12);

    sd->energiesAffineDOFs->addPotentialEnergy(sd->affineSpaceEnergies, 1.0);
    sd->energiesAffineDOFs->addPotentialEnergy(sd->smoothnessEnergyAffineDOF, params.smooth_w_init);
    sd->energiesAffineDOFs->init();

    sd->energiesXDOFs->setEnergyCoeffs(int(OptEnergyType::SMOOTHNESS_ENERGY), 0.0);

    // exit(1);
  }

#if defined(USE_KNITRO)
  if (params.numPositionDOFs == NumPositionDOFs::TWELVE) {
    sd->icpSolver.init(sd->energiesAffineDOFs, sd->constraints, sd->q, sd->qlow, sd->qhi, sd->clow, sd->chi);
  }
  else if (params.numPositionDOFs == NumPositionDOFs::THREE) {
    sd->icpSolver.init(sd->energiesXDOFs, sd->constraints, sd->q, sd->qlow, sd->qhi, sd->clow, sd->chi);
  }
  else {
    throw std::invalid_argument("unsupported num position dofs");
  }

  sd->icpSolver.solver->setVerbose(params.solverVerbose);
  sd->icpSolver.solver->setMaxIter(params.solverNumIter);
#endif
}

SphereICPSolver::~SphereICPSolver()
{
  delete sd;
}

int SphereICPSolver::icp(const double center[3], double radius,
  const TriMeshGeo &targetMesh, const TriMeshBVTree &targetMeshBVTree, const TriMeshPseudoNormal &targetMeshNormals,
  const std::string &prefix, TriMeshGeo &meshOut)
{
  int dumpMesh = 0;
  if (prefix.length()) {
    dumpMesh = 1;
  }

  int nv = sphereTemplateMesh.numVertices();
  int nv3 = nv * 3;
  int ntri = sphereTemplateMesh.numTriangles();

  sd->sphereMeshRest = sphereTemplateMesh;

  ES::V3d sphereCenter(center[0], center[1], center[2]);
  for (int i = 0; i < nv; i++) {
    sd->sphereMeshRest.pos(i) = norm(sphereTemplateMesh.pos(i)) * radius + Vec3d(center);
  }

  if (dumpMesh) {
    sd->sphereMeshRest.save(fmt::format("{}-si.obj", prefix));
  }

  // initial sphere volume
  double initialVol = 4.0 / 3.0 * M_PI * radius * radius * radius;
  double initialArea = 4.0 * M_PI * radius * radius;

  // update rest positions
  for (int vi = 0; vi < nv; vi++) {
    sd->sphereRestPostions.row(vi) = ES::toESV3d(sd->sphereMeshRest.pos(vi));
  }

  //// bi laplacian matrix
  // sd->biL.setZero();
  // biLaplacian(sd->biL, sd->sphereRestPostions, sd->sphereTriangles);

  //// mass matrix
  // sd->M.setZero();
  // massMatrix(sd->M, sd->sphereRestPostions, sd->sphereTriangles);
  // ES::SpMatD L1 = sd->biL, M1 = sd->M;
  sd->updateLaplacianAndMass();
  // std::cout << (L1 - sd->biL).norm() << ',' << (M1 - sd->M).norm() << std::endl;

  // normals
  sd->sphereMeshCur = sd->sphereMeshRest;
  sd->sphereMeshNormals.updateVertexPositions(sd->sphereMeshCur);

  // update lower bound
  if (params.addConstraints) {
    double lb = std::max(initialVol / ntri * 1e-3, 1e-10);

    sd->clow.setConstant(lb);
    sd->icpSolver.solver->setcRange(sd->clow.data(), sd->chi.data());
  }

  // x
  // for (int vi = 0; vi < nv; vi++) {
  tbb::parallel_for(0, nv, [&](int vi) {
    sd->restPosition.segment<3>(vi * 3) = sd->sphereRestPostions.row(vi);
    sd->curPosition.segment<3>(vi * 3) = sd->sphereRestPostions.row(vi);
    sd->curPosition1.segment<3>(vi * 3) = sd->sphereRestPostions.row(vi);
    sd->curVel.segment<3>(vi * 3).setZero();

    if (params.numPositionDOFs == NumPositionDOFs::THREE) {
      sd->q.segment<3>(vi * 3) = sd->restPosition.segment<3>(vi * 3);
    }
    else if (params.numPositionDOFs == NumPositionDOFs::TWELVE) {
      Eigen::Matrix<double, 3, 4> A;
      A.block<3, 3>(0, 0).setIdentity();
      A.col(3).setZero();

      sd->q.segment<12>(vi * 12) = ES::Mp<ES::V12d>(A.data());
    }
  });

  if (params.inertial_w > 0) {
    // update A_inertia
    // M 1/h(1/h(x - x0) -v0)
    // M (1/h2 x - 1/h2 x0 -1/h v0)
    ES::Mp<ES::VXd>(sd->A_inertial.valuePtr(), sd->A_inertial.nonZeros()) = ES::Mp<const ES::VXd>(sd->M3.valuePtr(), sd->M3.nonZeros()) / (params.timestep * params.timestep);
  }

  // update smoothness
  // 1/2 (x - x0)^T biL (x - x0)
  // 1/2 xT biL x - x0^T biL x + x0^T BiL x0
  ES::Mv(sd->biL3, sd->restPosition, sd->b_smooth);
  sd->b_smooth *= -1.0;

  if (sd->smoothnessEnergyAbsHXDOF)
    sd->smoothnessEnergyAbsHXDOF->updateRestInfo();

  if (sd->affineSpaceEnergies)
    sd->affineSpaceEnergies->computeAToX();

  // update translation
  sd->sphereTranslationEnergy->setRestCenter();
  ES::V3d restCenter = sd->sphereTranslationEnergy->getRestCenter();

  // update barrier
  tbb::parallel_for(0, ntri, [&](int t) {
    ES::M3d M;
    for (int v = 0; v < 3; v++) {
      M.col(v) = ES::toESV3d(sd->sphereMeshCur.pos(t, v)) - sphereCenter;
    }
    sd->triangleRestVols[t] = VegaFEM::NonlinearOptimization::Determinant::Dim3::det(M.data());
    sd->triangleRestAreas[t] = getTriangleArea(sd->sphereMeshCur.pos(t, 0), sd->sphereMeshCur.pos(t, 1), sd->sphereMeshCur.pos(t, 2));
  });
  sd->inversionBarrierEnergy->setRadius(radius);
  sd->inversionBarrierEnergy->setCenter(center);

  // VegaFEM::NonlinearOptimization::FiniteDifference fd(VegaFEM::NonlinearOptimization::FiniteDifference::M_FIVE_POINT, 1e-7);
  // ES::VXd xtest = x;
  // for (int i = 0; i < sphereMeshCur.numVertices(); i++) {
  //   xtest.segment<3>(i * 3) = (x.segment<3>(i * 3) - sphereCenter) * 1e-2 + sphereCenter;
  // }
  // double err[2];
  // fd.testEnergy(barrierEnergy, true, true, -1.0, xtest.data(), -1, err, err + 1);
  // std::cout << err[0] << ',' << err[1] << std::endl;
  // std::cout << barrierEnergy->func(xtest) << std::endl;
  // exit(1);

  sd->energiesXDOFs->setEnergyCoeffs(int(OptEnergyType::ICP_ENERGY), params.icp_w);
  sd->energiesXDOFs->setEnergyCoeffs(int(OptEnergyType::SMOOTHNESS_ENERGY), params.smooth_w_init);
  sd->energiesXDOFs->setEnergyCoeffs(int(OptEnergyType::EXPANSION_ENERGY), params.expansion_w);
  sd->energiesXDOFs->setEnergyCoeffs(int(OptEnergyType::TRANSLATION_ENERGY), params.fixTranslation_w);
  sd->energiesXDOFs->setEnergyCoeffs(int(OptEnergyType::INVERSION_ENERGY), params.barrier_w);

  if (params.numPositionDOFs == NumPositionDOFs::TWELVE) {
    sd->energiesXDOFs->setEnergyCoeffs(int(OptEnergyType::SMOOTHNESS_ENERGY), 0.0);
    sd->energiesAffineDOFs->setEnergyCoeffs(1, params.smooth_w_init);
  }

  if (params.inertial_w > 0) {
    sd->energiesXDOFs->setEnergyCoeffs(int(OptEnergyType::INERTIA_ENERGY), params.inertial_w);
  }

  int numDynamicIterations = int(params.numIter * 0.8);
  for (int iter = 1; iter <= params.numIter; iter++) {
    double ratio = (double)(iter - 1) / (params.numWeightIter - 1);
    ratio = std::clamp(ratio, 0.0, 1.0);

    double w_cur = params.smooth_w_init * (1.0 - ratio) + params.smooth_w_final * ratio;
    w_cur = std::pow(10, w_cur);
    // double w_cur = std::exp(std::log(params.smooth_w_init) * (1 - ratio) + std::log(params.smooth_w_final) * ratio);

    fmt::print("Iter: {}\n", iter);
    fmt::print("Alpha: {}\n", w_cur);

    // compute icp constraints and distances
    // tbb::parallel_for(0, nv, [&](int i) {
    for (int i = 0; i < nv; i++) {
      Vec3d sp_qi = sd->sphereMeshCur.pos(i);
      auto ret = targetMeshBVTree.closestTriangleQuery(targetMesh, sp_qi);
      Vec3d dir = ret.closestPosition - sp_qi;
      Vec3d target_n = targetMeshNormals.getPseudoNormal(targetMesh.triangles().data(), ret.triID, ret.feature);
      Vec3d source_n = sd->sphereMeshNormals.vtxNormal(i);
      double dist = len(dir);
      sd->distAll[i] = dist;

      double cosAngle = std::abs(dot(target_n, source_n));
      if (i == 3990) {
        std::cout << cosAngle << std::endl;
      }

      // double distToCenter = len(ret.closestPosition - Vec3d(center));
      if (cosAngle >= std::cos(90.0 / 180.0 * M_PI)) {
        // && distToCenter <= radius * 2
        // || dot(target_n, trans) < 0.0) {  // && dist <= 0.8 * r) { // smaller than 60 degrees, cos(60), or outside of the mesh
        double angleWeight = 1.0;
        double c1 = M_PI / 3, c2 = M_PI / 2;
        if (cosAngle < std::cos(c1)) {
          double angle = std::acos(cosAngle);
          angleWeight = 1 - std::pow(angle - c1, 2.0) / std::pow(c2 - c1, 2.0);
        }
        angleWeight = std::clamp(angleWeight, 0.0, 1.0);

        double distWeight = 1.0;
        // double c3 = radius * 1.6, c4 = radius * 2;
        // if (distToCenter > c3) {
        //   distWeight = 1 - std::pow(distToCenter - c3, 2.0) / std::pow(c4 - c3, 2.0);
        // }
        // distWeight = std::clamp(distWeight, 0.0, 1.0);

        sd->icpData.coeffs(i) = sd->vtxAreas[i] * angleWeight * distWeight;
        sd->expansionData.coeffs(i) = 0;
        sd->distAll[i] = dist;
      }
      else {
        sd->icpData.coeffs(i) = 0.0;
        sd->expansionData.coeffs(i) = sd->vtxAreas[i];
        sd->distAll[i] = 0;
      }

      sd->icpData.targetPositions.segment<3>(i * 3) = ES::toESV3d(ret.closestPosition);
      sd->icpData.targetNormals.segment<3>(i * 3) = ES::toESV3d(target_n);
    }  // );

    // compute expansion constraints
    double dMax = sd->distAll.maxCoeff();
    double dMid = 0.0;
    int dcount = 0;

    for (int i = 0; i < (int)sd->distAll.size(); i++) {
      if (sd->distAll[i] > 0) {
        dMid += sd->distAll[i];
        dcount += 1;
      }
    }
    dMid /= dcount;

    tbb::parallel_for(0, nv, [&](int i) {
      // for(int i = 0; i < num_v; i++) {
      Vec3d sp_qi = sd->sphereMeshCur.pos(i);
      // expansion
      Vec3d ldiff = sp_qi - Vec3d(sphereCenter.data());
      double ldiff_len = len(ldiff);
      sd->expansionData.targetNormals.segment<3>(i * 3) = ES::toESV3d(ldiff / ldiff_len);

      double defaultIncLength = std::min(ldiff_len * 0.1, 0.02);
      double icpIncLength = dMid;
      double incLength = std::max(icpIncLength, defaultIncLength);  //      *(1 - ratio);
      // double inc_len = std::max(dMax, std::min(ldiff_len * 0.1, 0.02)) * (1 - ratio);

      // double finalLength = std::max(ldiff_len + inc_len, radius);
      double finalLength = incLength + ldiff_len;

      sd->expansionData.targetPositions.segment<3>(i * 3) = sd->expansionData.targetNormals.segment<3>(i * 3) * finalLength + sphereCenter;

      if (params.expansionConstraintType == ExpansionConstraintType::LINEAR) {
        sd->expansionForce.segment<3>(i * 3) = -sd->expansionData.targetNormals.segment<3>(i * 3) * sd->vtxAreas[i];
      }
    });

    if (0) {
      // filtering constraints
      // sd->x1.noalias() = sd->x;
      // for (int i = 0; i < sd->sphereMeshCur.numVertices(); i++) {
      //  if (sd->icpData.coeffs[i] > 0)
      //    sd->x1.segment<3>(i * 3) = sd->icpData.targetPositions.segment<3>(i * 3);
      //}
      // sd->constraints->func(sd->x1, sd->g);

      // sd->invalidConstraintsSet.clear();
      // for (int i = 0; i < (int)sd->g.size(); i++) {
      //   if (sd->g[i] < sd->clow[0]) {
      //     for (int j = 0; j < 3; j++) {
      //       int vid = sd->sphereMeshCur.triVtxID(i, j);
      //       sd->icpData.coeffs[vid] = sd->vtxAreas[vid] * 0.0;
      //       sd->expansionData.coeffs[vid] = sd->vtxAreas[vid];
      //       sd->invalidConstraintsSet.emplace(vid);
      //     }
      //   }
      // }
      // fmt::print("drop: {:4d}/{:d}\n", sd->invalidConstraintsSet.size(), sd->sphereMeshCur.numVertices());
    }

    if (dumpMesh) {
      sd->dumpConstraints().save(fmt::format("{}-exp-pair{:02d}.obj", prefix, iter));
      LibIGLInterface::saveOffMesh(fmt::format("{}-pair-w{:02d}.off", prefix, iter).c_str(), sd->sphereMeshCur, sd->icpData.coeffs.data());
    }

    // update energies
    if (sd->icpEnergyQuadratic) {
      sd->icpEnergyQuadratic->setTargetPos(sd->icpData.targetPositions.data());
      sd->icpEnergyQuadratic->setNormals(sd->icpData.targetNormals.data());
      sd->icpEnergyQuadratic->setCoeff(sd->icpData.coeffs.data());
      sd->icpEnergyQuadratic->computeHessian();
    }

    if (sd->icpEnergyPOrder) {
      sd->icpEnergyPOrder->setTargetPos(sd->icpData.targetPositions.data());
      sd->icpEnergyPOrder->setNormals(sd->icpData.targetNormals.data());
      sd->icpEnergyPOrder->setCoeff(sd->icpData.coeffs.data());

      // VegaFEM::NonlinearOptimization::FiniteDifference fd(VegaFEM::NonlinearOptimization::FiniteDifference::M_FIVE_POINT, 1e-7);
      // double err[2];
      // fd.testEnergy(sd->icpEnergyPOrder, true, true, -1.0, sd->x.data(), -1, err, err + 1);
      // std::cout << err[0] << ',' << err[1] << std::endl;
      // std::cout << barrierEnergy->func(xtest) << std::endl;
      // exit(1);
    }

    // expansion
    if (params.expansionConstraintType == ExpansionConstraintType::QUADRATIC) {
      sd->sphereExpansionEnergyQuad->setTargetPos(sd->expansionData.targetPositions.data());
      sd->sphereExpansionEnergyQuad->setNormals(sd->expansionData.targetNormals.data());
      sd->sphereExpansionEnergyQuad->setCoeff(sd->expansionData.coeffs.data());
      sd->sphereExpansionEnergyQuad->computeHessian();
    }
    else if (params.expansionConstraintType == ExpansionConstraintType::LINEAR) {
      // linear one
      double curVol = 0;
      for (int ti = 0; ti < ntri; ti++) {
        ES::M3d M;
        M.col(0) = ES::toESV3d(sd->sphereMeshCur.pos(ti, 0)) - sphereCenter;
        M.col(1) = ES::toESV3d(sd->sphereMeshCur.pos(ti, 1)) - sphereCenter;
        M.col(2) = ES::toESV3d(sd->sphereMeshCur.pos(ti, 2)) - sphereCenter;
        curVol += std::abs(VegaFEM::NonlinearOptimization::Determinant::Dim3::det(M.data()));
      }
      curVol /= 6.0;
      sd->expansionForce *= std::min(initialVol / curVol, 1.0);
    }

    // update center
    sd->sphereTranslationEnergy->setCenter(restCenter.data());

    if (params.addConstraints) {
      double lb = std::max(initialVol / sd->sphereMeshRest.numTriangles() * 1e-2, 1e-10);

      sd->clow.setConstant(lb);
      sd->constraints->setCenter(center);
      sd->icpSolver.solver->setcRange(sd->clow.data(), sd->chi.data());
    }

    // update inertial
    if (params.inertial_w > 0) {
      // M (1/h2 x - 1/h2 x0 -1/h v0)
      sd->temp.noalias() = sd->curPosition / (params.timestep * params.timestep);
      sd->temp.noalias() += sd->curVel / params.timestep;
      sd->temp *= -1.0;

      ES::Mv(sd->M3, sd->temp, sd->b_inertial);

      // dynamics weights
      if (iter > numDynamicIterations) {
        sd->energiesXDOFs->setEnergyCoeffs(int(OptEnergyType::INERTIA_ENERGY), 0);
      }
    }

    // smoothness weight
    if (params.numPositionDOFs == NumPositionDOFs::THREE) {
      sd->energiesXDOFs->setEnergyCoeffs(int(OptEnergyType::SMOOTHNESS_ENERGY), w_cur);
    }
    else if (params.numPositionDOFs == NumPositionDOFs::TWELVE) {
      sd->energiesXDOFs->setEnergyCoeffs(int(OptEnergyType::SMOOTHNESS_ENERGY), 0);
      sd->energiesAffineDOFs->setEnergyCoeffs(1, w_cur);
    }

    // sd->energies->setEnergyCoeffs(int(OptEnergyType::INERTIA_ENERGY), 0);
    // energies->setEnergyCoeffs(int(EnergyType::EXPANSION_ENERGY), 0);

    if (params.centerFlag == CenterOptimizationFlag::OPTIM) {
      if (params.numPositionDOFs == NumPositionDOFs::THREE) {
        sd->q.segment<3>(nv3) = sphereCenter;
        // sd->xlow.segment<3>(nv3) = sphereCenter;
        // sd->xhi.segment<3>(nv3) = sphereCenter;

        ES::Mp<ES::M3d>(sd->q.data() + nv3 + 3).setIdentity();
        ES::Mp<ES::M3d>(sd->qlow.data() + nv3 + 3).setIdentity();
        ES::Mp<ES::M3d>(sd->qhi.data() + nv3 + 3).setIdentity();
      }
      else if (params.numPositionDOFs == NumPositionDOFs::TWELVE) {
        sd->q.segment<3>(nv * 12) = sphereCenter;

        ES::Mp<ES::M3d>(sd->q.data() + nv * 12 + 3).setIdentity();
        ES::Mp<ES::M3d>(sd->qlow.data() + nv * 12 + 3).setIdentity();
        ES::Mp<ES::M3d>(sd->qhi.data() + nv * 12 + 3).setIdentity();
      }
    }

    std::cout << "Before Opt:\n";
    if (params.numPositionDOFs == NumPositionDOFs::THREE) {
      sd->energiesXDOFs->printEnergy(sd->q);
    }
    else if (params.numPositionDOFs == NumPositionDOFs::TWELVE) {
      sd->energiesXDOFs->printEnergy(sd->curPosition);
      sd->energiesAffineDOFs->printEnergy(sd->q);

      // sd->energiesXDOFs->setEnergyCoeffs(int(OptEnergyType::SMOOTHNESS_ENERGY), 0.0);

      // VegaFEM::NonlinearOptimization::FiniteDifference fd(VegaFEM::NonlinearOptimization::FiniteDifference::M_FIVE_POINT, 1e-7);
      // double err[2];
      // fd.testEnergy(sd->affineSpaceEnergies, true, true, -1.0, sd->q.data(), -1, err, err + 1);
      // std::cout << err[0] << ',' << err[1] << std::endl;
      // exit(1);
    }

#if defined(USE_KNITRO)
    sd->curPosition1.noalias() = sd->curPosition;

    if (params.centerFlag == CenterOptimizationFlag::OPTIM) {
      sd->icpSolver.solver->setxRange(sd->qlow.data(), sd->qhi.data());
    }

    sd->icpSolver.solver->setxInit(sd->q.data());
    int solverRet = sd->icpSolver.solver->solve();
    sd->icpSolver.getx(sd->q);

    if (solverRet != 0) {
      std::cerr << "Solver non zero ret: " << solverRet << std::endl;
      /*meshOut = sd->sphereMeshCur;

      if (dumpMesh) {
        sd->sphereMeshCur.save(fmt::format("{}-icp-mesh{:02d}.obj", prefix, iter));
      }
      return 1;*/
    }
#endif

    std::cout << "After Opt:\n";
    if (params.numPositionDOFs == NumPositionDOFs::THREE) {
      sd->energiesXDOFs->printEnergy(sd->q);

      sd->curPosition.noalias() = sd->q.head(nv3);
    }
    else if (params.numPositionDOFs == NumPositionDOFs::TWELVE) {
      tbb::parallel_for(
        0, nv, [&](int vi) {
          Eigen::Matrix<double, 3, 4> A = ES::Mp<Eigen::Matrix<double, 3, 4>>(sd->q.data() + vi * 12);
          ES::V3d x0 = sd->restPosition.segment<3>(vi * 3);

          sd->curPosition.segment<3>(vi * 3) = A.block<3, 3>(0, 0) * x0 + A.col(3);
        },
        tbb::static_partitioner());

      sd->energiesXDOFs->printEnergy(sd->curPosition);
      sd->energiesAffineDOFs->printEnergy(sd->q);
    }

    sd->curVel.noalias() = (sd->curPosition - sd->curPosition1) / params.timestep;

    if (params.centerFlag == CenterOptimizationFlag::OPTIM) {
      std::cout << "Old center: " << sphereCenter.transpose() << std::endl;
      if (params.numPositionDOFs == NumPositionDOFs::THREE) {
        std::cout << "New center: " << sd->q.segment<3>(nv3).transpose() << std::endl;
      }
      else if (params.numPositionDOFs == NumPositionDOFs::TWELVE) {
        std::cout << "New center: " << sd->q.segment<3>(nv * 12).transpose() << std::endl;
      }
    }

    for (int i = 0; i < sd->sphereMeshCur.numVertices(); i++) {
      sd->sphereMeshCur.pos(i) = Vec3d(sd->curPosition.data() + i * 3);
    }

    sd->sphereMeshNormals.updateVertexPositions(sd->sphereMeshCur);

    if (dumpMesh) {
      sd->sphereMeshCur.save(fmt::format("{}-icp-mesh{:02d}.obj", prefix, iter));
    }

    fmt::print("Iter {} done.\n", iter);
    // aaa
  }

  meshOut = sd->sphereMeshCur;

  return 0;
}

int SphereICPSolver::sim(const double center[3], double radius,
  const TriMeshGeo &targetMesh, const TriMeshBVTree &targetMeshBVTree, const TriMeshPseudoNormal &targetMeshNormals,
  const std::string &prefix, TriMeshGeo &meshOut, Vec3d *updatedCenter)
{
  int dumpMesh = 0;
  if (prefix.length()) {
    dumpMesh = 1;

    if (!std::filesystem::exists(prefix))
      std::filesystem::create_directories(prefix);
  }

  BoundingBox bb(targetMesh.positions());
  Vec3d side = bb.sides();
  double refLength = std::max({ side[0], side[1], side[2] });
  double expansionMaxStepSize = std::min(refLength * params.expansionRelMaxStepSize, radius * 0.5);

  int nv = sphereTemplateMesh.numVertices();
  int nv3 = nv * 3;
  int ntri = sphereTemplateMesh.numTriangles();

  sd->sphereMeshRest = sphereTemplateMesh;

  ES::V3d sphereCenter(center[0], center[1], center[2]);
  for (int i = 0; i < nv; i++) {
    sd->sphereMeshRest.pos(i) = norm(sphereTemplateMesh.pos(i)) * radius + Vec3d(center);
  }

  if (dumpMesh) {
    sd->sphereMeshRest.save(fmt::format("{}/si.obj", prefix));
  }

  // initial sphere volume
  double initialVol = 4.0 / 3.0 * M_PI * radius * radius * radius;
  double initialArea = 4.0 * M_PI * radius * radius;

  // update rest positions
  for (int vi = 0; vi < nv; vi++) {
    sd->sphereRestPostions.row(vi) = ES::toESV3d(sd->sphereMeshRest.pos(vi));
  }

  //// bi laplacian matrix
  // sd->biL.setZero();
  // biLaplacian(sd->biL, sd->sphereRestPostions, sd->sphereTriangles);

  //// mass matrix
  // sd->M.setZero();
  // massMatrix(sd->M, sd->sphereRestPostions, sd->sphereTriangles);
  // ES::SpMatD L1 = sd->biL, M1 = sd->M;
  sd->updateLaplacianAndMass();
  // std::cout << (L1 - sd->biL).norm() << ',' << (M1 - sd->M).norm() << std::endl;

  // normals
  sd->sphereMeshCur = sd->sphereMeshRest;
  sd->sphereMeshNormals.updateVertexPositions(sd->sphereMeshCur);

  // update lower bound
  if (params.addConstraints) {
    double lb = std::max(initialVol / ntri * 1e-3, 1e-10);

    sd->clow.setConstant(lb);
    sd->icpSolver.solver->setcRange(sd->clow.data(), sd->chi.data());
  }

  // x
  // for (int vi = 0; vi < nv; vi++) {
  tbb::parallel_for(0, nv, [&](int vi) {
    sd->restPosition.segment<3>(vi * 3) = sd->sphereRestPostions.row(vi);
    sd->curPosition.segment<3>(vi * 3) = sd->sphereRestPostions.row(vi);
    sd->curPosition1.segment<3>(vi * 3) = sd->sphereRestPostions.row(vi);
    sd->curVel.segment<3>(vi * 3).setZero();

    if (params.numPositionDOFs == NumPositionDOFs::THREE) {
      sd->q.segment<3>(vi * 3) = sd->restPosition.segment<3>(vi * 3);
    }
    else if (params.numPositionDOFs == NumPositionDOFs::TWELVE) {
      Eigen::Matrix<double, 3, 4> A;
      A.block<3, 3>(0, 0).setIdentity();
      A.col(3).setZero();

      sd->q.segment<12>(vi * 12) = ES::Mp<ES::V12d>(A.data());
    }
  });

  double dm = 0.01;
  if (params.inertial_w > 0) {
    // update A_inertia
    // M acc + D vel
    // M 1/h(v-v0) + D 1/h (x - x0)
    // M 1/h(1/h(x - x0) -v0) + aM 1/h(x - x0)
    // M (1/h2 x - 1/h2 x0 -1/h v0)
    // (1/h2 M + a/h M) x - M (1/h2 x0 + 1/h v0) - a/h M x0
    double s = 1.0 / (params.timestep * params.timestep) + dm / params.timestep;
    ES::Mp<ES::VXd>(sd->A_inertial.valuePtr(), sd->A_inertial.nonZeros()) = ES::Mp<const ES::VXd>(sd->M3.valuePtr(), sd->M3.nonZeros()) * s;
  }

  // update smoothness
  // 1/2 (x - x0)^T biL (x - x0)
  // 1/2 xT biL x - x0^T biL x + x0^T BiL x0
  ES::Mv(sd->biL3, sd->restPosition, sd->b_smooth);
  sd->b_smooth *= -1.0;

  if (sd->smoothnessEnergyAbsHXDOF)
    sd->smoothnessEnergyAbsHXDOF->updateRestInfo();

  if (sd->smoothnessEnergyTriangleStrainXDOF)
    sd->smoothnessEnergyTriangleStrainXDOF->updateRestInfo();

  if (sd->affineSpaceEnergies)
    sd->affineSpaceEnergies->computeAToX();

  // update translation
  sd->sphereTranslationEnergy->setRestCenter();
  ES::V3d restCenter = sd->sphereTranslationEnergy->getRestCenter();

  // update barrier
  tbb::parallel_for(0, ntri, [&](int t) {
    ES::M3d M;
    for (int v = 0; v < 3; v++) {
      M.col(v) = ES::toESV3d(sd->sphereMeshCur.pos(t, v)) - sphereCenter;
    }
    sd->triangleRestVols[t] = VegaFEM::NonlinearOptimization::Determinant::Dim3::det(M.data());
    sd->triangleRestAreas[t] = getTriangleArea(sd->sphereMeshCur.pos(t, 0), sd->sphereMeshCur.pos(t, 1), sd->sphereMeshCur.pos(t, 2));
  });
  sd->inversionBarrierEnergy->setRadius(radius);
  sd->inversionBarrierEnergy->setCenter(center);

  // VegaFEM::NonlinearOptimization::FiniteDifference fd(VegaFEM::NonlinearOptimization::FiniteDifference::M_FIVE_POINT, 1e-7);
  // ES::VXd xtest = x;
  // for (int i = 0; i < sphereMeshCur.numVertices(); i++) {
  //   xtest.segment<3>(i * 3) = (x.segment<3>(i * 3) - sphereCenter) * 1e-2 + sphereCenter;
  // }
  // double err[2];
  // fd.testEnergy(barrierEnergy, true, true, -1.0, xtest.data(), -1, err, err + 1);
  // std::cout << err[0] << ',' << err[1] << std::endl;
  // std::cout << barrierEnergy->func(xtest) << std::endl;
  // exit(1);

  sd->energiesXDOFs->setEnergyCoeffs(int(OptEnergyType::ICP_ENERGY), params.icp_w);
  sd->energiesXDOFs->setEnergyCoeffs(int(OptEnergyType::SMOOTHNESS_ENERGY), params.smooth_w_init);
  sd->energiesXDOFs->setEnergyCoeffs(int(OptEnergyType::EXPANSION_ENERGY), params.expansion_w);
  sd->energiesXDOFs->setEnergyCoeffs(int(OptEnergyType::TRANSLATION_ENERGY), params.fixTranslation_w);
  sd->energiesXDOFs->setEnergyCoeffs(int(OptEnergyType::INVERSION_ENERGY), params.barrier_w);

  if (params.numPositionDOFs == NumPositionDOFs::TWELVE) {
    sd->energiesXDOFs->setEnergyCoeffs(int(OptEnergyType::SMOOTHNESS_ENERGY), 0.0);
    sd->energiesAffineDOFs->setEnergyCoeffs(1, params.smooth_w_init);
  }

  if (params.inertial_w > 0) {
    sd->energiesXDOFs->setEnergyCoeffs(int(OptEnergyType::INERTIA_ENERGY), params.inertial_w);
  }

  if (params.sampleCount <= 1) {
    sd->icpEnergyQuadratic->setCheckContact(1);
  }

  sd->lastExpansionTargetPosition.resize(0);
  sd->lastExpansionTargetSamplePosition.resize(0);

  int numDynamicIterations = int(params.numIter * 0.8);
  for (int iter = 1; iter <= params.numIter; iter++) {
    double ratio = (double)(iter - 1) / (params.numWeightIter - 1);
    ratio = std::clamp(ratio, 0.0, 1.0);

    double w_cur = params.smooth_w_init * (1.0 - ratio) + params.smooth_w_final * ratio;
    w_cur = std::pow(10, w_cur);
    // double w_cur = std::exp(std::log(params.smooth_w_init) * (1 - ratio) + std::log(params.smooth_w_final) * ratio);

    fmt::print("Iter: {}\n", iter);
    fmt::print("Alpha: {}\n", w_cur);

    // compute contact constraints and distances
    if (params.sampleCount > 1) {
      tbb::parallel_for(0, (int)sd->sampleRestCoeffs.size(), [&](int si) {
        // for (int si = 0; si < (int)sd->sampleRestCoeffs.size(); si++) {
        Vec3d q(0.0);
        for (int vi = 0; vi < 3; vi++) {
          q += sd->sphereMeshCur.pos(sd->sampleVertexIndices[si * 3 + vi]) * sd->sampleVertexWeights[si * 3 + vi];
        }

        auto ret = targetMeshBVTree.closestTriangleQuery(targetMesh, q);
        Vec3d dir = ret.closestPosition - q;
        Vec3d target_n = targetMeshNormals.getPseudoNormal(targetMesh.triangles().data(), ret.triID, ret.feature);
        double dist = len(dir);
        dir /= dist;

        double dotValue = dot(dir, target_n);

        // in contact
        if (dotValue <= 0 || dist < 1e-3) {
          sd->icpSampleData.coeffs[si] = sd->sampleRestCoeffs[si];
          sd->icpSampleData.targetPositions.segment<3>(si * 3) = ES::toESV3d(ret.closestPosition);
          sd->icpSampleData.targetNormals.segment<3>(si * 3) = -ES::toESV3d(target_n);
        }
        else {
          sd->icpSampleData.coeffs[si] = 0;
        }
      });

      tbb::parallel_for(0, (int)sd->sampleRestCoeffs.size(), [&](int si) {
        // for (int si = 0; si < (int)sd->sampleRestCoeffs.size(); si++) {
        Vec3d q(0.0);
        Vec3d n(0.0);
        for (int vi = 0; vi < 3; vi++) {
          q += sd->sphereMeshCur.pos(sd->sampleVertexIndices[si * 3 + vi]) * sd->sampleVertexWeights[si * 3 + vi];
          n += sd->sphereMeshNormals.vtxNormal(sd->sampleVertexIndices[si * 3 + vi]) * sd->sampleVertexWeights[si * 3 + vi];
        }

        // expansion
        Vec3d ldiff = q - Vec3d(sphereCenter.data());
        double ldiff_len = len(ldiff);

        ES::V3d tgtDir, tgtP;
        tgtDir = ES::toESV3d(n).normalized();
        // ES::toESV3d(ldiff).normalized();

        // does not increase after numDynamicIterations
        // does not increase if there is contact
        if (iter >= numDynamicIterations || sd->icpSampleData.coeffs[si] > 0) {
          if (sd->lastExpansionTargetSamplePosition.size()) {
            tgtP = sd->lastExpansionTargetSamplePosition.segment<3>(si * 3);
          }
          else {
            tgtP = ES::toESV3d(q);
          }
        }
        else {
          double defaultIncLength = std::min(ldiff_len * params.expansionRelStepSize, expansionMaxStepSize);
          double finalLength = defaultIncLength + ldiff_len;
          finalLength = std::min(finalLength, radius * 2.5);

          // sd->expansionData.targetNormals.segment<3>(i * 3) * finalLength + sphereCenter;
          tgtP = std::max(finalLength - ldiff_len, 0.0) * tgtDir + ES::toESV3d(q);
        }

        sd->expansionSampleData.targetNormals.segment<3>(si * 3) = tgtDir;
        sd->expansionSampleData.targetPositions.segment<3>(si * 3) = tgtP;
        sd->expansionSampleData.coeffs[si] = sd->sampleRestCoeffs[si];
      });

      if (dumpMesh) {
        sd->dumpConstraints(true).save(fmt::format("{}/{:03d}exp-pair.obj", prefix, iter));
      }

      sd->icpEnergyQuadraticDenseSampled->setTargetPos(sd->icpSampleData.targetPositions.data());
      sd->icpEnergyQuadraticDenseSampled->setNormals(sd->icpSampleData.targetNormals.data());
      sd->icpEnergyQuadraticDenseSampled->setCoeff(sd->icpSampleData.coeffs.data());

      sd->sphereExpansionEnergyQuadDenseSampled->setTargetPos(sd->expansionSampleData.targetPositions.data());
      sd->sphereExpansionEnergyQuadDenseSampled->setNormals(sd->expansionSampleData.targetNormals.data());
      sd->sphereExpansionEnergyQuadDenseSampled->setCoeff(sd->expansionSampleData.coeffs.data());
      sd->sphereExpansionEnergyQuadDenseSampled->computeHessian();
      sd->sphereExpansionEnergyQuadDenseSampled->printErrorInfo(sd->q);

      sd->lastExpansionTargetSamplePosition = sd->expansionSampleData.targetPositions;
    }
    else {
      tbb::parallel_for(0, nv, [&](int i) {
        // for (int i = 0; i < nv; i++) {
        Vec3d sp_qi = sd->sphereMeshCur.pos(i);
        auto ret = targetMeshBVTree.closestTriangleQuery(targetMesh, sp_qi);
        Vec3d dir = ret.closestPosition - sp_qi;
        Vec3d target_n = targetMeshNormals.getPseudoNormal(targetMesh.triangles().data(), ret.triID, ret.feature);
        double dist = len(dir);
        sd->distAll[i] = dist;

        dir /= dist;
        double dotValue = dot(dir, target_n);

        // in contact
        if (dotValue <= 0 || dist < 1e-3) {
          sd->icpData.coeffs[i] = sd->vtxAreas[i];
        }
        else {
          sd->icpData.coeffs[i] = 0;
        }

        sd->icpData.targetPositions.segment<3>(i * 3) = ES::toESV3d(ret.closestPosition);
        sd->icpData.targetNormals.segment<3>(i * 3) = -ES::toESV3d(target_n);
      });

      tbb::parallel_for(0, nv, [&](int i) {
        // for(int i = 0; i < num_v; i++) {

        const Vec3d &sp_qi = sd->sphereMeshCur.pos(i);
        // expansion
        Vec3d ldiff = sp_qi - Vec3d(sphereCenter.data());
        double ldiff_len = len(ldiff);

        // sd->expansionData.targetNormals.segment<3>(i * 3) = ES::toESV3d(ldiff).normalized();
        ES::V3d tgtDir = ES::toESV3d(sd->sphereMeshNormals.vtxNormal(i));

        ES::V3d tgtP;
        // in both case, we did not change target position
        if (iter >= numDynamicIterations || sd->icpData.coeffs[i] > 0) {
          if (sd->lastExpansionTargetPosition.size()) {
            tgtP = sd->lastExpansionTargetPosition.segment<3>(i * 3);
          }
          else {
            tgtP = ES::toESV3d(sp_qi);
          }
        }
        // otherwise
        else {
          double defaultIncLength = std::min(ldiff_len * params.expansionRelStepSize, expansionMaxStepSize);
          double finalLength = defaultIncLength + ldiff_len;
          finalLength = std::min(finalLength, radius * 2.5);
          tgtP = tgtDir * std::max(finalLength - ldiff_len, 0.0) + ES::toESV3d(sp_qi);
          // sd->expansionData.targetNormals.segment<3>(i * 3) * finalLength + sphereCenter;
        }

        sd->expansionData.targetNormals.segment<3>(i * 3) = tgtDir;
        sd->expansionData.targetPositions.segment<3>(i * 3) = tgtP;
        sd->expansionData.coeffs[i] = sd->vtxAreas[i];

        if (params.expansionConstraintType == ExpansionConstraintType::LINEAR) {
          sd->expansionForce.segment<3>(i * 3) = -sd->expansionData.targetNormals.segment<3>(i * 3) * sd->vtxAreas[i];
        }
      });

      if (dumpMesh) {
        sd->dumpConstraints().save(fmt::format("{}/{:03d}exp-pair.obj", prefix, iter));
        LibIGLInterface::saveOffMesh(fmt::format("{}/{:03d}pair-w.off", prefix, iter).c_str(), sd->sphereMeshCur, sd->icpData.coeffs.data());
      }

      // update energies
      if (sd->icpEnergyQuadratic) {
        sd->icpEnergyQuadratic->setTargetPos(sd->icpData.targetPositions.data());
        sd->icpEnergyQuadratic->setNormals(sd->icpData.targetNormals.data());
        sd->icpEnergyQuadratic->setCoeff(sd->icpData.coeffs.data());
      }

      if (sd->icpEnergyPOrder) {
        sd->icpEnergyPOrder->setTargetPos(sd->icpData.targetPositions.data());
        sd->icpEnergyPOrder->setNormals(sd->icpData.targetNormals.data());
        sd->icpEnergyPOrder->setCoeff(sd->icpData.coeffs.data());

        // VegaFEM::NonlinearOptimization::FiniteDifference fd(VegaFEM::NonlinearOptimization::FiniteDifference::M_FIVE_POINT, 1e-7);
        // double err[2];
        // fd.testEnergy(sd->icpEnergyPOrder, true, true, -1.0, sd->x.data(), -1, err, err + 1);
        // std::cout << err[0] << ',' << err[1] << std::endl;
        // std::cout << barrierEnergy->func(xtest) << std::endl;
        // exit(1);
      }

      // expansion
      if (params.expansionConstraintType == ExpansionConstraintType::QUADRATIC) {
        sd->sphereExpansionEnergyQuad->setTargetPos(sd->expansionData.targetPositions.data());
        sd->sphereExpansionEnergyQuad->setNormals(sd->expansionData.targetNormals.data());
        sd->sphereExpansionEnergyQuad->setCoeff(sd->expansionData.coeffs.data());
        sd->sphereExpansionEnergyQuad->computeHessian();

        sd->lastExpansionTargetPosition = sd->expansionData.targetPositions;
        sd->sphereExpansionEnergyQuad->printErrorInfo(sd->q);
      }
      else if (params.expansionConstraintType == ExpansionConstraintType::LINEAR) {
        // linear one
        double curVol = 0;
        for (int ti = 0; ti < ntri; ti++) {
          ES::M3d M;
          M.col(0) = ES::toESV3d(sd->sphereMeshCur.pos(ti, 0)) - sphereCenter;
          M.col(1) = ES::toESV3d(sd->sphereMeshCur.pos(ti, 1)) - sphereCenter;
          M.col(2) = ES::toESV3d(sd->sphereMeshCur.pos(ti, 2)) - sphereCenter;
          curVol += std::abs(VegaFEM::NonlinearOptimization::Determinant::Dim3::det(M.data()));
        }
        curVol /= 6.0;
        sd->expansionForce *= std::min(initialVol / curVol, 1.0);
      }
    }

    if (params.addConstraints) {
      double lb = std::max(initialVol / sd->sphereMeshRest.numTriangles() * 1e-2, 1e-10);

      sd->clow.setConstant(lb);
      sd->constraints->setCenter(center);
      sd->icpSolver.solver->setcRange(sd->clow.data(), sd->chi.data());
    }

    // update inertial
    if (params.inertial_w > 0) {
      // update A_inertia
      // M acc + D vel
      // M 1/h(v-v0) + D 1/h (x - x0)
      // M 1/h(1/h(x - x0) -v0) + aM 1/h(x - x0)
      // M (1/h2 x - 1/h2 x0 -1/h v0)+ aM 1/h(x - x0)
      // (1/h2 M + a/h M) x - M (1/h2 x0 + 1/h v0 + a/h x0)
      double s = 1.0 / (params.timestep * params.timestep) + dm / params.timestep;
      sd->temp.noalias() = sd->curPosition * s;
      sd->temp.noalias() += sd->curVel / params.timestep;
      sd->temp *= -1.0;

      ES::Mv(sd->M3, sd->temp, sd->b_inertial);

      // dynamics weights
      // if (iter > numDynamicIterations) {
      //   sd->energiesXDOFs->setEnergyCoeffs(int(OptEnergyType::INERTIA_ENERGY), 0);
      // }
    }

    // smoothness weight
    if (params.numPositionDOFs == NumPositionDOFs::THREE) {
      sd->energiesXDOFs->setEnergyCoeffs(int(OptEnergyType::SMOOTHNESS_ENERGY), w_cur);
    }
    else if (params.numPositionDOFs == NumPositionDOFs::TWELVE) {
      sd->energiesXDOFs->setEnergyCoeffs(int(OptEnergyType::SMOOTHNESS_ENERGY), 0);
      sd->energiesAffineDOFs->setEnergyCoeffs(1, w_cur);
    }

    // sd->energies->setEnergyCoeffs(int(OptEnergyType::INERTIA_ENERGY), 0);
    // energies->setEnergyCoeffs(int(EnergyType::EXPANSION_ENERGY), 0);

    if (params.centerFlag == CenterOptimizationFlag::OPTIM) {
      if (params.numPositionDOFs == NumPositionDOFs::THREE) {
        sd->q.segment<3>(nv3) = sphereCenter;
        // sd->xlow.segment<3>(nv3) = sphereCenter;
        // sd->xhi.segment<3>(nv3) = sphereCenter;

        ES::Mp<ES::M3d>(sd->q.data() + nv3 + 3).setIdentity();
        ES::Mp<ES::M3d>(sd->qlow.data() + nv3 + 3).setIdentity();
        ES::Mp<ES::M3d>(sd->qhi.data() + nv3 + 3).setIdentity();
      }
      else if (params.numPositionDOFs == NumPositionDOFs::TWELVE) {
        sd->q.segment<3>(nv * 12) = sphereCenter;

        ES::Mp<ES::M3d>(sd->q.data() + nv * 12 + 3).setIdentity();
        ES::Mp<ES::M3d>(sd->qlow.data() + nv * 12 + 3).setIdentity();
        ES::Mp<ES::M3d>(sd->qhi.data() + nv * 12 + 3).setIdentity();
      }
    }
    else if (params.centerFlag == CenterOptimizationFlag::FOLLOW) {
      // update translation
      sd->sphereTranslationEnergy->setRestCenter();
      restCenter = sd->sphereTranslationEnergy->getRestCenter();

      // update barrier
      tbb::parallel_for(0, ntri, [&](int t) {
        ES::M3d M;
        for (int v = 0; v < 3; v++) {
          M.col(v) = ES::toESV3d(sd->sphereMeshCur.pos(t, v)) - restCenter;
        }
        sd->triangleRestVols[t] = VegaFEM::NonlinearOptimization::Determinant::Dim3::det(M.data());
        sd->triangleRestAreas[t] = getTriangleArea(sd->sphereMeshCur.pos(t, 0), sd->sphereMeshCur.pos(t, 1), sd->sphereMeshCur.pos(t, 2));
      });

      double radiusCur = 0.0;
      for (int vi = 0; vi < nv; vi++) {
        radiusCur += (ES::toESV3d(sd->sphereMeshCur.pos(vi)) - restCenter).norm();
      }
      radiusCur /= nv;

      sd->inversionBarrierEnergy->setRadius(radiusCur);
      sd->inversionBarrierEnergy->setCenter(restCenter.data());

      // update center
      sd->sphereTranslationEnergy->setCenter(restCenter.data());
      sphereCenter = restCenter;
    }
    else if (params.centerFlag == CenterOptimizationFlag::FIXED) {
      // update center
      sd->sphereTranslationEnergy->setCenter(restCenter.data());
    }

    std::cout << "Before Opt:\n";
    if (params.numPositionDOFs == NumPositionDOFs::THREE) {
      sd->energiesXDOFs->printEnergy(sd->q);
    }
    else if (params.numPositionDOFs == NumPositionDOFs::TWELVE) {
      sd->energiesXDOFs->printEnergy(sd->curPosition);
      sd->energiesAffineDOFs->printEnergy(sd->q);

      // sd->energiesXDOFs->setEnergyCoeffs(int(OptEnergyType::SMOOTHNESS_ENERGY), 0.0);

      // VegaFEM::NonlinearOptimization::FiniteDifference fd(VegaFEM::NonlinearOptimization::FiniteDifference::M_FIVE_POINT, 1e-7);
      // double err[2];
      // fd.testEnergy(sd->affineSpaceEnergies, true, true, -1.0, sd->q.data(), -1, err, err + 1);
      // std::cout << err[0] << ',' << err[1] << std::endl;
      // exit(1);
    }

#if defined(USE_KNITRO)
    sd->curPosition1.noalias() = sd->curPosition;

    if (params.centerFlag == CenterOptimizationFlag::OPTIM) {
      sd->icpSolver.solver->setxRange(sd->qlow.data(), sd->qhi.data());
    }

    sd->icpSolver.solver->setxInit(sd->q.data());
    int solverRet = sd->icpSolver.solver->solve();
    sd->icpSolver.getx(sd->q);

    if (solverRet != 0) {
      std::cerr << "Solver non zero ret: " << solverRet << std::endl;
      /*meshOut = sd->sphereMeshCur;

      if (dumpMesh) {
        sd->sphereMeshCur.save(fmt::format("{}-icp-mesh{:02d}.obj", prefix, iter));
      }
      return 1;*/
    }
#endif

    std::cout << "After Opt:\n";
    if (params.numPositionDOFs == NumPositionDOFs::THREE) {
      sd->energiesXDOFs->printEnergy(sd->q);

      if (sd->sphereExpansionEnergyQuadDenseSampled)
        sd->sphereExpansionEnergyQuadDenseSampled->printErrorInfo(sd->q);
      else if (sd->sphereExpansionEnergyQuad)
        sd->sphereExpansionEnergyQuad->printErrorInfo(sd->q);

      sd->curPosition.noalias() = sd->q.head(nv3);
    }
    else if (params.numPositionDOFs == NumPositionDOFs::TWELVE) {
      tbb::parallel_for(
        0, nv, [&](int vi) {
          Eigen::Matrix<double, 3, 4> A = ES::Mp<Eigen::Matrix<double, 3, 4>>(sd->q.data() + vi * 12);
          ES::V3d x0 = sd->restPosition.segment<3>(vi * 3);

          sd->curPosition.segment<3>(vi * 3) = A.block<3, 3>(0, 0) * x0 + A.col(3);
        },
        tbb::static_partitioner());

      sd->energiesXDOFs->printEnergy(sd->curPosition);
      sd->energiesAffineDOFs->printEnergy(sd->q);
    }

    sd->curVel.noalias() = (sd->curPosition - sd->curPosition1) / params.timestep;

    if (params.centerFlag == CenterOptimizationFlag::OPTIM) {
      std::cout << "Old center: " << sphereCenter.transpose() << std::endl;
      if (params.numPositionDOFs == NumPositionDOFs::THREE) {
        std::cout << "New center: " << sd->q.segment<3>(nv3).transpose() << std::endl;
      }
      else if (params.numPositionDOFs == NumPositionDOFs::TWELVE) {
        std::cout << "New center: " << sd->q.segment<3>(nv * 12).transpose() << std::endl;
      }
    }

    for (int i = 0; i < sd->sphereMeshCur.numVertices(); i++) {
      sd->sphereMeshCur.pos(i) = Vec3d(sd->curPosition.data() + i * 3);
    }

    sd->sphereMeshNormals.updateVertexPositions(sd->sphereMeshCur);

    // update smoothness
    // 1/2 (x - x0)^T biL (x - x0)
    // 1/2 xT biL x - x0^T biL x + x0^T BiL x0
    // TriMeshGeo aa = sd->sphereMeshRest;
    // sd->sphereMeshRest = sd->sphereMeshCur;
    // sd->updateLaplacianAndMass();

    // sd->sphereMeshRest = aa;

    // ES::Mv(sd->biL3, sd->curPosition, sd->b_smooth);
    // sd->b_smooth *= -1.0;

    if (dumpMesh) {
      sd->sphereMeshCur.save(fmt::format("{}/{:03d}icp-mesh.obj", prefix, iter));
    }

    fmt::print("Iter {} done.\n", iter);
    // exit(1);
    // aaa
  }

  fmt::print("Final sphere center: {}, {}, {}\n", sphereCenter[0], sphereCenter[1], sphereCenter[2]);
  std::cout << std::flush;

  if (updatedCenter) {
    *updatedCenter = Vec3d(sphereCenter.data());
  }

  // sphereCorefine(sd->sphereMeshCur, sphereCenter.data(), targetMesh, targetMeshBVTree, targetMeshNormals, meshOut);

  meshOut = sd->sphereMeshCur;

  return 0;
}

#if 0
int MedialAxisRepresentation::sphereICP(const TriMeshGeo &sphereMesh,
  const TriMeshGeo &targetMesh, const TriMeshBVTree &targetMeshBVTree, const TriMeshPseudoNormal &targetMeshNormals,
  double w_init, double w_final, int numIter, double timestep, const std::string &prefix, TriMeshGeo &meshOut)
{
  int dumpMesh = 0;
  if (prefix.length()) {
    dumpMesh = 1;
  }

  // compute center
  ES::V3d sphereCenter(0.0, 0.0, 0.0);
  for (int i = 0; i < sphereMesh.numVertices(); i++) {
    sphereCenter += ES::toESV3d(sphereMesh.pos(i));
  }
  sphereCenter /= sphereMesh.numVertices();

  // compute sphere radius
  double radius = 0.0;
  for (int i = 0; i < sphereMesh.numVertices(); i++) {
    radius += (ES::toESV3d(sphereMesh.pos(i)) - sphereCenter).norm();
  }
  radius /= sphereMesh.numVertices();

  // initial sphere volume
  double initialVol = 4.0 / 3.0 * M_PI * radius * radius * radius;
  double initialArea = 4.0 * M_PI * radius * radius;

  std::vector<int> sphereSamples;
  while (sphereSamples.size() < 100) {
    int furthestOne = -1;
    double maxDist = 0;
    for (int i = 0; i < sphereMesh.numVertices(); i++) {
      double closestSampleDist = 1e100;
      for (int vid : sphereSamples) {
        Vec3d diff = sphereMesh.pos(i) - sphereMesh.pos(vid);
        closestSampleDist = std::min(closestSampleDist, len2(diff));
      }
      if (closestSampleDist > maxDist) {
        maxDist = closestSampleDist;
        furthestOne = i;
      }
    }

    sphereSamples.emplace_back(furthestOne);
  }
  TriMeshGeo kk;
  for (int vid : sphereSamples) {
    kk.addPos(sphereMesh.pos(vid));
  }
  kk.save("zz.obj");

  // get matrix from sphere
  ES::MXd sphereRestPostions;
  ES::MXi sphereTriangles;
  BasicUtilities::triMeshGeoToMatrices(sphereMesh, sphereRestPostions, sphereTriangles);

  // bi laplacian matrix
  ES::SpMatD biL;
  biLaplacian(biL, sphereRestPostions, sphereTriangles);

  // mass matrix
  ES::SpMatD M;
  massMatrix(M, sphereRestPostions, sphereTriangles);

  // normals
  TriMeshGeo sphereMeshCur = sphereMesh;
  TriMeshPseudoNormal sphereMeshNormal;
  sphereMeshNormal.buildPseudoNormals(sphereMeshCur);

  // optimization variables
  ES::VXd xlow, xhi;
  ES::VXd g, lambda, clow, chi;
  xlow.setConstant(sphereMesh.numVertices() * 3, -1e20);
  xhi.setConstant(sphereMesh.numVertices() * 3, 1e20);
  g.setZero(sphereMesh.numTriangles());
  lambda.setZero(sphereMesh.numTriangles());

  double lb = std::max(initialVol / sphereMesh.numTriangles() * 1e-3, 1e-10);
  clow.setConstant(sphereMesh.numTriangles(), lb);
  chi.setConstant(sphereMesh.numTriangles(), 1e20);

  // x
  ES::VXd x(sphereMesh.numVertices() * 3), xvel(sphereMesh.numVertices() * 3), x1(sphereMesh.numVertices() * 3);
  for (int vi = 0; vi < sphereMesh.numVertices(); vi++) {
    x.segment<3>(vi * 3) = sphereRestPostions.row(vi);
    xvel.segment<3>(vi * 3).setZero();
    x1.segment<3>(vi * 3) = x.segment<3>(vi * 3);
  }
  ES::VXd xrest = x;

  struct ICPData
  {
    ES::VXd coeffs;
    ES::VXd targetPositions;
    ES::VXd targetNormals;
    std::vector<int> vertexIndices;

    ICPData(int n):
      coeffs(n), targetPositions(n * 3), targetNormals(n * 3), vertexIndices(n)
    {
      std::iota(vertexIndices.begin(), vertexIndices.end(), 0);
    }
  };

  ICPData icpData(sphereMesh.numVertices()), expansionData(sphereMesh.numVertices());
  ES::VXd vtxAreas(sphereMesh.numVertices());
  ES::VXd distAll(sphereMeshCur.numVertices());
  ES::VXd temp = x;
  std::vector<tbb::spin_mutex> lock(sphereMesh.numVertices());
  std::unordered_set<int> invalidConstraintsSet;

  // vtx area
  vtxAreas.setZero();
  tbb::parallel_for(0, sphereMeshCur.numTriangles(), [&](int t) {
    double area_t = getTriangleArea(sphereMeshCur.pos(t, 0), sphereMeshCur.pos(t, 1), sphereMeshCur.pos(t, 2));

    for (int ti = 0; ti < 3; ti++) {
      int v_ti = sphereMeshCur.triVtxID(t, ti);

      lock.at(v_ti).lock();
      vtxAreas(v_ti) += area_t / 3.0;
      lock.at(v_ti).unlock();
    }
  });
  vtxAreas /= vtxAreas.maxCoeff();

  ES::SpMatD baseHessian(sphereMesh.numVertices() * 3, sphereMesh.numVertices() * 3);
  std::vector<ES::TripletD> baseHessianEntries;
  baseHessianEntries.reserve(81 * sphereMesh.numTriangles());
  for (int t = 0; t < sphereMesh.numTriangles(); t++) {
    for (int vi = 0; vi < 3; vi++) {
      for (int vj = 0; vj < 3; vj++) {
        int vid_i = sphereMesh.triVtxID(t, vi);
        int vid_j = sphereMesh.triVtxID(t, vj);
        for (int i = 0; i < 3; i++) {
          for (int j = 0; j < 3; j++) {
            baseHessianEntries.emplace_back(vid_i * 3 + i, vid_j * 3 + j, 0.0);
          }
        }
      }
    }
  }
  baseHessian.setFromTriplets(baseHessianEntries.begin(), baseHessianEntries.end());

  // inertial
  double inertiaWeights = 0.0;
  if (timestep > 0) {
    inertiaWeights = 1.0;
  }
  else {
    timestep = 1.0;
    inertiaWeights = 0.0;
  }

  // inertia
  // M 1/h(1/h(x - x0) -v0)
  // M (1/h2 x - 1/h2 x0 -1/h v0)
  ES::SpMatD A_inertial = M / (timestep * timestep);
  ES::VXd b_inertial(sphereMesh.numVertices() * 3);
  std::shared_ptr<VegaFEM::NonlinearOptimization::QuadraticPotentialEnergy> inertiaEnergy = std::make_shared<VegaFEM::NonlinearOptimization::QuadraticPotentialEnergy>(A_inertial, b_inertial);

  // icp energy
  std::shared_ptr<VegaFEM::PredefinedPotentialEnergies::MultipleVertexSliding> icpEnergy =
    std::make_shared<VegaFEM::PredefinedPotentialEnergies::MultipleVertexSliding>(
      baseHessian, xrest, sphereMesh.numVertices(), icpData.coeffs.data(), icpData.targetPositions.data(), icpData.targetNormals.data(),
      icpData.vertexIndices.data(), 1.0, 1, 0, 0);

  // smoothness energy
  // 1/2 (x - x0)^T biL (x - x0)
  // 1/2 xT biL x - x0^T biL x + x0^T BiL x0
  ES::VXd b_smooth(sphereMesh.numVertices() * 3);
  ES::Mv(biL, xrest, b_smooth);
  b_smooth *= -1.0;
  std::shared_ptr<VegaFEM::NonlinearOptimization::QuadraticPotentialEnergy> smoothnessEnergy = std::make_shared<VegaFEM::NonlinearOptimization::QuadraticPotentialEnergy>(biL, b_smooth);

  // expansion
  std::shared_ptr<VegaFEM::PredefinedPotentialEnergies::MultipleVertexSliding> sphereExpansionEnergyQuad =
    std::make_shared<VegaFEM::PredefinedPotentialEnergies::MultipleVertexSliding>(
      baseHessian, xrest, sphereMesh.numVertices(), expansionData.coeffs.data(), expansionData.targetPositions.data(), expansionData.targetNormals.data(),
      expansionData.vertexIndices.data(), 1.0, 1, 0, 0);

  ES::VXd expansionForce(sphereMesh.numVertices() * 3);
  std::shared_ptr<VegaFEM::NonlinearOptimization::LinearPotentialEnergy> sphereExpansionEnergyLinear = std::make_shared<VegaFEM::NonlinearOptimization::LinearPotentialEnergy>(expansionForce);

  // translation
  std::shared_ptr<VegaFEM::PredefinedPotentialEnergies::MultipleVertexConstrainedRigidMotion> sphereTranslationEnergy =
    std::make_shared<VegaFEM::PredefinedPotentialEnergies::MultipleVertexConstrainedRigidMotion>(
      xrest, sphereSamples, nullptr, 0.0, 1.0, false, false);

  // barrier
  ES::VXd triangleRestVols(sphereMesh.numTriangles());
  ES::VXd triangleRestAreas(sphereMesh.numTriangles());

  tbb::parallel_for(0, sphereMesh.numTriangles(), [&](int t) {
    ES::M3d M;
    for (int v = 0; v < 3; v++) {
      M.col(v) = ES::toESV3d(sphereMeshCur.pos(t, v)) - sphereCenter;
    }
    triangleRestVols[t] = VegaFEM::NonlinearOptimization::Determinant::Dim3::det(M.data());
    triangleRestAreas[t] = getTriangleArea(sphereMeshCur.pos(t, 0), sphereMeshCur.pos(t, 1), sphereMeshCur.pos(t, 2));
  });

  std::shared_ptr<SphereInversionBarrierEnergy> inversionBarrierEnergy = std::make_shared<SphereInversionBarrierEnergy>(
    sphereCenter, sphereMeshCur, triangleRestVols, triangleRestAreas, radius, false, 1e-1, 1e-1, 1e-4, 4);

  // VegaFEM::NonlinearOptimization::FiniteDifference fd(VegaFEM::NonlinearOptimization::FiniteDifference::M_FIVE_POINT, 1e-7);
  // ES::VXd xtest = x;
  // for (int i = 0; i < sphereMeshCur.numVertices(); i++) {
  //   xtest.segment<3>(i * 3) = (x.segment<3>(i * 3) - sphereCenter) * 1e-2 + sphereCenter;
  // }
  // double err[2];
  // fd.testEnergy(barrierEnergy, true, true, -1.0, xtest.data(), -1, err, err + 1);
  // std::cout << err[0] << ',' << err[1] << std::endl;
  // std::cout << barrierEnergy->func(xtest) << std::endl;
  // exit(1);

  std::shared_ptr<VegaFEM::NonlinearOptimization::PotentialEnergies> energies = std::make_shared<VegaFEM::NonlinearOptimization::PotentialEnergies>(sphereMesh.numVertices() * 3);
  energies->addPotentialEnergy(icpEnergy, 1e2);
  energies->addPotentialEnergy(smoothnessEnergy, w_init);
  energies->addPotentialEnergy(sphereExpansionEnergyQuad, 1);
  energies->addPotentialEnergy(inertiaEnergy, inertiaWeights);
  energies->addPotentialEnergy(sphereTranslationEnergy, 1e7);
  energies->addPotentialEnergy(inversionBarrierEnergy, 1);

  energies->init();

  std::shared_ptr<SphereTriangleFlippingConstraints> constraints = std::make_shared<SphereTriangleFlippingConstraints>(sphereCenter, sphereTriangles, sphereMesh.numVertices());

#  if defined(USE_KNITRO)
  struct ICPSolver
  {
    ICPSolver(VegaFEM::NonlinearOptimization::PotentialEnergy_const_p energy, VegaFEM::NonlinearOptimization::ConstraintFunctions_const_p constraints,
      ES::ConstRefVecXd x, ES::ConstRefVecXd xlow, ES::ConstRefVecXd xhi, ES::ConstRefVecXd clow, ES::ConstRefVecXd chi)
    {
      if (constraints) {
        problem = std::make_unique<VegaFEM::NonlinearOptimization::KnitroProblem>(energy, constraints);
      }
      else {
        problem = std::make_unique<VegaFEM::NonlinearOptimization::KnitroProblem>(energy);
      }

      problem->setInit(x);
      problem->setRange(xlow, xhi);

      if (constraints) {
        problem->setConstraintsRange(clow, chi);
      }

      solver = std::make_unique<VegaFEM::NonlinearOptimization::KnitroOptimizer>(problem.get());
      solver->setConfigFile("config.opt");

      solver->setMaxIter(300);

      solver->setOptTol(1e-6);
      solver->setFeasTol(1e-8);
      solver->setVerbose(0);

      solver->enableMultiEvaluation(0);
      solver->init();
    }

    void getx(ES::VXd &x) const
    {
      x.noalias() = Eigen::Map<const ES::VXd>(solver->getx(), problem->getn());
    }

    std::unique_ptr<VegaFEM::NonlinearOptimization::KnitroProblem> problem;
    std::unique_ptr<VegaFEM::NonlinearOptimization::KnitroOptimizer> solver;
  };

  ICPSolver freeICPSolver(energies, nullptr, x, xlow, xhi, clow, chi);
  ICPSolver constrainedICPSolver(energies, constraints, x, xlow, xhi, clow, chi);

#  endif

  auto dumpConstraints = [&]() -> TriMeshGeo {
    // dump constraints
    TriMeshGeo abc;
    for (int i = 0; i < sphereMeshCur.numVertices(); i++) {
      Vec3d p0 = sphereMeshCur.pos(i);

      if (expansionData.coeffs[i] > 0) {
        Vec3d p1(expansionData.targetPositions.data() + i * 3);
        abc.addMesh(createSingleTriangleMesh(p0, p1, p0 + Vec3d(1e-5)));
      }
      else {
        Vec3d p1(icpData.targetPositions.data() + i * 3);
        abc.addMesh(createSingleTriangleMesh(p0, p1, p0 + Vec3d(1e-5)));
      }
    }
    return abc;
  };

  int numDynamicIterations = int(numIter * 0.8);
  for (int iter = 1; iter <= numIter; iter++) {
    double ratio = (double)(iter - 1) / (numDynamicIterations - 1);
    ratio = std::clamp(ratio, 0.0, 1.0);

    double w_cur = w_init * (1.0 - ratio) + w_final * ratio;
    // double w_cur = std::exp(std::log(w_init) * (1 - ratio) + std::log(w_final) * ratio);

    fmt::print("Iter: {}\n", iter);
    fmt::print("Alpha: {}\n", w_cur);

    // compute icp constraints and distances
    tbb::parallel_for(0, sphereMeshCur.numVertices(), [&](int i) {
      // for(int i = 0; i < num_v; i++) {
      Vec3d sp_qi = sphereMeshCur.pos(i);
      auto ret = targetMeshBVTree.closestTriangleQuery(targetMesh, sp_qi);
      Vec3d dir = ret.closestPosition - sp_qi;
      Vec3d target_n = targetMeshNormals.getPseudoNormal(targetMesh.triangles().data(), ret.triID, ret.feature);
      Vec3d source_n = sphereMeshNormal.vtxNormal(i);
      double dist = len(dir);
      distAll[i] = dist;

      if (dot(target_n, source_n) > std::cos(60.0 / 180.0 * M_PI)) {
        // || dot(target_n, trans) < 0.0) {  // && dist <= 0.8 * r) { // smaller than 60 degrees, cos(60), or outside of the mesh
        icpData.coeffs(i) = vtxAreas[i];
        expansionData.coeffs(i) = 0;
        distAll[i] = dist;
      }
      else {
        icpData.coeffs(i) = 0.0;
        expansionData.coeffs(i) = vtxAreas[i];
        distAll[i] = 0;
      }

      icpData.targetPositions.segment<3>(i * 3) = ES::toESV3d(ret.closestPosition);
      icpData.targetNormals.segment<3>(i * 3) = ES::toESV3d(target_n);
    });

    // compute expansion constraints
    double dMax = distAll.maxCoeff();
    double dMid = 0.0;
    int dcount = 0;

    for (int i = 0; i < (int)distAll.size(); i++) {
      if (distAll[i] > 0) {
        dMid += distAll[i];
        dcount += 1;
      }
    }
    dMid /= dcount;

    tbb::parallel_for(0, sphereMeshCur.numVertices(), [&](int i) {
      // for(int i = 0; i < num_v; i++) {
      Vec3d sp_qi = sphereMeshCur.pos(i);
      // expansion
      Vec3d ldiff = sp_qi - Vec3d(sphereCenter.data());
      double ldiff_len = len(ldiff);
      expansionData.targetNormals.segment<3>(i * 3) = ES::toESV3d(ldiff / ldiff_len);

      double defaultIncLength = std::min(ldiff_len * 0.1, 0.02);
      double icpIncLength = dMid;
      double incLength = std::max(icpIncLength, defaultIncLength);  //      *(1 - ratio);
      // double inc_len = std::max(dMax, std::min(ldiff_len * 0.1, 0.02)) * (1 - ratio);

      // double finalLength = std::max(ldiff_len + inc_len, radius);
      double finalLength = incLength + ldiff_len;

      expansionData.targetPositions.segment<3>(i * 3) = expansionData.targetNormals.segment<3>(i * 3) * finalLength + sphereCenter;

      expansionForce.segment<3>(i * 3) = -expansionData.targetNormals.segment<3>(i * 3) * vtxAreas[i];
    });

    if (0) {
      // filtering constraints
      x1.noalias() = x;
      for (int i = 0; i < sphereMeshCur.numVertices(); i++) {
        if (icpData.coeffs[i] > 0)
          x1.segment<3>(i * 3) = icpData.targetPositions.segment<3>(i * 3);
      }
      constraints->func(x1, g);

      invalidConstraintsSet.clear();
      for (int i = 0; i < (int)g.size(); i++) {
        if (g[i] < lb) {
          for (int j = 0; j < 3; j++) {
            int vid = sphereMeshCur.triVtxID(i, j);
            icpData.coeffs[vid] = vtxAreas[vid] * 0.0;
            expansionData.coeffs[vid] = vtxAreas[vid];
            invalidConstraintsSet.emplace(vid);
          }
        }
      }
      fmt::print("drop: {:4d}/{:d}\n", invalidConstraintsSet.size(), sphereMeshCur.numVertices());
    }

    if (dumpMesh) {
      dumpConstraints().save(fmt::format("{}-exp-pair{:02d}.obj", prefix, iter));
    }

    // update energies
    icpEnergy->setTargetPos(icpData.targetPositions.data());
    icpEnergy->setNormals(icpData.targetNormals.data());
    icpEnergy->setCoeff(icpData.coeffs.data());
    icpEnergy->computeHessian();

    // expansion
    sphereExpansionEnergyQuad->setTargetPos(expansionData.targetPositions.data());
    sphereExpansionEnergyQuad->setNormals(expansionData.targetNormals.data());
    sphereExpansionEnergyQuad->setCoeff(expansionData.coeffs.data());
    sphereExpansionEnergyQuad->computeHessian();

    // linear one
    double curVol = 0;
    for (int ti = 0; ti < sphereMeshCur.numTriangles(); ti++) {
      ES::M3d M;
      M.col(0) = ES::toESV3d(sphereMeshCur.pos(ti, 0)) - sphereCenter;
      M.col(1) = ES::toESV3d(sphereMeshCur.pos(ti, 1)) - sphereCenter;
      M.col(2) = ES::toESV3d(sphereMeshCur.pos(ti, 2)) - sphereCenter;
      curVol += std::abs(VegaFEM::NonlinearOptimization::Determinant::Dim3::det(M.data()));
    }
    curVol /= 6.0;

    expansionForce *= std::min(initialVol / curVol, 1.0);

    // update center
    sphereTranslationEnergy->setCenter(sphereCenter.data());

    // update inertial
    // M (1/h2 x - 1/h2 x0 -1/h v0)
    ES::VXd temp = x / (timestep * timestep) + xvel / timestep;
    temp *= -1.0;
    ES::Mv(M, temp, b_inertial);

    // smoothness weight
    energies->setEnergyCoeffs(int(OptEnergyType::SMOOTHNESS_ENERGY), w_cur);
    // energies->setEnergyCoeffs(2, alpha * 10);
    // energies->setEnergyCoeffs(2, alpha * 10 * std::exp(-5 * r));
    // energies->setEnergyCoeffs(2, initialVol / curVol * 0.5 * (1 - ratio));
    // energies->setEnergyCoeffs(2, alpha / ((r / min_bbx_edge) * (r / min_bbx_edge)));

    // dynamics weights
    if (iter > numDynamicIterations) {
      energies->setEnergyCoeffs(int(OptEnergyType::INERTIA_ENERGY), 0);
    }
    energies->setEnergyCoeffs(int(OptEnergyType::INERTIA_ENERGY), 0);
    // energies->setEnergyCoeffs(int(EnergyType::EXPANSION_ENERGY), 0);

    energies->printEnergy(x);
#  if defined(USE_KNITRO)
    x1.noalias() = x;

    freeICPSolver.solver->setxInit(x.data());
    int solverRet = freeICPSolver.solver->solve();
    freeICPSolver.getx(x);

    if (solverRet != 0) {
      std::cerr << "Solver non zero ret: " << solverRet << std::endl;
      meshOut = sphereMeshCur;

      return 1;
    }
    else {
      std::cout << "Solved.\n";
    }

    energies->printEnergy(x);
    xvel.noalias() = (x - x1) / timestep;
#  endif

    for (int i = 0; i < sphereMeshCur.numVertices(); i++) {
      sphereMeshCur.pos(i) = Vec3d(x(i * 3), x(i * 3 + 1), x(i * 3 + 2));
    }

    sphereMeshNormal.updateVertexPositions(sphereMeshCur);

    if (dumpMesh) {
      sphereMeshCur.save(fmt::format("{}-icp-mesh{:02d}.obj", prefix, iter));
    }

    fmt::print("Iter {} done.\n", iter);
    // aaa
  }

  meshOut = sphereMeshCur;

  return 0;
}
#endif