#include "medialSkeletalDiagramOpt.h"
#include "skeletonRVD.h"

#include "primitives/primitiveICP.h"

#include "pgoLogging.h"
#include "EigenSupport.h"
#include "basicAlgorithms.h"
#include "geometryQuery.h"
#include "tetMeshGeo.h"
#include "tetMeshManifold.h"
#include "tetgenInterface.h"
#include "configFileJSON.h"
#include "libiglInterface.h"

#include <fmt/format.h>
#include <CGAL/AABB_triangle_primitive.h>
#include <CGAL/AABB_traits.h>

#include <nlopt.h>

#include <tbb/parallel_for.h>

#include <queue>
#include <filesystem>

namespace MedialAxisRepresentation
{

namespace ES = pgo::EigenSupport;
using CGALPrimitive = CGAL::AABB_triangle_primitive<EK, EK::Triangle_3 *>;

inline ES::V3d toVec3(const EK::Point_3 &pt)
{
  return ES::V3d(CGAL::to_double(pt.x()), CGAL::to_double(pt.y()), CGAL::to_double(pt.z()));
}

struct IterationData
{
  IterationData(const pgo::Mesh::TriMeshGeo &targetMesh_, const pgo::Mesh::TriMeshBVTree &targetMeshBVTree_, const pgo::Mesh::TriMeshPseudoNormal &targetNormals_,
    const pgo::Mesh::TriMeshGeo &targetMeshSmall_, const pgo::Mesh::TriMeshBVTree &targetMeshBVTreeSmall_, const pgo::Mesh::TriMeshPseudoNormal &targetNormalsSmall_):
    targetMesh(targetMesh_),
    targetMeshBVTree(targetMeshBVTree_), targetNormals(targetNormals_),
    targetMeshSmall(targetMeshSmall_), targetMeshBVTreeSmall(targetMeshBVTreeSmall_), targetNormalsSmall(targetNormalsSmall_) {}

  const pgo::Mesh::TriMeshGeo &targetMesh;
  const pgo::Mesh::TriMeshBVTree &targetMeshBVTree;
  const pgo::Mesh::TriMeshPseudoNormal &targetNormals;
  std::vector<double> targetVertexWeights;

  pgo::Mesh::TetMeshGeo targetMeshTetMesh;
  std::vector<std::vector<int>> targetMeshTetMeshTetNeighbors;
  std::vector<ES::V3d> targetMeshTetMeshTetCenters;

  const pgo::Mesh::TriMeshGeo &targetMeshSmall;
  const pgo::Mesh::TriMeshBVTree &targetMeshBVTreeSmall;
  const pgo::Mesh::TriMeshPseudoNormal &targetNormalsSmall;

  double avgLength = 0;
  int iter = 0;
  int computeUnmaskedRegion = 0;
  int maxIter = 20;
  int nAddedPts = 4;

  std::vector<ES::V3d> samplePoints, samplePoints0;
  std::vector<EK::Point_3> medialAxisVerticesInitial;
  std::vector<std::vector<int>> medialAxisFacetsInitial;
  std::vector<ES::V3d> medialAxisVerticesInitialIK;

  std::vector<EK::Point_3> medialAxisVerticesCur;
  std::vector<std::vector<int>> medialAxisFacetsCur;

  std::vector<ES::V3d> finalSkeletonPoints;
  std::vector<std::pair<int, int>> finalSkeletonEdges;
  std::vector<std::tuple<int, int, int>> finalSkeletonTriangles;
  std::vector<pgo::Mesh::TriMeshGeo> finalFitMeshes;
  std::vector<ES::V3d> cellCenters;

  std::vector<ES::V3d> newCenters;
  std::vector<ES::V3d> oldCenters;

  std::vector<int> finalSelectedMAVertices;
};

double funcMinimizeCoverage(unsigned n, const double *x, double *grad, void *my_func_data)
{
  constexpr int dumpMesh = 0;

  IterationData *itData = reinterpret_cast<IterationData *>(my_func_data);
  if (itData == nullptr) {
    abort();
  }
  int maxIter = itData->maxIter;
  int sc = 1e-1;

  // this is the good one
  // int sc = 1e-1
  // int nExpIter = 10;

  PGO_ALOG(n % 3 == 0);

  if constexpr (0) {
    itData->samplePoints.assign(n / 3, ES::V3d::Zero());

    pgo::Mesh::VertexBVTree bvTree;
    bvTree.buildByInertiaPartition(itData->medialAxisVerticesInitialIK);

    itData->samplePoints.assign(n / 3 + itData->oldCenters.size(), ES::V3d::Zero());

    for (int i = 0; i < n / 3; i++) {
      ES::V3d pp(x[i * 3], x[i * 3 + 1], x[i * 3 + 2]);
      int vid = bvTree.getClosestVertex(itData->medialAxisVerticesInitialIK, pp);

      // itData->samplePoints[i] = ES::V3d(x[i * 3], x[i * 3 + 1], x[i * 3 + 2]);
      itData->samplePoints[i] = itData->medialAxisVerticesInitialIK[vid];
    }
  }
  else {
    itData->samplePoints.reserve(n / 3 + itData->oldCenters.size());
    itData->samplePoints.clear();
    for (const auto &p : itData->oldCenters) {
      itData->samplePoints.emplace_back(p);
    }

    std::cout << "rr: ";
    for (int i = 0; i < n / 3; i++) {
      ES::V3d pp(x[i * 3], x[i * 3 + 1], x[i * 3 + 2]);
      itData->samplePoints.emplace_back(pp);
    }

    std::cout << n << ',' << itData->samplePoints.size() << std::endl;

    std::vector<ES::V3d> samplePoints2;
    samplePoints2.reserve(n / 3 + itData->oldCenters.size());

    std::set<int> visitedTriangles;

    for (int i = 0; i < (int)itData->samplePoints.size(); i++) {
      ES::V3d p = itData->samplePoints[i];

      double minDist = 1e100;
      ES::V3d closestPt;
      int triID = -1;
      for (int ti = 0; ti < (int)itData->medialAxisFacetsInitial.size(); ti++) {
        if (itData->medialAxisFacetsInitial[ti].size() == 2ull) {
          ES::V2d w;
          ES::V3d cpt = pgo::Mesh::getClosestPointToLineSegment(p,
            itData->medialAxisVerticesInitialIK[itData->medialAxisFacetsInitial[ti][0]],
            itData->medialAxisVerticesInitialIK[itData->medialAxisFacetsInitial[ti][1]], w);

          double d2 = (cpt - p).squaredNorm();
          if (d2 < minDist) {
            minDist = d2;
            closestPt = cpt;
            triID = ti;
          }
        }
        else if (itData->medialAxisFacetsInitial[ti].size() == 3ull) {
          ES::V3d w;
          ES::V3d cpt;
          int feature;
          double d2 = pgo::Mesh::getSquaredDistanceToTriangle(p,
            itData->medialAxisVerticesInitialIK[itData->medialAxisFacetsInitial[ti][0]],
            itData->medialAxisVerticesInitialIK[itData->medialAxisFacetsInitial[ti][1]],
            itData->medialAxisVerticesInitialIK[itData->medialAxisFacetsInitial[ti][2]],
            feature, cpt, w);

          if (d2 < minDist) {
            minDist = d2;
            closestPt = cpt;
            triID = ti;
          }
        }
      }
      PGO_ALOG(triID >= 0);

      if (visitedTriangles.find(triID) != visitedTriangles.end()) {
        return 1;
      }
      else {
        if (itData->medialAxisFacetsInitial[triID].size() == 2ull) {
          ES::V3d c = (itData->medialAxisVerticesInitialIK[itData->medialAxisFacetsInitial[triID][0]] +
                        itData->medialAxisVerticesInitialIK[itData->medialAxisFacetsInitial[triID][1]]) *
            0.5;
          samplePoints2.emplace_back(c);
        }
        else if (itData->medialAxisFacetsInitial[triID].size() == 3ull) {
          ES::V3d c = (itData->medialAxisVerticesInitialIK[itData->medialAxisFacetsInitial[triID][0]] +
                        itData->medialAxisVerticesInitialIK[itData->medialAxisFacetsInitial[triID][1]] +
                        itData->medialAxisVerticesInitialIK[itData->medialAxisFacetsInitial[triID][2]]) /
            3;
          samplePoints2.emplace_back(c);
        }
        else {
          abort();
        }

        visitedTriangles.emplace(triID);
      }
    }

    itData->samplePoints = samplePoints2;
  }

  pgo::Mesh::TriMeshGeo aa;
  for (int i = 0; i < itData->samplePoints.size(); i++) {
    aa.addPos(itData->samplePoints[i]);
  }
  // aa.save("pts.obj");

  double minDist = 1e100;
  for (int si = 0; si < (int)itData->samplePoints.size(); si++) {
    for (int sj = si + 1; sj < (int)itData->samplePoints.size(); sj++) {
      minDist = std::min(minDist, (itData->samplePoints[si] - itData->samplePoints[sj]).norm());
    }
  }

  if (minDist < 1e-2) {
    std::cout << "too close" << std::endl;
    return 1;
  }

  itData->medialAxisVerticesCur = itData->medialAxisVerticesInitial;
  itData->medialAxisFacetsCur = itData->medialAxisFacetsInitial;

  itData->finalSkeletonPoints.clear();
  itData->finalSkeletonEdges.clear();
  itData->finalSkeletonTriangles.clear();
  itData->finalFitMeshes.clear();
  itData->samplePoints0 = itData->samplePoints;

  computeRVD_internal_Dijkstras_withTri_refine_withEdge(itData->targetMesh, itData->targetMeshBVTree,
    itData->targetMeshSmall, itData->targetMeshBVTreeSmall, itData->samplePoints,
    itData->medialAxisVerticesCur, itData->medialAxisFacetsCur,
    itData->finalSkeletonPoints, itData->finalSkeletonEdges, itData->finalSkeletonTriangles,
    &itData->cellCenters);

  // to increase stability,
  // we can redo RVD
  if constexpr (1) {
    itData->medialAxisVerticesCur = itData->medialAxisVerticesInitial;
    itData->medialAxisFacetsCur = itData->medialAxisFacetsInitial;

    itData->samplePoints = itData->finalSkeletonPoints;
    itData->finalSkeletonPoints.clear();
    itData->finalSkeletonEdges.clear();
    itData->finalSkeletonTriangles.clear();
    itData->finalFitMeshes.clear();

    if constexpr (1) {
      std::vector<ES::V3d> samplePoints2;
      samplePoints2.reserve(n / 3 + itData->oldCenters.size());

      std::set<int> visitedTriangles;

      for (int i = 0; i < (int)itData->samplePoints.size(); i++) {
        ES::V3d p(itData->samplePoints[i].data());

        double minDist = 1e100;
        ES::V3d closestPt;
        int triID = -1;
        for (int ti = 0; ti < (int)itData->medialAxisFacetsInitial.size(); ti++) {
          if (itData->medialAxisFacetsInitial[ti].size() == 2ull) {
            ES::V2d w;
            ES::V3d cpt = pgo::Mesh::getClosestPointToLineSegment(p,
              itData->medialAxisVerticesInitialIK[itData->medialAxisFacetsInitial[ti][0]],
              itData->medialAxisVerticesInitialIK[itData->medialAxisFacetsInitial[ti][1]], w);

            double d2 = (cpt - p).squaredNorm();
            if (d2 < minDist) {
              minDist = d2;
              closestPt = cpt;
              triID = ti;
            }
          }
          else if (itData->medialAxisFacetsInitial[ti].size() == 3ull) {
            ES::V3d w;
            ES::V3d cpt;
            int feature;
            double d2 = pgo::Mesh::getSquaredDistanceToTriangle(p,
              itData->medialAxisVerticesInitialIK[itData->medialAxisFacetsInitial[ti][0]],
              itData->medialAxisVerticesInitialIK[itData->medialAxisFacetsInitial[ti][1]],
              itData->medialAxisVerticesInitialIK[itData->medialAxisFacetsInitial[ti][2]],
              feature, cpt, w);

            if (d2 < minDist) {
              minDist = d2;
              closestPt = cpt;
              triID = ti;
            }
          }
        }
        PGO_ALOG(triID >= 0);

        if (visitedTriangles.find(triID) != visitedTriangles.end()) {
          continue;
        }
        else {
          if (itData->medialAxisFacetsInitial[triID].size() == 2ull) {
            ES::V3d c = (itData->medialAxisVerticesInitialIK[itData->medialAxisFacetsInitial[triID][0]] +
                          itData->medialAxisVerticesInitialIK[itData->medialAxisFacetsInitial[triID][1]]) *
              0.5;
            samplePoints2.emplace_back(c);
          }
          else if (itData->medialAxisFacetsInitial[triID].size() == 3ull) {
            ES::V3d c = (itData->medialAxisVerticesInitialIK[itData->medialAxisFacetsInitial[triID][0]] +
                          itData->medialAxisVerticesInitialIK[itData->medialAxisFacetsInitial[triID][1]] +
                          itData->medialAxisVerticesInitialIK[itData->medialAxisFacetsInitial[triID][2]]) /
              3;
            samplePoints2.emplace_back(c);
          }
          else {
            abort();
          }

          visitedTriangles.emplace(triID);
        }
      }

      itData->samplePoints = samplePoints2;
    }

    computeRVD_internal_Dijkstras_withTri_refine_withEdge(itData->targetMesh, itData->targetMeshBVTree,
      itData->targetMeshSmall, itData->targetMeshBVTreeSmall, itData->samplePoints,
      itData->medialAxisVerticesCur, itData->medialAxisFacetsCur,
      itData->finalSkeletonPoints, itData->finalSkeletonEdges, itData->finalSkeletonTriangles,
      nullptr);
  }

  if constexpr (1) {
    std::vector<EK::Point_3> finalSkeletonPtEK(itData->finalSkeletonPoints.size());
    for (int i = 0; i < (int)itData->finalSkeletonPoints.size(); i++) {
      finalSkeletonPtEK[i] = EK::Point_3(itData->finalSkeletonPoints[i][0], itData->finalSkeletonPoints[i][1], itData->finalSkeletonPoints[i][2]);
    }

    std::set<std::pair<int, int>> edgeSet;
    for (auto e : itData->finalSkeletonEdges) {
      if (e.first > e.second) {
        std::swap(e.first, e.second);
      }
      edgeSet.emplace(e);
    }

    std::set<std::array<int, 3>, ES::IntArrayLess<3>> triSet;
    for (int ti = 0; ti < (int)itData->finalSkeletonTriangles.size(); ti++) {
      std::array<int, 3> vids{
        std::get<0>(itData->finalSkeletonTriangles[ti]),
        std::get<1>(itData->finalSkeletonTriangles[ti]),
        std::get<2>(itData->finalSkeletonTriangles[ti]),
      };
      std::sort(vids.begin(), vids.end());

      triSet.emplace(vids);
    }

    int edgeDel = 0;
    for (const auto &tri : triSet) {
      for (int j = 0; j < 3; j++) {
        std::pair edgeID{ tri[j], tri[(j + 1) % 3] };
        if (edgeID.first > edgeID.second) {
          std::swap(edgeID.first, edgeID.second);
        }

        auto it = edgeSet.find(edgeID);
        if (it != edgeSet.end()) {
          edgeSet.erase(it);
          edgeDel++;
        }
      }
    }

    itData->finalSkeletonEdges.assign(edgeSet.begin(), edgeSet.end());
    itData->finalSkeletonTriangles.clear();
    for (const auto &tri : triSet) {
      itData->finalSkeletonTriangles.emplace_back(tri[0], tri[1], tri[2]);
    }
    std::cout << "EE:" << edgeDel << std::endl;

    SkeletonExact maSk;
    maSk.initFromInput(finalSkeletonPtEK, itData->finalSkeletonEdges, itData->finalSkeletonTriangles);

    // remove redundant vertices that in the middle of a almost straight line.
    while (1) {
      bool modified = false;
      for (const auto &pr : maSk.vertices) {
        if (pr.second.neighboringEdges.size() == 2ull &&
          pr.second.neighboringTriangles.size() == 0ull) {
          EK::Point_3 p0 = pr.second.pos;
          ES::V3d pt = toVec3(p0);

          bool isSample = false;
          for (int i = 0; i < (int)itData->samplePoints.size(); i++) {
            if ((pt - itData->samplePoints[i]).norm() < 1e-6) {
              isSample = true;
              break;
            }
          }

          if (isSample) {
            continue;
          }

          const auto &e0 = maSk.getE(pr.second.neighboringEdges[0]);
          const auto &e1 = maSk.getE(pr.second.neighboringEdges[1]);

          int vids[2];
          if (e0.vidx[0] == pr.first) {
            vids[0] = e0.vidx[1];
          }
          else {
            vids[0] = e0.vidx[0];
          }

          if (e1.vidx[0] == pr.first) {
            vids[1] = e1.vidx[1];
          }
          else {
            vids[1] = e1.vidx[0];
          }

          EK::Point_3 p1 = maSk.getV(vids[0]).pos;
          EK::Point_3 p2 = maSk.getV(vids[1]).pos;

          EK::Vector_3 d0 = p1 - p0;
          EK::Vector_3 d1 = p2 - p0;

          double d0_norm = std::sqrt(CGAL::to_double(d0.squared_length()));
          double d1_norm = std::sqrt(CGAL::to_double(d1.squared_length()));
          double dot_d0_d1 = CGAL::to_double(CGAL::scalar_product(d0, d1));
          double cosAngle = dot_d0_d1 / (d0_norm * d1_norm);

          if (cosAngle < std::cos(179.0 / 180.0 * M_PI)) {
            continue;
          }

          ES::V3d v0 = toVec3(p1);
          ES::V3d v1 = toVec3(p2);

          if (itData->targetMeshBVTree.hasLineSegmentIntersectionExact(itData->targetMesh, v0, v1)) {
            continue;
          }

          if (maSk.canCollapseEdge(pr.second.neighboringEdges[0]) == false) {
            continue;
          }

          maSk.collapseEdge(pr.second.neighboringEdges[0], p1);
          modified = true;
          break;
        }
      }

      if (modified == false)
        break;
    }
  }

  // compute vertex pos
  std::vector<ES::V3d> maVertices(itData->medialAxisVerticesCur.size());
  for (int vi = 0; vi < (int)itData->medialAxisVerticesCur.size(); vi++) {
    maVertices[vi] = toVec3(itData->medialAxisVerticesCur[vi]);
  }

  // compute vertex weights
  std::vector<double> maVertexWeights(itData->medialAxisVerticesCur.size(), 0.0);
  for (int ti = 0; ti < (int)itData->medialAxisFacetsCur.size(); ti++) {
    double aaa = 0;
    if (itData->medialAxisFacetsCur[ti].size() == 3ull) {
      aaa = pgo::Mesh::getTriangleArea(maVertices[itData->medialAxisFacetsCur[ti][0]], maVertices[itData->medialAxisFacetsCur[ti][1]], maVertices[itData->medialAxisFacetsCur[ti][2]]);
      maVertexWeights[itData->medialAxisFacetsCur[ti][0]] += aaa / 3.0;
      maVertexWeights[itData->medialAxisFacetsCur[ti][1]] += aaa / 3.0;
      maVertexWeights[itData->medialAxisFacetsCur[ti][2]] += aaa / 3.0;
    }
    else if (itData->medialAxisFacetsCur[ti].size() == 2ull) {
      aaa = (maVertices[itData->medialAxisFacetsCur[ti][0]] - maVertices[itData->medialAxisFacetsCur[ti][1]]).norm();
      maVertexWeights[itData->medialAxisFacetsCur[ti][0]] += aaa / 2.0;
      maVertexWeights[itData->medialAxisFacetsCur[ti][1]] += aaa / 2.0;
    }
  }

  std::vector<int> covered(itData->targetMesh.numVertices(), 0);
  std::vector<int> skeletonCovered(itData->medialAxisVerticesCur.size(), 0);
  std::vector<int> tetCovered(itData->targetMeshTetMeshTetCenters.size(), 0);

  pgo::Mesh::TriMeshGeo fitAll;
  bool isAllSucceed = true;
  int fitMeshID = 0;

  // fit sphere
  for (int si = 0; si < (int)itData->finalSkeletonPoints.size(); si++) {
    // fit sphere
    ES::MXd centers(1, 3);
    ES::VXd centerRadii(1);
    centers.row(0) = itData->finalSkeletonPoints[si];

    ES::V3d p = centers.row(0);
    auto ret = itData->targetMeshBVTree.closestTriangleQuery(itData->targetMesh, p);
    centerRadii(0) = std::max(std::sqrt(ret.dist2), 0.01);

    MedialAxisRepresentation::PrimitiveICPSolverParameters primitiveICPParams;
    primitiveICPParams.sphereRadius = centerRadii(0);
    primitiveICPParams.useCorefine = false;
    primitiveICPParams.verbose = false;
    primitiveICPParams.maxPDNumIter = 10;
    primitiveICPParams.maxNumIter = maxIter;
    primitiveICPParams.smoothnessCoeff = sc;
    primitiveICPParams.expansionCoeff = 2;
    primitiveICPParams.contactCoeff = 1e4;

    // std::ofstream(fmt::format("{}/curr.txt", saveDir)) << centers(0, 0) << ',' << centers(0, 1) << ',' << centers(0, 2) << std::endl;

    pgo::Mesh::TriMeshGeo fittingMesh;
    MedialAxisRepresentation::PrimitiveICPSolver solver(centers, centerRadii, primitiveICPParams, 1);
    int solverRet = solver.sim(centers, itData->targetMesh, itData->targetMeshBVTree, itData->targetNormals, fittingMesh);
    if (solverRet != 0) {
      isAllSucceed = false;
      continue;
    }

    pgo::Mesh::TriMeshBVTree fittingMeshBVTree;
    fittingMeshBVTree.buildByInertiaPartition(fittingMesh);

    // compute surface vertex coverage
    for (int vi = 0; vi < itData->targetMesh.numVertices(); vi++) {
      ES::V3d diff = itData->targetMesh.pos(vi) - itData->finalSkeletonPoints[si];
      double diff_len = diff.norm();
      double shrunk_len = std::max(diff_len - itData->avgLength, 1e-4);

      if (diff_len < 1e-6) {
        covered[vi] = 1;
        continue;
      }

      diff = shrunk_len / diff_len * diff;
      ES::V3d start = itData->finalSkeletonPoints[si];
      ES::V3d end = start + diff;

      if (fittingMeshBVTree.hasLineSegmentIntersectionExact(fittingMesh, start, end) == false) {
        covered[vi] = 1;
      }
    }

    // fittingMesh.save(fmt::format("{}/{}.obj", saveDir, fitMeshID));
    itData->finalFitMeshes.push_back(fittingMesh);
    fitMeshID++;
    fitAll.addMesh(fittingMesh);

    if (itData->computeUnmaskedRegion) {
      // compute skeleton vertex coverage
      for (int vi = 0; vi < (int)maVertices.size(); vi++) {
        ES::V3d pt = maVertices[vi];
        ES::V3d diff = pt - itData->finalSkeletonPoints[si];
        double diff_len = diff.norm();
        double shrunk_len = std::max(diff_len - itData->avgLength, 1e-4);

        if (diff_len < 1e-6) {
          skeletonCovered[vi] = 1;
          continue;
        }

        diff = shrunk_len / diff_len * diff;
        ES::V3d start = itData->finalSkeletonPoints[si];
        ES::V3d end = start + diff;

        if (fittingMeshBVTree.hasLineSegmentIntersectionExact(fittingMesh, start, end) == false) {
          skeletonCovered[vi] = 1;
        }
      }

      // compute tet mesh coverage
      for (int vi = 0; vi < (int)itData->targetMeshTetMeshTetCenters.size(); vi++) {
        ES::V3d pt = itData->targetMeshTetMeshTetCenters[vi];
        ES::V3d diff = pt - itData->finalSkeletonPoints[si];
        double diff_len = diff.norm();
        double shrunk_len = std::max(diff_len - itData->avgLength, 1e-4);

        if (diff_len < 1e-6) {
          tetCovered[vi] = 1;
          continue;
        }

        diff = shrunk_len / diff_len * diff;
        ES::V3d start = itData->finalSkeletonPoints[si];
        ES::V3d end = start + diff;

        if (fittingMeshBVTree.hasLineSegmentIntersectionExact(fittingMesh, start, end) == false) {
          tetCovered[vi] = 1;
        }
      }
    }
  }

  // fit edge
  for (int ei = 0; ei < (int)itData->finalSkeletonEdges.size(); ei++) {
    ES::MXd centers(2, 3);
    ES::VXd centerRadii(2);

    centers.row(0) = itData->finalSkeletonPoints[itData->finalSkeletonEdges[ei].first];
    centers.row(1) = itData->finalSkeletonPoints[itData->finalSkeletonEdges[ei].second];

    ES::V3d diff = centers.row(0) - centers.row(1);
    if (diff.norm() < 1e-4) {
      continue;
    }

    ES::V3d p0 = itData->finalSkeletonPoints[itData->finalSkeletonEdges[ei].first];
    ES::V3d p1 = itData->finalSkeletonPoints[itData->finalSkeletonEdges[ei].second];

    // std::cout << p0 << std::endl;
    // std::cout << p1 << std::endl;
    // std::cout << len(p0 - p1) << std::endl;

    if (itData->targetMeshBVTree.hasLineSegmentIntersectionExact(itData->targetMesh, p0, p1)) {
      continue;
    }

    ES::V3d p = centers.row(0);
    auto ret = itData->targetMeshBVTree.closestTriangleQuery(itData->targetMesh, p);
    centerRadii(0) = std::sqrt(ret.dist2);

    p = centers.row(1);
    ret = itData->targetMeshBVTree.closestTriangleQuery(itData->targetMesh, p);
    centerRadii(1) = std::sqrt(ret.dist2);

    // std::ofstream(fmt::format("{}/curr.txt", saveDir)) << centers(0, 0) << ',' << centers(0, 1) << ',' << centers(0, 2) << std::endl
    //                                                    << centers(1, 0) << ',' << centers(1, 1) << ',' << centers(1, 2) << std::endl;

    MedialAxisRepresentation::PrimitiveICPSolverParameters primitiveICPParams;
    primitiveICPParams.sphereRadius = (centerRadii(0) + centerRadii(1)) * 0.5;
    primitiveICPParams.useCorefine = false;
    primitiveICPParams.verbose = false;
    primitiveICPParams.maxPDNumIter = 15;
    primitiveICPParams.maxNumIter = maxIter;
    primitiveICPParams.smoothnessCoeff = sc;
    primitiveICPParams.expansionCoeff = 2;
    primitiveICPParams.contactCoeff = 1e4;

    pgo::Mesh::TriMeshGeo fittingMesh;
    MedialAxisRepresentation::PrimitiveICPSolver solver(centers, centerRadii, primitiveICPParams, 2);
    int solverRet = solver.sim(centers, itData->targetMesh, itData->targetMeshBVTree, itData->targetNormals, fittingMesh);
    if (solverRet != 0) {
      isAllSucceed = false;
      continue;
    }

    // fittingMesh.save("fit.obj");

    pgo::Mesh::TriMeshBVTree fittingMeshBVTree;
    fittingMeshBVTree.buildByInertiaPartition(fittingMesh);

    // check surface coverage
    // dT (x0 + td - x1) = 0
    // t dTd + dT (x0 - x1) = 0
    for (int vi = 0; vi < itData->targetMesh.numVertices(); vi++) {
      ES::V3d x0 = itData->finalSkeletonPoints[itData->finalSkeletonEdges[ei].first];
      ES::V3d x1 = itData->finalSkeletonPoints[itData->finalSkeletonEdges[ei].second];
      ES::V3d dir = x1 - x0;

      ES::V3d pt = itData->targetMesh.pos(vi);
      double a = dir.dot(dir);
      double b = (pt - x0).dot(dir);
      if (std::abs(a) < 1e-6) {
        continue;
      }

      double t = b / a;
      if (t < 0 || t > 1) {
        continue;
      }

      ES::V3d center = dir * t + x0;
      ES::V3d diff = itData->targetMesh.pos(vi) - center;
      double diff_len = diff.norm();
      double shrunk_len = std::max(diff_len - itData->avgLength, 1e-4);

      if (diff_len < 1e-6) {
        covered[vi] = 1;
        continue;
      }

      if (shrunk_len < 1e-6) {
        covered[vi] = 1;
        continue;
      }

      diff = shrunk_len / diff_len * diff;
      ES::V3d start = center;
      ES::V3d end = diff + start;

      if (fittingMeshBVTree.hasLineSegmentIntersectionExact(fittingMesh, start, end) == false) {
        covered[vi] = 1;
      }
    }

    // fittingMesh.save(fmt::format("{}/{}.obj", saveDir, fitMeshID));
    itData->finalFitMeshes.push_back(fittingMesh);
    fitMeshID++;
    fitAll.addMesh(fittingMesh);

    if (itData->computeUnmaskedRegion) {
      for (int vi = 0; vi < (int)itData->medialAxisVerticesCur.size(); vi++) {
        ES::V3d pt = maVertices[vi];
        ES::V3d x0 = itData->finalSkeletonPoints[itData->finalSkeletonEdges[ei].first];
        ES::V3d x1 = itData->finalSkeletonPoints[itData->finalSkeletonEdges[ei].second];
        ES::V3d dir = x1 - x0;
        double a = dir.dot(dir);
        double b = (pt - x0).dot(dir);
        if (std::abs(a) < 1e-6) {
          continue;
        }

        double t = b / a;
        if (t < 0 || t > 1) {
          continue;
        }

        ES::V3d center = dir * t + x0;
        ES::V3d diff = pt - center;
        double diff_len = diff.norm();
        double shrunk_len = std::max(diff_len - itData->avgLength, 1e-4);

        if (diff_len < 1e-6) {
          skeletonCovered[vi] = 1;
          continue;
        }

        if (shrunk_len < 1e-6) {
          skeletonCovered[vi] = 1;
          continue;
        }

        diff = shrunk_len / diff_len * diff;
        ES::V3d start = center;
        ES::V3d end = diff + start;

        if (fittingMeshBVTree.hasLineSegmentIntersectionExact(fittingMesh, start, end) == false) {
          skeletonCovered[vi] = 1;
        }
      }

      for (int vi = 0; vi < (int)itData->targetMeshTetMeshTetCenters.size(); vi++) {
        ES::V3d pt = itData->targetMeshTetMeshTetCenters[vi];
        ES::V3d x0 = itData->finalSkeletonPoints[itData->finalSkeletonEdges[ei].first];
        ES::V3d x1 = itData->finalSkeletonPoints[itData->finalSkeletonEdges[ei].second];
        ES::V3d dir = x1 - x0;
        double a = dir.dot(dir);
        double b = (pt - x0).dot(dir);
        if (std::abs(a) < 1e-6) {
          continue;
        }

        double t = b / a;
        if (t < 0 || t > 1) {
          continue;
        }

        ES::V3d center = dir * t + x0;
        ES::V3d diff = pt - center;
        double diff_len = diff.norm();
        double shrunk_len = std::max(diff_len - itData->avgLength, 1e-4);

        if (diff_len < 1e-6) {
          tetCovered[vi] = 1;
          continue;
        }

        if (shrunk_len < 1e-6) {
          tetCovered[vi] = 1;
          continue;
        }

        diff = shrunk_len / diff_len * diff;
        ES::V3d start = center;
        ES::V3d end = diff + start;

        if (fittingMeshBVTree.hasLineSegmentIntersectionExact(fittingMesh, start, end) == false) {
          tetCovered[vi] = 1;
        }
      }
    }
  }

  // std::cout << "****************************" << std::endl;

  // fit triangles
  for (int ei = 0; ei < (int)itData->finalSkeletonTriangles.size(); ei++) {
    ES::MXd centers(3, 3);
    ES::VXd centerRadii(3);

    centers.row(0) = itData->finalSkeletonPoints[std::get<0>(itData->finalSkeletonTriangles[ei])];
    centers.row(1) = itData->finalSkeletonPoints[std::get<1>(itData->finalSkeletonTriangles[ei])];
    centers.row(2) = itData->finalSkeletonPoints[std::get<2>(itData->finalSkeletonTriangles[ei])];

    ES::V3d e01 = centers.row(1) - centers.row(0);
    ES::V3d e02 = centers.row(2) - centers.row(0);
    if (e01.cross(e02).norm() < 1e-4) {
      continue;
    }

    ES::V3d p0 = itData->finalSkeletonPoints[std::get<0>(itData->finalSkeletonTriangles[ei])];
    ES::V3d p1 = itData->finalSkeletonPoints[std::get<1>(itData->finalSkeletonTriangles[ei])];
    ES::V3d p2 = itData->finalSkeletonPoints[std::get<2>(itData->finalSkeletonTriangles[ei])];

    std::vector<int> l;
    itData->targetMeshBVTree.triangleIntersectionExact(itData->targetMesh, p0, p1, p2, l);
    if (l.size()) {
      continue;
    }

    for (int j = 0; j < 3; j++) {
      ES::V3d p = centers.row(j);
      auto ret = itData->targetMeshBVTree.closestTriangleQuery(itData->targetMesh, p);
      centerRadii(j) = std::sqrt(ret.dist2);
    }

    // std::ofstream(fmt::format("{}/curr.txt", saveDir)) << centers(0, 0) << ',' << centers(0, 1) << ',' << centers(0, 2) << std::endl
    //                                                    << centers(1, 0) << ',' << centers(1, 1) << ',' << centers(1, 2) << std::endl
    //                                                    << centers(2, 0) << ',' << centers(2, 1) << ',' << centers(2, 2) << std::endl;

    MedialAxisRepresentation::PrimitiveICPSolverParameters primitiveICPParams;
    primitiveICPParams.sphereRadius = (centerRadii(0) + centerRadii(1) + centerRadii(2)) / 3;
    primitiveICPParams.useCorefine = false;
    primitiveICPParams.verbose = false;
    primitiveICPParams.maxPDNumIter = 15;
    primitiveICPParams.maxNumIter = maxIter;
    primitiveICPParams.smoothnessCoeff = sc;
    primitiveICPParams.expansionCoeff = 2;
    primitiveICPParams.contactCoeff = 1e4;

    pgo::Mesh::TriMeshGeo fittingMesh;
    MedialAxisRepresentation::PrimitiveICPSolver solver(centers, centerRadii, primitiveICPParams, 3);
    int solverRet = solver.sim(centers, itData->targetMesh, itData->targetMeshBVTree, itData->targetNormals, fittingMesh);
    if (solverRet != 0) {
      isAllSucceed = false;
      continue;
    }

    // fittingMesh.save("fit.obj");

    pgo::Mesh::TriMeshBVTree fittingMeshBVTree;
    fittingMeshBVTree.buildByInertiaPartition(fittingMesh);

    // dT (x0 + td - x1) = 0
    // t dTd + dT (x0 - x1) = 0
    for (int vi = 0; vi < itData->targetMesh.numVertices(); vi++) {
      ES::V3d x[3] = {
        itData->finalSkeletonPoints[std::get<0>(itData->finalSkeletonTriangles[ei])],
        itData->finalSkeletonPoints[std::get<1>(itData->finalSkeletonTriangles[ei])],
        itData->finalSkeletonPoints[std::get<2>(itData->finalSkeletonTriangles[ei])],
      };

      ES::V3d p = itData->targetMesh.pos(vi);

      int feature;
      ES::V3d closestPt;
      ES::V3d w;
      pgo::Mesh::getSquaredDistanceToTriangle(p, x[0], x[1], x[2], feature, closestPt, w);

      double dist = (p - closestPt).norm();
      if (dist < 1e-6) {
        covered[vi] = 1;
        continue;
      }

      ES::V3d dir = (p - closestPt).normalized();
      if (feature >= 0 && feature < 3) {
        ES::V3d e0 = (x[(feature + 1) % 3] - x[feature]).normalized();
        ES::V3d e1 = (x[(feature + 2) % 3] - x[feature]).normalized();

        if (std::abs(dir.dot(e0)) < 1e-6 ||
          std::abs(dir.dot(e1)) < 1e-6) {
          double shrunk_len = std::max(dist - itData->avgLength, 1e-4);
          if (shrunk_len < 1e-6) {
            covered[vi] = 1;
            continue;
          }

          ES::V3d start = closestPt;
          ES::V3d end = closestPt + dir * shrunk_len;

          if (fittingMeshBVTree.hasLineSegmentIntersectionExact(fittingMesh, start, end) == false) {
            covered[vi] = 1;
          }
        }
      }
      else {
        double shrunk_len = std::max(dist - itData->avgLength, 1e-4);
        if (shrunk_len < 1e-6) {
          covered[vi] = 1;
          continue;
        }

        ES::V3d start = closestPt;
        ES::V3d end = closestPt + dir * shrunk_len;

        if (fittingMeshBVTree.hasLineSegmentIntersectionExact(fittingMesh, start, end) == false) {
          covered[vi] = 1;
        }
      }
    }

    // fittingMesh.save(fmt::format("{}/{}.obj", saveDir, fitMeshID));
    itData->finalFitMeshes.push_back(fittingMesh);
    fitMeshID++;
    fitAll.addMesh(fittingMesh);

    if (itData->computeUnmaskedRegion) {
      for (int vi = 0; vi < (int)itData->medialAxisVerticesCur.size(); vi++) {
        ES::V3d p = maVertices[vi];
        ES::V3d x[3] = {
          itData->finalSkeletonPoints[std::get<0>(itData->finalSkeletonTriangles[ei])],
          itData->finalSkeletonPoints[std::get<1>(itData->finalSkeletonTriangles[ei])],
          itData->finalSkeletonPoints[std::get<2>(itData->finalSkeletonTriangles[ei])],
        };

        int feature;
        ES::V3d closestPt;
        ES::V3d w;
        pgo::Mesh::getSquaredDistanceToTriangle(p, x[0], x[1], x[2], feature, closestPt, w);

        double dist = (p - closestPt).norm();
        if (dist < 1e-6) {
          skeletonCovered[vi] = 1;
          continue;
        }

        ES::V3d dir = (p - closestPt).normalized();
        if (feature >= 0 && feature < 3) {
          ES::V3d e0 = (x[(feature + 1) % 3] - x[feature]).normalized();
          ES::V3d e1 = (x[(feature + 2) % 3] - x[feature]).normalized();

          if (std::abs(dir.dot(e0)) < 1e-6 ||
            std::abs(dir.dot(e1)) < 1e-6) {
            double shrunk_len = std::max(dist - itData->avgLength, 1e-4);
            if (shrunk_len < 1e-6) {
              skeletonCovered[vi] = 1;
              continue;
            }

            ES::V3d start = closestPt;
            ES::V3d end = closestPt + dir * shrunk_len;

            if (fittingMeshBVTree.hasLineSegmentIntersectionExact(fittingMesh, start, end) == false) {
              skeletonCovered[vi] = 1;
            }
          }
        }
        else {
          double shrunk_len = std::max(dist - itData->avgLength, 1e-4);
          if (shrunk_len < 1e-6) {
            skeletonCovered[vi] = 1;
            continue;
          }

          ES::V3d start = closestPt;
          ES::V3d end = closestPt + dir * shrunk_len;

          if (fittingMeshBVTree.hasLineSegmentIntersectionExact(fittingMesh, start, end) == false) {
            skeletonCovered[vi] = 1;
          }
        }
      }

      for (int vi = 0; vi < (int)itData->targetMeshTetMeshTetCenters.size(); vi++) {
        ES::V3d p = itData->targetMeshTetMeshTetCenters[vi];
        ES::V3d x[3] = {
          itData->finalSkeletonPoints[std::get<0>(itData->finalSkeletonTriangles[ei])],
          itData->finalSkeletonPoints[std::get<1>(itData->finalSkeletonTriangles[ei])],
          itData->finalSkeletonPoints[std::get<2>(itData->finalSkeletonTriangles[ei])],
        };

        int feature;
        ES::V3d closestPt;
        ES::V3d w;
        pgo::Mesh::getSquaredDistanceToTriangle(p, x[0], x[1], x[2], feature, closestPt, w);

        double dist = (p - closestPt).norm();
        if (dist < 1e-6) {
          tetCovered[vi] = 1;
          continue;
        }

        ES::V3d dir = (p - closestPt).normalized();
        if (feature >= 0 && feature < 3) {
          ES::V3d e0 = (x[(feature + 1) % 3] - x[feature]).normalized();
          ES::V3d e1 = (x[(feature + 2) % 3] - x[feature]).normalized();

          if (std::abs(dir.dot(e0)) < 1e-6 ||
            std::abs(dir.dot(e1)) < 1e-6) {
            double shrunk_len = std::max(dist - itData->avgLength, 1e-4);
            if (shrunk_len < 1e-6) {
              tetCovered[vi] = 1;
              continue;
            }

            ES::V3d start = closestPt;
            ES::V3d end = closestPt + dir * shrunk_len;

            if (fittingMeshBVTree.hasLineSegmentIntersectionExact(fittingMesh, start, end) == false) {
              tetCovered[vi] = 1;
            }
          }
        }
        else {
          double shrunk_len = std::max(dist - itData->avgLength, 1e-4);
          if (shrunk_len < 1e-6) {
            tetCovered[vi] = 1;
            continue;
          }

          ES::V3d start = closestPt;
          ES::V3d end = closestPt + dir * shrunk_len;

          if (fittingMeshBVTree.hasLineSegmentIntersectionExact(fittingMesh, start, end) == false) {
            tetCovered[vi] = 1;
          }
        }
      }
    }
  }

  std::cout << "=================================" << std::endl;

  std::cout << "Iter: " << itData->iter << ":\n";
  PGO_ALOG(fitMeshID == (int)itData->finalFitMeshes.size());

  if (isAllSucceed) {
    std::cout << "All fitting is working well." << std::endl;
  }
  else {
    std::cout << "Some fitting failed." << std::endl;
    // static int failed_id = 0;
    // pgo::Mesh::TriMeshGeo outSampleMesh;

    // for (const auto &p : itData->samplePoints) {
    //   outSampleMesh.addPos(p);
    // }
    // outSampleMesh.save(fmt::format("failed{:05d}.obj", failed_id++));
  }

  if (itData->computeUnmaskedRegion) {
    // pgo::Mesh::TriMeshGeo rr;
    // for (int vi = 0; vi < (int)itData->targetMeshTetMeshTetCenters.size(); vi++) {
    //   if (tetCovered[vi])
    //     continue;

    //  rr.addPos(itData->targetMeshTetMeshTetCenters[vi]);
    //}
    // rr.save("zz.obj");

    std::vector<std::vector<int>> unmaskedCCs;
    std::vector<int> visitedTets(itData->targetMeshTetMeshTetCenters.size(), 0);

    while (1) {
      int selID = -1;
      for (int vi = 0; vi < (int)itData->targetMeshTetMeshTetCenters.size(); vi++) {
        if (tetCovered[vi])
          continue;

        if (visitedTets[vi])
          continue;

        selID = vi;
        break;
      }

      if (selID < 0)
        break;

      std::vector<int> Q;
      Q.reserve(itData->targetMeshTetMeshTetCenters.size());
      Q.emplace_back(selID);
      visitedTets[selID] = 1;

      size_t start = 0;
      while (start < Q.size()) {
        int curID = Q[start++];

        for (int nti : itData->targetMeshTetMeshTetNeighbors[curID]) {
          if (visitedTets[nti])
            continue;

          if (tetCovered[nti])
            continue;

          Q.emplace_back(nti);
          visitedTets[nti] = 1;
        }
      }

      unmaskedCCs.emplace_back(std::move(Q));
    }

    std::cout << "# unmasked ccs: " << unmaskedCCs.size() << std::endl;
    std::vector<std::pair<double, int>> ccSpaces(unmaskedCCs.size());
    for (int i = 0; i < (int)unmaskedCCs.size(); i++) {
      double wAll = 0;
      for (int teti : unmaskedCCs[i]) {
        wAll += pgo::Mesh::getTetVolume(itData->targetMeshTetMesh.pos(teti, 0),
          itData->targetMeshTetMesh.pos(teti, 1),
          itData->targetMeshTetMesh.pos(teti, 2),
          itData->targetMeshTetMesh.pos(teti, 3));
      }

      ccSpaces[i].first = wAll;
      ccSpaces[i].second = i;
    }
    std::sort(ccSpaces.begin(), ccSpaces.end());

    itData->newCenters.clear();
    for (int i = (int)ccSpaces.size() - 1, j = itData->nAddedPts - 1; i >= 0 && j >= 0; i--, j--) {
      const auto &cc = unmaskedCCs[ccSpaces[i].second];
      ES::V3d center(0, 0, 0);

      pgo::Mesh::TriMeshGeo aa;
      for (int vi : cc) {
        center += itData->targetMeshTetMeshTetCenters[vi];
        aa.addPos(itData->targetMeshTetMeshTetCenters[vi]);
      }

      center /= (double)cc.size();
      itData->newCenters.emplace_back(center[0], center[1], center[2]);

      // static int z = 0;
      // aa.save(fmt::format("r{}.obj", z++));
    }

    std::cout << "# filtered ccs: " << itData->newCenters.size() << std::endl;
  }

  // int cc = std::count_if(covered.begin(), covered.end(), [](int v) { return v == 1; });
  double aAll = 0, aCovered = 0;
  for (int i = 0; i < (int)covered.size(); i++) {
    aAll += itData->targetVertexWeights[i];
    aCovered += itData->targetVertexWeights[i] * covered[i];
  }

  std::cout << aCovered << "/" << aAll << std::endl;

  // fitAll.save(fmt::format("{}/f{}.obj", saveDir, itData->iter));
  // fitAll.save("latest-fit.obj");

  itData->iter += 1;

  double cellDiff = 0;
  PGO_ALOG(n / 3 == int(itData->cellCenters.size() - itData->oldCenters.size()));
  for (int i = 0; i < n / 3; i++) {
    cellDiff += (itData->cellCenters[i + (int)itData->oldCenters.size()] - ES::V3d(x[i * 3], x[i * 3 + 1], x[i * 3 + 2])).squaredNorm();
  }
  cellDiff /= (n / 3);

  double E1 = -aCovered / aAll;
  double E2 = cellDiff * 1;
  double E3 = std::abs(double(itData->finalSkeletonPoints.size() - n / 3)) * 0.001;

  std::cout << "E1:" << E1 << "; E2:" << E2 << "; E3:" << E3 << std::endl;
  std::cout << itData->finalSkeletonPoints.size() << ',' << (n / 3) << std::endl;

  // static std::ofstream outfile("energy.txt");
  // outfile << itData->iter << ',' << itData->computeUnmaskedRegion << ',' << E1 << ',' << E2 << ',' << E3 << ',' << (E1 + E2 + E3) << std::endl;

  return E1 + E2 + E3;
}
}  // namespace MedialAxisRepresentation

void MedialAxisRepresentation::solveSkeleton(
  const pgo::Mesh::TriMeshGeo &targetMesh, const pgo::Mesh::TriMeshBVTree &targetMeshBVTree, const pgo::Mesh::TriMeshPseudoNormal &targetNormals,
  const pgo::Mesh::TriMeshGeo &targetMeshSmall, const pgo::Mesh::TriMeshBVTree &targetMeshBVTreeSmall, const pgo::Mesh::TriMeshPseudoNormal &targetNormalsSmall,
  const std::string &maInitFilename, int nPt, int nIt, int numAddedPt,
  std::vector<ES::V3d> &finalSkeletonPoints,
  std::vector<std::pair<int, int>> &finalSkeletonEdges,
  std::vector<std::tuple<int, int, int>> &finalSkeletonTriangles,
  const std::string &userInitFilename,
  std::vector<pgo::Mesh::TriMeshGeo> *finalFitMeshes)
{
  IterationData *data = new IterationData(targetMesh, targetMeshBVTree, targetNormals,
    targetMeshSmall, targetMeshBVTreeSmall, targetNormalsSmall);
  data->nAddedPts = numAddedPt;

  pgo::Mesh::BoundingBox targetMeshBB(targetMesh.positions());
  int numPt = nPt;

  // compute volume coverage mesh
  if constexpr (1) {
    ES::MXd vtx;
    ES::MXi tet;
    pgo::TetgenInterface::computeTetMesh(targetMesh, "pq1.05", vtx, tet);

    std::vector<ES::V3d> tetVtx;
    std::vector<ES::V4i> tetEles;

    for (int vi = 0; vi < (int)vtx.rows(); vi++) {
      tetVtx.emplace_back(vtx(vi, 0), vtx(vi, 1), vtx(vi, 2));
    }

    for (int ti = 0; ti < (int)tet.rows(); ti++) {
      tetEles.emplace_back(tet(ti, 0), tet(ti, 1), tet(ti, 2), tet(ti, 3));
    }

    data->targetMeshTetMesh = pgo::Mesh::TetMeshGeo(tetVtx, tetEles);
    data->targetMeshTetMeshTetCenters.resize(tetEles.size());

    for (int ti = 0; ti < (int)tetEles.size(); ti++) {
      data->targetMeshTetMeshTetCenters[ti].setZero();
      data->targetMeshTetMeshTetCenters[ti] = (tetVtx[tetEles[ti][0]] + tetVtx[tetEles[ti][1]] + tetVtx[tetEles[ti][2]] + tetVtx[tetEles[ti][3]]) * 0.25;
    }

    std::map<pgo::Mesh::UTriKey, std::array<int, 10>> tetFaces;
    data->targetMeshTetMeshTetNeighbors.resize(tetEles.size());
    for (int ti = 0; ti < (int)tetEles.size(); ti++) {
      int faceIdx[][3] = {
        { 0, 1, 2 },
        { 0, 1, 3 },
        { 0, 2, 3 },
        { 1, 2, 3 },
      };

      for (int fi = 0; fi < 4; fi++) {
        pgo::Mesh::UTriKey triIdx(tetEles[ti][faceIdx[fi][0]], tetEles[ti][faceIdx[fi][1]], tetEles[ti][faceIdx[fi][2]]);
        auto it = tetFaces.find(triIdx);
        if (it == tetFaces.end()) {
          std::array<int, 10> tets;
          tets[0] = ti;
          tets.back() = 1;

          tetFaces.emplace(triIdx, tets);
        }
        else {
          it->second[it->second.back()] = ti;
          it->second.back()++;
        }
      }
    }

    for (const auto &it : tetFaces) {
      for (int i = 0; i < it.second.back(); i++) {
        for (int j = i + 1; j < it.second.back(); j++) {
          data->targetMeshTetMeshTetNeighbors[it.second[i]].emplace_back(it.second[j]);
          data->targetMeshTetMeshTetNeighbors[it.second[j]].emplace_back(it.second[i]);
        }
      }
    }

    for (int i = 0; i < (int)data->targetMeshTetMeshTetNeighbors.size(); i++) {
      pgo::BasicAlgorithms::sortAndDeduplicateWithErase(data->targetMeshTetMeshTetNeighbors[i]);
    }

    std::cout << "done tet" << std::endl;
  }

  ES::VXd x(numPt * 3);
  ES::VXd xlow(numPt * 3), xhi(numPt * 3);

  ES::VXd xtol(numPt * 3);
  xtol.setConstant(1e-2);

  for (int vi = 0; vi < numPt; vi++) {
    xlow.segment<3>(vi * 3) = targetMeshBB.bmin();
    xhi.segment<3>(vi * 3) = targetMeshBB.bmax();
  }

  // load initial medial axis
  if constexpr (1) {
    SkeletonExact ma;
    if (ma.load(maInitFilename.c_str()) != 0) {
      std::cerr << "Cannot load MA " << maInitFilename << std::endl;
      exit(1);
    }

    std::map<int, int> vidMap;
    for (const auto &pr : ma.vertices) {
      data->medialAxisVerticesInitial.emplace_back(pr.second.pos);
      vidMap.emplace(pr.first, (int)data->medialAxisVerticesInitial.size() - 1);
    }

    for (const auto &pr : ma.edges) {
      auto it0 = vidMap.find(pr.second.vidx[0]);
      PGO_ALOG(it0 != vidMap.end());

      auto it1 = vidMap.find(pr.second.vidx[1]);
      PGO_ALOG(it1 != vidMap.end());
      data->medialAxisFacetsInitial.emplace_back(std::vector<int>{ it0->second, it1->second });
    }

    for (const auto &pr : ma.triangles) {
      auto it0 = vidMap.find(pr.second.vidx[0]);
      PGO_ALOG(it0 != vidMap.end());

      auto it1 = vidMap.find(pr.second.vidx[1]);
      PGO_ALOG(it1 != vidMap.end());

      auto it2 = vidMap.find(pr.second.vidx[2]);
      PGO_ALOG(it2 != vidMap.end());
      data->medialAxisFacetsInitial.emplace_back(std::vector<int>{ it0->second, it1->second, it2->second });
    }
  }

  std::cout << "Done loading ma." << std::endl;

  std::vector<EK::Point_3> v;
  std::vector<std::vector<int>> f;
  preprocessingMA(data->medialAxisVerticesInitial, data->medialAxisFacetsInitial, targetMesh, targetMeshBVTree, v, f);

  data->medialAxisVerticesInitial = v;
  data->medialAxisFacetsInitial = f;

  std::vector<ES::V3d> maVertices(data->medialAxisVerticesInitial.size());
  for (int vi = 0; vi < (int)data->medialAxisVerticesInitial.size(); vi++) {
    maVertices[vi] = toVec3(data->medialAxisVerticesInitial[vi]);
  }

  std::vector<std::pair<double, int>> curvatures(maVertices.size());
  std::vector<std::pair<double, int>> distances(maVertices.size());

  ES::VXd vtxH;
  pgo::libiglInterface::meanCuravtures(targetMesh, vtxH);

  tbb::parallel_for(0, (int)maVertices.size(), [&](int vi) {
    auto ret = targetMeshBVTree.closestTriangleQuery(targetMesh, maVertices[vi]);
    double threshold = std::sqrt(ret.dist2) + 0.1;
    double hAll = 0;
    int c = 0;
    for (int ti = 0; ti < targetMesh.numTriangles(); ti++) {
      int feature;
      ES::V3d cpt;
      ES::V3d w;
      double d2 = pgo::Mesh::getSquaredDistanceToTriangle(maVertices[vi],
        targetMesh.pos(ti, 0), targetMesh.pos(ti, 1), targetMesh.pos(ti, 2), feature, cpt, w);

      if (d2 < threshold * threshold) {
        hAll += vtxH(targetMesh.triVtxID(ti, 0)) * w[0] + vtxH(targetMesh.triVtxID(ti, 1)) * w[1] + vtxH(targetMesh.triVtxID(ti, 2)) * w[2];
        c += 1;
      }
    }
    PGO_ALOG(c > 0);
    hAll /= c;

    curvatures[vi].first = hAll;
    curvatures[vi].second = vi;

    distances[vi].first = ret.dist2;
    distances[vi].second = vi;
  });

  std::sort(curvatures.begin(), curvatures.end());
  std::sort(distances.begin(), distances.end());

  std::vector<std::vector<int>> vertexNeighboringVertices(data->medialAxisVerticesInitial.size());
  for (int ti = 0; ti < (int)data->medialAxisFacetsInitial.size(); ti++) {
    for (int j = 0; j < (int)data->medialAxisFacetsInitial[ti].size(); j++) {
      vertexNeighboringVertices[data->medialAxisFacetsInitial[ti][j]].emplace_back(data->medialAxisFacetsInitial[ti][(j + 1) % (int)data->medialAxisFacetsInitial[ti].size()]);
      vertexNeighboringVertices[data->medialAxisFacetsInitial[ti][(j + 1) % (int)data->medialAxisFacetsInitial[ti].size()]].emplace_back(data->medialAxisFacetsInitial[ti][j]);
    }
  }

  for (int vi = 0; vi < (int)vertexNeighboringVertices.size(); vi++) {
    pgo::BasicAlgorithms::sortAndDeduplicateWithErase(vertexNeighboringVertices[vi]);
  }

  std::vector<int> vertexVisited(data->medialAxisVerticesInitial.size(), 0);
  std::vector<std::vector<int>> ccs;
  while (1) {
    int selID = -1;
    for (int vi = 0; vi < (int)vertexVisited.size(); vi++) {
      if (vertexVisited[vi] == 0) {
        selID = vi;
        break;
      }
    }

    if (selID < 0)
      break;

    std::vector<int> Q;
    Q.reserve(vertexVisited.size());
    Q.emplace_back(selID);

    size_t start = 0;
    while (start < Q.size()) {
      int curID = Q[start++];

      for (int nvi : vertexNeighboringVertices[curID]) {
        if (vertexVisited[nvi]) {
          continue;
        }

        Q.emplace_back(nvi);
        vertexVisited[nvi] = 1;
      }
    }

    ccs.emplace_back(std::move(Q));
  }

  int maxID = -1;
  int maxCC = 0;
  for (int ci = 0; ci < (int)ccs.size(); ci++) {
    if ((int)ccs[ci].size() > maxCC) {
      maxCC = (int)ccs[ci].size();
      maxID = ci;
    }
  }

  std::set<int> selectedVertexFromH, selectedVertexFromD, selectedVertexFromCC;
  // std::max((int)curvatures.size() - 300, 0)
  for (int vi = 0; vi < (int)curvatures.size(); vi++) {
    selectedVertexFromH.emplace(curvatures[vi].second);
  }

  for (int vi = 0; vi < (int)distances.size(); vi++) {
    if (distances[vi].first > 0) {
      selectedVertexFromD.emplace(distances[vi].second);
    }
  }

  for (int vi : ccs[maxID]) {
    selectedVertexFromCC.emplace(vi);
  }

  std::set<int> temp1;
  std::set_intersection(selectedVertexFromH.begin(), selectedVertexFromH.end(),
    selectedVertexFromCC.begin(), selectedVertexFromCC.end(),
    std::inserter(temp1, temp1.begin()));

  std::set_intersection(selectedVertexFromD.begin(), selectedVertexFromD.end(),
    temp1.begin(), temp1.end(),
    std::back_inserter(data->finalSelectedMAVertices));

  std::cout << "#selected vtx: " << data->finalSelectedMAVertices.size() << "/" << maVertices.size() << std::endl;
  pgo::Mesh::TriMeshGeo ss;
  // data->finalSelectedMAVertices = ccs[maxID];
  // for (int vi : data->finalSelectedMAVertices) {
  //  ss.addPos(maVertices[vi]);
  //}
  // ss.save("ss.obj");

  std::vector<std::vector<double>> vertexDistances(maVertices.size());
  tbb::parallel_for(0, (int)maVertices.size(), [&](int vi) {
    std::vector<double> &dist = vertexDistances[vi];
    dist.assign(maVertices.size(), std::numeric_limits<double>::max());
    dist[vi] = 0;

    std::priority_queue<std::pair<double, int>, std::vector<std::pair<double, int>>, std::greater<std::pair<double, int>>> Q;
    Q.emplace(0, vi);

    while (!Q.empty()) {
      int curID = Q.top().second;
      ES::V3d centerCurID = maVertices[curID];

      Q.pop();

      for (const auto &ntri : vertexNeighboringVertices[curID]) {
        if (ntri == curID)
          continue;

        ES::V3d centerNextID = maVertices[ntri];
        double weight = (centerCurID - centerNextID).norm();

        // Vec3d midPt = (centerCurID + centerNextID) * 0.5;
        // auto queryRet = targetMeshBVTree.closestTriangleQuery(targetMesh, midPt);
        double w_scaled = 1.0;  // / (std::sqrt(queryRet.dist2) + 1e-4);

        weight *= w_scaled;

        if (dist[ntri] > dist[curID] + weight) {
          dist[ntri] = dist[curID] + weight;
          Q.emplace(dist[ntri], ntri);
        }
      }
    }
  });

  int startVtx = 0;
  std::vector<int> sampleIDs;
  std::vector<ES::V3d> userInputSamplePoints;

  if (userInitFilename.length() > 0 && std::filesystem::exists(userInitFilename)) {
    try {
      nlohmann::json juser;
      std::ifstream infile(userInitFilename);
      if (infile && (infile >> juser)) {
        std::vector<std::array<double, 3>> pos = juser["vtx"].get<std::vector<std::array<double, 3>>>();
        for (const auto &p : pos) {
          userInputSamplePoints.emplace_back(p.data());
        }
      }
    }
    catch (...) {
    }
  }

  if (userInputSamplePoints.size()) {
    for (const auto &p : userInputSamplePoints) {
      for (int vi = 0; vi < (int)maVertices.size(); vi++) {
      }
    }
  }
  else {
    sampleIDs.emplace_back(startVtx);
  }

  while ((int)sampleIDs.size() < numPt) {
    double maxDist = 0;
    int sel = -1;
    for (int vi : data->finalSelectedMAVertices) {
      double minDist = 1e100;
      for (int si = 0; si < (int)sampleIDs.size(); si++) {
        double dist = vertexDistances[vi][sampleIDs[si]];
        if (dist < minDist) {
          minDist = dist;
        }
      }

      if (minDist == 1e100) {
        std::cout << vi << ',' << sampleIDs[0] << ':' << vertexDistances[vi][sampleIDs[0]] << vertexDistances[sampleIDs[0]][vi] << std::endl;
      }

      if (minDist > maxDist) {
        maxDist = minDist;
        sel = vi;
      }
    }

    sampleIDs.emplace_back(sel);
  }

  for (int si = 0; si < (int)sampleIDs.size(); si++) {
    x.segment<3>(si * 3) = maVertices[sampleIDs[si]];
  }

  data->medialAxisVerticesInitialIK = maVertices;
  data->avgLength = 0;
  int ct = 0;
  for (int ti = 0; ti < targetMesh.numTriangles(); ti++) {
    for (int j = 0; j < 3; j++) {
      double l = (targetMesh.pos(ti, j) - targetMesh.pos(ti, (j + 1) % 3)).norm();
      data->avgLength += l;
      ct += 1;
    }
  }
  data->avgLength /= ct;
  std::cout << "Avg length: " << data->avgLength << std::endl;

  data->targetVertexWeights.assign(targetMesh.numVertices(), 0);
  for (int ti = 0; ti < targetMesh.numTriangles(); ti++) {
    double a = pgo::Mesh::getTriangleArea(targetMesh.pos(ti, 0), targetMesh.pos(ti, 1), targetMesh.pos(ti, 2));
    data->targetVertexWeights[targetMesh.triVtxID(ti, 0)] += a / 3.0;
    data->targetVertexWeights[targetMesh.triVtxID(ti, 1)] += a / 3.0;
    data->targetVertexWeights[targetMesh.triVtxID(ti, 2)] += a / 3.0;
  }

  funcMinimizeCoverage(x.size(), x.data(), nullptr, data);
  // debug
  // pgo::Mesh::TriMeshGeo finalSkeleton;
  // for (int i = 0; const auto &e : data->finalSkeletonEdges) {
  //  finalSkeleton.addPos(data->finalSkeletonPoints[e.first]);
  //  finalSkeleton.addPos(data->finalSkeletonPoints[e.second]);
  //  finalSkeleton.addPos(data->finalSkeletonPoints[e.second] + pgo::asVec3d(1e-6));
  //  finalSkeleton.addTri(ES::V3i(3 * i, 3 * i + 1, 3 * i + 2));
  //  i++;
  //}
  // finalSkeleton.save("rrr0.obj");

  data->oldCenters.clear();

  for (int giter = 0; giter < nIt; giter++) {
    data->computeUnmaskedRegion = 0;
    nlopt_opt opt = nlopt_create(NLOPT_LN_NELDERMEAD, (int)x.size());
    nlopt_set_lower_bounds(opt, xlow.data());
    nlopt_set_upper_bounds(opt, xhi.data());
    nlopt_set_min_objective(opt, funcMinimizeCoverage, data);
    nlopt_set_maxeval(opt, std::min(50, (int)x.size() * 3));
    nlopt_set_xtol_abs(opt, xtol.data());
    // nlopt_set_stopval(opt, 1e-2);

    nlopt_result r;
    double minf = 0;
    if ((r = nlopt_optimize(opt, x.data(), &minf)) < 0) {
      printf("nlopt failed! Error: %d\n", r);
    }

    nlopt_destroy(opt);

    data->computeUnmaskedRegion = 1;
    funcMinimizeCoverage(x.size(), x.data(), nullptr, data);

    std::cout << "Last iteration: adding points. " << std::endl;
    std::cout << "#old pts: " << data->oldCenters.size() << std::endl;
    std::cout << "#existing pts: " << data->cellCenters.size() << std::endl;
    std::cout << "#new pts: " << data->newCenters.size() << std::endl;

    int nNewPts = (int)data->cellCenters.size() - (int)data->oldCenters.size();
    for (int i = 0; i < nNewPts; i++) {
      data->oldCenters.emplace_back(x[i * 3], x[i * 3 + 1], x[i * 3 + 2]);
    }

    x.resize(data->newCenters.size() * 3ull);
    for (int si = 0; si < (int)data->newCenters.size(); si++) {
      x.segment<3>(si * 3) = data->newCenters[si];
    }

    xtol.setConstant(x.size(), 1e-2);
    xlow.resize(x.size());
    xhi.resize(x.size());

    for (int vi = 0; vi < (int)x.size() / 3; vi++) {
      xlow.segment<3>(vi * 3) = targetMeshBB.bmin();
      xhi.segment<3>(vi * 3) = targetMeshBB.bmax();
    }

    // std::cin.get();
  }

  // funcMinimizeCoverage(x.size(), x.data(), nullptr, data);

  std::vector<EK::Point_3> vtxEK;
  for (const auto &p : data->finalSkeletonPoints) {
    vtxEK.emplace_back(p[0], p[1], p[2]);
  }

  for (auto &e : data->finalSkeletonEdges) {
    if (e.first > e.second) {
      std::swap(e.first, e.second);
    }
  }

  for (auto &t : data->finalSkeletonTriangles) {
    std::array<int, 3> tAry{ std::get<0>(t), std::get<1>(t), std::get<2>(t) };
    std::sort(tAry.begin(), tAry.end());

    t = std::tuple(tAry[0], tAry[1], tAry[2]);
  }

  SkeletonExact maSk;
  maSk.initFromInput(vtxEK, data->finalSkeletonEdges, data->finalSkeletonTriangles);
  maSk.save("opt.ma.json");

  if constexpr (1) {
    finalSkeletonPoints.clear();
    finalSkeletonEdges.clear();
    finalSkeletonTriangles.clear();

    std::map<int, int> vidMap;
    for (const auto &pr : maSk.vertices) {
      finalSkeletonPoints.emplace_back(toVec3(pr.second.pos));
      vidMap.emplace(pr.first, (int)finalSkeletonPoints.size() - 1);
    }

    for (const auto &pr : maSk.edges) {
      auto it0 = vidMap.find(pr.second.vidx[0]);
      PGO_ALOG(it0 != vidMap.end());

      auto it1 = vidMap.find(pr.second.vidx[1]);
      PGO_ALOG(it1 != vidMap.end());
      finalSkeletonEdges.emplace_back(it0->second, it1->second);
    }

    for (const auto &pr : maSk.triangles) {
      auto it0 = vidMap.find(pr.second.vidx[0]);
      PGO_ALOG(it0 != vidMap.end());

      auto it1 = vidMap.find(pr.second.vidx[1]);
      PGO_ALOG(it1 != vidMap.end());

      auto it2 = vidMap.find(pr.second.vidx[2]);
      PGO_ALOG(it2 != vidMap.end());
      finalSkeletonTriangles.emplace_back(it0->second, it1->second, it2->second);
    }
  }

  if (finalFitMeshes) {
    finalFitMeshes = new std::vector<pgo::Mesh::TriMeshGeo>(data->finalFitMeshes.size());
    for (int i = 0; i < (int)data->finalFitMeshes.size(); i++) {
      finalFitMeshes->at(i) = data->finalFitMeshes[i];
    }
  }

  delete data;
}
