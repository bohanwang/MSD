#include "cylinderCorefinement.h"
#include "corefinementUtilities.h"

#include "pgoLogging.h"
#include "libiglInterface.h"

#include "basicAlgorithms.h"
#include "createTriMesh.h"
#include "triMeshNeighbor.h"
#include "EigenSupport.h"

#include <tbb/enumerable_thread_specific.h>
#include <tbb/parallel_for.h>
#include <tbb/concurrent_hash_map.h>
#include <tbb/spin_mutex.h>

#include <chrono>
#include <unordered_map>
#include <unordered_set>

namespace CGAL
{
double to_double(const MedialAxisRepresentation::IK::FT &v)
{
  return MedialAxisRepresentation::tod(v);
}
}  // namespace CGAL

namespace MedialAxisRepresentation
{
std::tuple<K::Point_3, int> computeProjectionPt(const K::Point_3 &pt, const K::Point_3 primitiveCenter[2])
{
  K::Vector_3 vec = primitiveCenter[1] - primitiveCenter[0];
  K::Vector_3 dir = pt - primitiveCenter[0];

  // vecT(t vec + v0 - x)
  K::FT b = vec.squared_length();
  K::FT a = CGAL::scalar_product(vec, dir);
  K::FT t = a / b;
  K::Point_3 closestPt = primitiveCenter[0] + vec * t;

  if (t < 0) {
    return std::tuple(closestPt, 0);
  }
  else if (t > 1) {
    return std::tuple(closestPt, 1);
  }
  else {
    return std::tuple(closestPt, 2);
  }
};

V3 computeCylinderCoordinate(const K::Point_3 &pt, const std::tuple<K::Point_3, int> &projPt, const V3 &ax_x, const V3 &ax_y, const V3 primitiveCenter[2])
{
  V3 ptNew = toV3(pt);
  V3 origin = toV3(std::get<0>(projPt));

  real x = (ptNew - origin).dot(ax_x);
  real y = (ptNew - origin).dot(ax_y);
  real r2 = x * x + y * y;
  real r = sqrt(r2);

  real cosAngle = x / r;
  real sinAngle = y / r;
  real angle = acos(cosAngle);

  if (sinAngle < 0) {
    angle = -angle;
  }

  V3 diff0 = primitiveCenter[1] - primitiveCenter[0];
  V3 diffAbs(abs(diff0[0]), abs(diff0[1]), abs(diff0[2]));
  int maxIDX;
  diffAbs.maxCoeff(&maxIDX);

  V3 diff1 = origin - primitiveCenter[0];
  return V3(angle, diff1[maxIDX] / diff0[maxIDX], r);
};

void cylinderVertexProjection(Poly &finalMesh, const pgo::Mesh::TriMeshGeo &targetMesh, const Poly &targetMeshPoly, const FaceGraphTree &targetMeshBVTreeExact,
  const K::Point_3 primitiveCenter[2], double maxAllowedDist,
  const std::unordered_map<int, int> &vtxIDNew2Old, const std::map<int, std::array<int, 2>> &vertexIsOnTargetMesh,
  const std::vector<std::pair<int, int>> &selectedEdges)
{
  pgo::Mesh::TriMeshNeighbor targetMeshNeighbor(targetMesh);

  for (auto it = finalMesh.vertices_begin(); it != finalMesh.vertices_end(); ++it) {
    K::Point_3 pEK = it->point();
    auto projPtRet = computeProjectionPt(pEK, primitiveCenter);
    K::Point_3 projPtEK = std::get<0>(projPtRet);

    K::Vector_3 dir = pEK - projPtEK;
    double len = std::sqrt(tod(dir.squared_length()));

    // shoot a ray
    K::Ray_3 ray(projPtEK, pEK);
    const auto ret = targetMeshBVTreeExact.first_intersection(ray);

    if (ret) {
    }
    else {
      std::cerr << "No intersection happens." << std::endl;
      abort();
    }

    K::Point_3 closestPtEK = boost::get<K::Point_3>(ret->first);
    K::Vector_3 dirClosestEK = closestPtEK - projPtEK;
    double curLength = std::sqrt(tod(dirClosestEK.squared_length()));

    if (curLength <= len + maxAllowedDist) {
      it->point() = closestPtEK;
    }
    else {
      it->point() = projPtEK + dir * (len + maxAllowedDist) / len;
    }

    int newVID = it->id();
    auto vid_itt = vtxIDNew2Old.find(newVID);
    if (vid_itt != vtxIDNew2Old.end()) {
      int oldVID = vid_itt->second;
      auto vit = vertexIsOnTargetMesh.find(oldVID);
      if (vit != vertexIsOnTargetMesh.end()) {
        if (vit->second[0] == 0) {
          K::Point_3 ptgt = toP3_EK(targetMesh.pos(vit->second[1]));

          const auto &triIDs = targetMeshNeighbor.getVtxNearbyTriangles(vit->second[1]);
          K::FT maxAngle = 0;
          for (int tri : triIDs) {
            K::Plane_3 triPlane(toP3_EK(targetMesh.pos(tri, 0)), toP3_EK(targetMesh.pos(tri, 1)), toP3_EK(targetMesh.pos(tri, 2)));
            K::Point_3 origin = ptgt;
            K::Point_3 v0 = projPtEK;
            K::Point_3 v1 = triPlane.projection(v0);
            K::FT a = CGAL::approximate_angle(v0, origin, v1);
            maxAngle = CGAL::max(CGAL::abs(a), maxAngle);
          }

          // minAngle < 20
          if (maxAngle < 20) {
            std::cout << "tgt vtx: Max angle too small a=" << tod(maxAngle) << std::endl;
            if (curLength > len) {
              it->point() = CGAL::midpoint(ptgt, pEK);
            }
          }
          else {
            std::cout << "tgt vtx: Max angle ok a=" << tod(maxAngle) << std::endl;
            it->point() = ptgt;
          }

          K::Vector_3 dir = it->point() - projPtEK;
          double lenNew = std::sqrt(tod(dir.squared_length()));
          if (lenNew > curLength + 1e-7) {
            it->point() = closestPtEK;
          }
          else if (lenNew > len + maxAllowedDist) {
            it->point() = projPtEK + dir * (len + maxAllowedDist) / lenNew;
          }
        }
        else if (vit->second[0] == 1) {
          int vids[2] = { selectedEdges[vit->second[1]].first, selectedEdges[vit->second[1]].second };
          K::Point_3 pos[2] = { toP3_EK(targetMesh.pos(vids[0])), toP3_EK(targetMesh.pos(vids[1])) };
          K::Point_3 cpt;
          K::FT t;
          if constexpr (1) {
            // nT x + m = 0
            K::Vector_3 n = primitiveCenter[1] - primitiveCenter[0];
            K::FT m = -CGAL::scalar_product(n, K::Vector_3(projPtEK[0], projPtEK[1], projPtEK[2]));
            // nT (t d + p0) + m = 0
            // t (nT d) + nT p0 + m = 0
            K::Vector_3 dir = pos[1] - pos[0];
            K::FT b = CGAL::scalar_product(n, dir);
            K::FT a = CGAL::scalar_product(n, K::Vector_3(pos[0][0], pos[0][1], pos[0][2])) + m;
            t = -a / b;
            if (t < 0)
              t = 0;
            else if (t > 1)
              t = 1;

            cpt = pos[0] + dir * t;
          }
          else {
            K::Vector_3 dir = pos[1] - pos[0];
            // dirT (dir t + p0 - x) = 0
            K::Vector_3 diff = it->point() - pos[0];
            t = CGAL::scalar_product(diff, dir) / dir.squared_length();
            cpt = pos[0] + dir * t;
          }

          // double dist = std::sqrt(tod((cpt - it->point()).squared_length()));
          // std::cout << "Dist: " << dist << std::endl;
          // std::cout << "t: " << tod(t) << std::endl;

          K::FT maxAngle = 0;
          for (int ci = 0; ci < 2; ci++) {
            const auto &triIDs = targetMeshNeighbor.getVtxNearbyTriangles(vids[ci]);
            for (int tri : triIDs) {
              // if the selected edge is shared by the triangle
              const ES::V3i &triVtxIDs = targetMesh.tri(tri);
              if ((vids[0] == triVtxIDs[0] || vids[0] == triVtxIDs[1] || vids[0] == triVtxIDs[2]) &&
                (vids[1] == triVtxIDs[0] || vids[1] == triVtxIDs[1] || vids[1] == triVtxIDs[2])) {
                K::Plane_3 triPlane(toP3_EK(targetMesh.pos(tri, 0)), toP3_EK(targetMesh.pos(tri, 1)), toP3_EK(targetMesh.pos(tri, 2)));
                K::Point_3 origin = cpt;
                K::Point_3 v0 = projPtEK;
                K::Point_3 v1 = triPlane.projection(v0);
                K::FT a = CGAL::approximate_angle(v0, origin, v1);
                maxAngle = CGAL::max(CGAL::abs(a), maxAngle);
              }
            }
          }

          // minAngle < 20
          if (maxAngle < 20) {
            std::cout << "tgt vtx: Min angle too small a=" << tod(maxAngle) << std::endl;
            if (curLength > len) {
              it->point() = CGAL::midpoint(cpt, pEK);
            }
          }
          else {
            std::cout << "tgt edge: Min angle ok a=" << tod(maxAngle) << std::endl;
            it->point() = cpt;
          }

          K::Vector_3 dir = cpt - projPtEK;
          double lenNew = std::sqrt(tod(dir.squared_length()));
          if (lenNew > curLength + 1e-7) {
            it->point() = closestPtEK;
          }
          else if (lenNew > len + maxAllowedDist) {
            it->point() = projPtEK + dir * (len + maxAllowedDist) / lenNew;
          }
        }
      }
    }

    // pgo::asVec3d zz(-0.181634, 0.151783, 0.00689945);
    // pgo::asVec3d zz1 = toVec3(it->point());
    // if (std::sqrt(dot(zz1 - zz, zz1 - zz)) < 1e-4) {
    //   std::cout << "catched" << std::endl;
    //   std::cout << it->id() << std::endl;
    //   ;
    // }
  }
}
}  // namespace MedialAxisRepresentation

int MedialAxisRepresentation::corefineCylinderMeshWithTarget(const pgo::Mesh::TriMeshGeo &cylinderMesh, const ES::V3d centers[2],
  const pgo::Mesh::TriMeshGeo &targetMesh, const pgo::Mesh::TriMeshBVTree &targetMeshBVTree, const pgo::Mesh::TriMeshPseudoNormal &targetMeshNormals,
  double maxConfidentDist, double maxAllowedRadiusDist, pgo::Mesh::TriMeshGeo &meshOut, const char *filename)
{
  // double maxConfidentDist = 1e-3;
  // double maxAllowedRadiusDist = 0.1;
  constexpr int dumpMesh = 0;
  constexpr int removeSphereMeshTopology = 1;
  constexpr int triangulationOnly = 0;

  K::Point_3 primitiveCenterEK[2] = { toP3_EK(centers[0]), toP3_EK(centers[1]) };

  std::vector<int> targetMeshChoosenVertex(targetMesh.numVertices(), 0);

  HClockPt t1 = HClock::now();

  pgo::Mesh::TriMeshBVTree cylinderMeshBVTree;
  cylinderMeshBVTree.buildByInertiaPartition(cylinderMesh);

  pgo::Mesh::TriMeshPseudoNormal cylinderMeshNormal;
  cylinderMeshNormal.buildPseudoNormals(cylinderMesh);

  // compute the selected vertices for the target mesh to project
  // tbb::parallel_for(0, targetMesh.numVertices(), [&](int vi) {
  for (int vi = 0; vi < targetMesh.numVertices(); vi++) {
    // project target vertex to the skeleton
    // if the projection is successful, then keep going
    // otherwise, it is done
    ES::V3d p = targetMesh.pos(vi);
    auto projPt = computeProjectionPt(toP3_EK(p), primitiveCenterEK);
    if (std::get<1>(projPt) != 2) {
      targetMeshChoosenVertex[vi] = 0;
      continue;
      // return;
    }

    thread_local std::vector<std::tuple<double, double, int>> closestDistDueryStack;
    thread_local std::stack<int> lineSegmentQueryStack;

    // compute distance to the cylinder mesh
    // close to the sphere
    closestDistDueryStack.clear();
    auto ret = cylinderMeshBVTree.closestTriangleQuery(cylinderMesh, p, closestDistDueryStack);
    if (ret.dist2 < maxConfidentDist * maxConfidentDist) {
      targetMeshChoosenVertex[vi] = 1;
    }

    // in contact with the sphere
    ES::V3d dir = (ret.closestPosition - p).normalized();
    ES::V3d n = cylinderMeshNormal.getPseudoNormal(cylinderMesh.triangles().data(), ret.triID, ret.feature);
    if (dir.dot(n) >= 0) {
      targetMeshChoosenVertex[vi] = 1;
    }

    // compute segment from skeleton to the target position
    ES::V3d segStart = toVec3(std::get<0>(projPt));
    ES::V3d segEnd = p;

    double segLength = (segEnd - segStart).norm();
    ES::V3d segEndSelf = (segEnd - segStart) / segLength * (segLength - 1e-5) + segStart;

    // check  if it intersects with itself
    while (!lineSegmentQueryStack.empty())
      lineSegmentQueryStack.pop();

    double segW[2];
    int retID = targetMeshBVTree.lineSegmentFirstIntersectionPointSemiExact(targetMesh, segStart, segEndSelf, segW, nullptr, &lineSegmentQueryStack);
    // if it passes through itself
    if (retID >= 0) {
      targetMeshChoosenVertex[vi] = 0;
      // return;
      continue;
    }

    while (!lineSegmentQueryStack.empty())
      lineSegmentQueryStack.pop();

    // intersect with sphere mesh
    retID = cylinderMeshBVTree.lineSegmentFirstIntersectionPointSemiExact(cylinderMesh, segStart, segEnd, segW, nullptr, &lineSegmentQueryStack);

    // if there is a intersection point
    if (retID >= 0) {
      ES::V3d pt = segStart * segW[0] + segEnd * segW[1];
      double dist = (pt - p).norm();

      // if dist is small enough
      if (dist <= maxAllowedRadiusDist) {
        targetMeshChoosenVertex[vi] = 1;
      }
    }

    while (!lineSegmentQueryStack.empty())
      lineSegmentQueryStack.pop();

    segEnd = (p - segStart) * 5.0 + segStart;
    retID = cylinderMeshBVTree.lineSegmentFirstIntersectionPointSemiExact(cylinderMesh, segStart, segEnd, segW, nullptr, &lineSegmentQueryStack);

    // if there is a intersection point,
    // it means the target mesh is in contact
    if (retID >= 0) {
      ES::V3d pt = segStart * segW[0] + segEnd * segW[1];
      double dist0 = (pt - segStart).norm();
      double dist1 = (p - segStart).norm();

      if (dist1 < dist0) {
        targetMeshChoosenVertex[vi] = 1;
      }
    }
  }  // );

  if (dumpMesh) {
    pgo::Mesh::TriMeshGeo m;
    for (int i = 0; i < targetMesh.numVertices(); i++) {
      if (targetMeshChoosenVertex[i]) {
        m.addPos(targetMesh.pos(i));
      }
    }
    m.save("zra.obj");
  }

  for (int ti = 0; ti < targetMesh.numTriangles(); ti++) {
    for (int j = 0; j < 3; j++) {
      if (targetMeshChoosenVertex[targetMesh.triVtxID(ti, j)] &&
        targetMeshChoosenVertex[targetMesh.triVtxID(ti, (j + 1) % 3)] == 0) {
        ES::V3d p0 = targetMesh.pos(ti, j);
        ES::V3d p1 = targetMesh.pos(ti, (j + 1) % 3);

        auto projPt0 = computeProjectionPt(toP3_EK(p0), primitiveCenterEK);
        auto projPt1 = computeProjectionPt(toP3_EK(p1), primitiveCenterEK);
        if (std::get<1>(projPt0) == 2 && std::get<1>(projPt1) != 2) {
          targetMeshChoosenVertex[targetMesh.triVtxID(ti, (j + 1) % 3)] = 1;
        }
      }
      else if (targetMeshChoosenVertex[targetMesh.triVtxID(ti, j)] == 0 &&
        targetMeshChoosenVertex[targetMesh.triVtxID(ti, (j + 1) % 3)]) {
        ES::V3d p0 = targetMesh.pos(ti, j);
        ES::V3d p1 = targetMesh.pos(ti, (j + 1) % 3);

        auto projPt0 = computeProjectionPt(toP3_EK(p0), primitiveCenterEK);
        auto projPt1 = computeProjectionPt(toP3_EK(p1), primitiveCenterEK);
        if (std::get<1>(projPt0) != 2 && std::get<1>(projPt1) == 2) {
          targetMeshChoosenVertex[targetMesh.triVtxID(ti, j)] = 1;
        }
      }
    }
  }

  if (dumpMesh) {
    pgo::Mesh::TriMeshGeo m;
    for (int i = 0; i < targetMesh.numVertices(); i++) {
      if (targetMeshChoosenVertex[i]) {
        m.addPos(targetMesh.pos(i));
      }
    }
    m.save("zra.obj");

    ES::V3d p0 = targetMesh.pos(6968);
    ES::V3d p1 = targetMesh.pos(6573);
    auto projPt0 = computeProjectionPt(toP3_EK(p0), primitiveCenterEK);
    auto projPt1 = computeProjectionPt(toP3_EK(p1), primitiveCenterEK);

    pgo::Mesh::TriMeshGeo cc;
    cc.addMesh(pgo::Mesh::createSingleTriangleMesh(toVec3(std::get<0>(projPt0)), p0, p0 + pgo::asVec3d(1e-6)));
    cc.addMesh(pgo::Mesh::createSingleTriangleMesh(toVec3(std::get<0>(projPt1)), p1, p1 + pgo::asVec3d(1e-6)));
    cc.save("zra1.obj");
  }

  std::map<std::pair<int, int>, int> selectedEdgesMap;
  using EdgeMapIterator = std::map<std::pair<int, int>, int>::iterator;
  std::map<int, std::array<EdgeMapIterator, 3>> selectedTrianglesWithEdges;

  // compute the selected triangles for the target mesh
  // for each triangle edge, find necessary edge
  for (int ti = 0; ti < targetMesh.numTriangles(); ti++) {
    std::array<EdgeMapIterator, 3> triangleEdgeIDArray{};
    int edgeCount = 0;

    // for each triangle edge
    for (int j = 0; j < 3; j++) {
      // if the edge is selected
      if (targetMeshChoosenVertex[targetMesh.triVtxID(ti, j)] &&
        targetMeshChoosenVertex[targetMesh.triVtxID(ti, (j + 1) % 3)]) {
        // check if it has been visited
        std::pair<int, int> edgeID = { targetMesh.triVtxID(ti, j), targetMesh.triVtxID(ti, (j + 1) % 3) };
        if (edgeID.first > edgeID.second)
          std::swap(edgeID.first, edgeID.second);

        // if yes
        auto it = selectedEdgesMap.lower_bound(edgeID);
        if (it != selectedEdgesMap.end() && it->first == edgeID) {
        }
        else {
          it = selectedEdgesMap.emplace_hint(it, edgeID, 0);
        }

        // increment the counter
        triangleEdgeIDArray[edgeCount++] = it;
      }
    }

    // if a triangle is selected
    if (edgeCount == 3) {
      selectedTrianglesWithEdges.emplace(ti, triangleEdgeIDArray);
    }
  }

  // gather selected edges
  std::vector<std::pair<int, int>> selectedEdges;
  selectedEdges.reserve(selectedEdgesMap.size());
  for (auto &pr : selectedEdgesMap) {
    selectedEdges.emplace_back(pr.first);
    pr.second = (int)selectedEdges.size() - 1;
  }

  double avgEdgeLength = 0;
  for (const auto &pr : selectedEdges) {
    avgEdgeLength += (targetMesh.pos(pr.first) - targetMesh.pos(pr.second)).norm();
  }
  avgEdgeLength /= double(selectedEdges.size());

  // gather selected triangles
  std::map<std::array<int, 3>, int, ES::IntArrayLess<3>> selectedTrianglesQueryByEdgeIDs;
  for (const auto &pr : selectedTrianglesWithEdges) {
    int triID = pr.first;
    std::array<int, 3> edgeIDs = { pr.second[0]->second, pr.second[1]->second, pr.second[2]->second };
    std::sort(edgeIDs.begin(), edgeIDs.end());

    auto ret = selectedTrianglesQueryByEdgeIDs.emplace(edgeIDs, triID);
    PGO_ALOG(ret.second == true);
  }

  HClockPt t2 = HClock::now();
  SPDLOG_LOGGER_INFO(pgo::Logging::lgr(), "S0 done. t: {:.4f}s", duraSecond(t1, t2));

  t1 = HClock::now();

  // create target mesh
  // construct target mesh BV tree
  std::vector<K::Point_3> targetMeshVerticesER(targetMesh.numVertices());
  tbb::parallel_for(
    0, targetMesh.numVertices(), [&](int vi) {
      const ES::V3d &v_d = targetMesh.pos(vi);
      targetMeshVerticesER[vi] = K::Point_3(v_d[0], v_d[1], v_d[2]);
    },
    tbb::static_partitioner());

  Poly targetMeshPoly;
  pgo::CGALInterface::triangleMesh2Polyhedron(targetMesh, targetMeshPoly);
  FaceGraphTree targetMeshBVTreeExact(CGAL::faces(targetMeshPoly).first, CGAL::faces(targetMeshPoly).second, targetMeshPoly);

  // compute axis
  V3 ax_x = V3::UnitX();
  V3 diff = toV3(primitiveCenterEK[1]) - toV3(primitiveCenterEK[0]);
  real diff_len = diff.squaredNorm();
  V3 ax_z = diff / sqrt(diff_len);

  if (abs(ax_x.dot(ax_z)) < 1e-6) {
    ax_x = V3::UnitY();
  }

  V3 ax_y = ax_z.cross(ax_x);
  ax_y /= sqrt(ax_y.squaredNorm());

  ax_x = ax_y.cross(ax_z);
  ax_x /= sqrt(ax_x.squaredNorm());

  std::vector<K::Point_3> cylinderMeshVerticesER(cylinderMesh.numVertices());
  for (int vi = 0; vi < cylinderMesh.numVertices(); vi++) {
    cylinderMeshVerticesER[vi] = K::Point_3(cylinderMesh.pos(vi)[0], cylinderMesh.pos(vi)[1], cylinderMesh.pos(vi)[2]);
  }

  Poly cylinderMeshPoly;
  pgo::CGALInterface::triangleMesh2Polyhedron(cylinderMesh, cylinderMeshPoly);
  FaceGraphTree cylinderMeshBVTreeExact(CGAL::faces(cylinderMeshPoly).first, CGAL::faces(cylinderMeshPoly).second, cylinderMeshPoly);

  struct CylinderTriangleInfo
  {
    std::array<K::Triangle_2, 2> triangle2D;
    std::array<int, 2> isProjectionValid;
    std::array<real, 3> vertexRadii;
    int projID;
  };

  std::vector<CylinderTriangleInfo> cylinderMeshProjectedTriangles(cylinderMesh.numTriangles());
  K::FT K_pi(M_PI);
  K::FT K_pi_half = K_pi / 2;
  K::FT K_pi_threshold = K_pi * 0.99999;

  auto rotateAngle = [&K_pi](const K::FT &val) -> K::FT {
    if (val < 0) {
      return K_pi + val;
    }
    else {
      return val - K_pi;
    }
  };

  auto pt2DEqual = [&K_pi, &rotateAngle](const K::Point_2 &p0, const K::Point_2 &p1) -> bool {
    bool ret1 = p0[1] == p1[1];
    bool ret0 = (p0[0] == p1[0] || rotateAngle(p0[0]) == p1[0] || rotateAngle(p1[0]) == p0[0]);
    /*
    std::cout << "Check equal:\n";
    std::cout << "  " << tod(p0[1]) << ',' << tod(p1[1]) << std::endl;
    std::cout << "  " << tod(p0[0]) << ',' << tod(p1[0]) << std::endl;
    std::cout << "  " << tod(rotateAngle(p0[0])) << ',' << tod(p1[0]) << std::endl;
    std::cout << "  " << tod(p0[0]) << ',' << tod(rotateAngle(p1[0])) << std::endl;
    std::cout << "  " << ret0 << ',' << ret1 << std::endl;
    */

    return ret0 && ret1;
  };

  K::Vector_3 axisDiff_EK = primitiveCenterEK[1] - primitiveCenterEK[0];
  V3 axisDiff(EK_to_IK<real>(axisDiff_EK[0]), EK_to_IK<real>(axisDiff_EK[1]), EK_to_IK<real>(axisDiff_EK[2]));
  V3 primitiveCenterIK[2] = {
    toV3(primitiveCenterEK[0]),
    toV3(primitiveCenterEK[1])
  };

  auto pt2D_to_pt3D = [&](const K::Point_2 &pt, int triID) -> IK::Point_3 {
    int projID = cylinderMeshProjectedTriangles[triID].projID;
    const auto &tri = cylinderMeshProjectedTriangles[triID].triangle2D[projID];
    K::FT w[3];
    CGAL::Barycentric_coordinates::triangle_coordinates_2(tri[0], tri[1], tri[2], pt, w);

    /*
    std::cout << "angle: " << tod(pt[0]) << ',' << tod(pt[1]) << std::endl;
    std::cout << tod(w[0]) << ',' << tod(w[1]) << ',' << tod(w[2]) << std::endl;
    */

    IK::FT r = 0;
    for (int j = 0; j < 3; j++) {
      r += cylinderMeshProjectedTriangles[triID].vertexRadii[j] * EK_to_IK<real>(w[j]);
    }

    V3 x;
    if (projID == 0) {
      real angle = EK_to_IK<real>(pt[0]);
      x = (ax_x * cos(angle) + ax_y * sin(angle)) * r;
    }
    else {
      real angle = EK_to_IK<real>(pt[0]);
      x = -(ax_x * cos(angle) + ax_y * sin(angle)) * r;
    }

    V3 vert = axisDiff * EK_to_IK<real>(pt[1]) + primitiveCenterIK[0];
    x += vert;

    return IK::Point_3(x[0], x[1], x[2]);
  };

  for (int ti = 0; ti < cylinderMesh.numTriangles(); ti++) {
    // tbb::parallel_for(0, cylinderMesh.numTriangles(), [&](int ti) {
    K::Point_2 p2d[2][3];
    for (int j = 0; j < 3; j++) {
      K::Point_3 pt = toP3_EK(cylinderMesh.pos(ti, j));
      auto projPt = computeProjectionPt(pt, primitiveCenterEK);
      V3 uvw = computeCylinderCoordinate(pt, projPt, ax_x, ax_y, primitiveCenterIK);

      p2d[0][j] = K::Point_2(IK_to_EK(uvw[0]), IK_to_EK(uvw[1]));
      p2d[1][j] = K::Point_2(rotateAngle(p2d[0][j][0]), p2d[0][j][1]);

      cylinderMeshProjectedTriangles[ti].vertexRadii[j] = uvw[2];
      // V2 uv1 = computeCylinderCoordinate(toV3(cylinderMesh.pos(ti, j)), -ax_x, -ax_y);
      // fmt::print("{},{},{}\n", uv[0], uv1[0], tod(p2d[1][j][0]));
      // std::cout << std::flush;
    }

    cylinderMeshProjectedTriangles[ti].triangle2D[0] = K::Triangle_2(p2d[0][0], p2d[0][1], p2d[0][2]);
    cylinderMeshProjectedTriangles[ti].triangle2D[1] = K::Triangle_2(p2d[1][0], p2d[1][1], p2d[1][2]);

    int count = 0;
    for (int j = 0; j < 3; j++) {
      if (CGAL::abs(p2d[0][j][0]) < K_pi_half) {
        count++;
      }
    }

    int selID = 0;
    if (count > 3 - count) {
      selID = 0;
    }
    else {
      selID = 1;
    }

    K::FT diff[2][3] = {
      { CGAL::abs(p2d[0][0][0] - p2d[0][1][0]),
        CGAL::abs(p2d[0][1][0] - p2d[0][2][0]),
        CGAL::abs(p2d[0][2][0] - p2d[0][0][0]) },
      { CGAL::abs(p2d[1][0][0] - p2d[1][1][0]),
        CGAL::abs(p2d[1][1][0] - p2d[1][2][0]),
        CGAL::abs(p2d[1][2][0] - p2d[1][0][0]) }
    };

    // std::cout << ti << '\n';
    // fmt::print("{},{},{}\n", tod(p2d[selID][0][0]), tod(p2d[selID][1][0]), tod(p2d[selID][2][0]));
    // fmt::print("{},{},{}\n", tod(diff[0]), tod(diff[1]), tod(diff[2]));
    // std::cout << std::flush;

    PGO_ALOG(diff[selID][0] < K_pi_threshold);
    PGO_ALOG(diff[selID][1] < K_pi_threshold);
    PGO_ALOG(diff[selID][2] < K_pi_threshold);

    for (int ci = 0; ci < 2; ci++) {
      if (diff[ci][0] < K_pi_threshold && diff[ci][1] < K_pi_threshold && diff[ci][2] < K_pi_threshold) {
        cylinderMeshProjectedTriangles[ti].isProjectionValid[ci] = 1;
      }
      else {
        cylinderMeshProjectedTriangles[ti].isProjectionValid[ci] = 0;
      }
    }

    cylinderMeshProjectedTriangles[ti].projID = selID;
  }  // );

  // pgo::Mesh::TriMeshGeo mm;
  // for (int ti = 0; ti < (int)cylinderMeshProjectedTriangles.size(); ti++) {
  //   int projID = cylinderMeshProjectedTriangles[ti].projID;
  //   IK::Point_3 p[3];
  //   for (int j = 0; j < 3; j++) {
  //     p[j] = pt2D_to_pt3D(cylinderMeshProjectedTriangles[ti].triangle2D[projID][j], ti);
  //   }

  //   ES::V3d pV3d[3] = {
  //     ES::V3d(tod(p[0][0]), tod(p[0][1]), tod(p[0][2])),
  //     ES::V3d(tod(p[1][0]), tod(p[1][1]), tod(p[1][2])),
  //     ES::V3d(tod(p[2][0]), tod(p[2][1]), tod(p[2][2])),
  //   };

  //   mm.addMesh(pgo::Mesh::createSingleTriangleMesh(pV3d[0], pV3d[1], pV3d[2]));
  // }
  // mm.save("rarrr.obj");

  // build cutting mesh
  std::vector<K::Segment_2> cuttingSegments[2];

  cuttingSegments[0].resize(selectedEdges.size());
  cuttingSegments[1].resize(selectedEdges.size());

  std::vector<real> targetMeshVertexRadius(targetMesh.numVertices());
  std::vector<tbb::spin_mutex> targetMeshVertexLocks(targetMesh.numVertices());

  // for each edge find all intersected triangles
  std::vector<std::vector<std::tuple<int, std::array<int, 2>>>> edgeIntersectedTriangles(selectedEdges.size());
  tbb::parallel_for(0, (int)selectedEdges.size(), [&](int ei) {
    // for (int ei = 0; ei < (int)selectedEdges.size(); ei++) {
    int vids[2] = {
      selectedEdges[ei].first,
      selectedEdges[ei].second
    };

    K::Segment_2 seg[2];

    K::Point_3 pt0 = toP3_EK(targetMesh.pos(vids[0]));
    auto projPt0 = computeProjectionPt(pt0, primitiveCenterEK);
    V3 coord0 = computeCylinderCoordinate(pt0, projPt0, ax_x, ax_y, primitiveCenterIK);

    K::Point_3 pt1 = toP3_EK(targetMesh.pos(vids[1]));
    auto projPt1 = computeProjectionPt(pt1, primitiveCenterEK);
    V3 coord1 = computeCylinderCoordinate(pt1, projPt1, ax_x, ax_y, primitiveCenterIK);

    K::Vector_3 dir0 = pt0 - std::get<0>(projPt0);
    double dir0_length = std::sqrt(tod(dir0.squared_length()));

    K::Vector_3 dir1 = pt1 - std::get<0>(projPt1);
    double dir1_length = std::sqrt(tod(dir1.squared_length()));

    K::Point_3 pt0plus = std::get<0>(projPt0) + dir0 * std::max(0.2 / dir0_length, 2.0);
    K::Point_3 pt1plus = std::get<0>(projPt1) + dir1 * std::max(0.2 / dir1_length, 2.0);

    K::Point_3 pt01plus = std::get<0>(projPt1) + dir0 * std::max(0.2 / dir0_length, 2.0);
    K::Point_3 pt10plus = std::get<0>(projPt0) + dir1 * std::max(0.2 / dir1_length, 2.0);

    K::Point_3 pts[6] = {
      pt0plus,
      pt1plus,

      std::get<0>(projPt0),
      std::get<0>(projPt1),

      pt01plus,
      pt10plus,
    };

    // CGAL::oriented_bounding_box(sm, obb_points,  CGAL::parameters::use_convex_hull(true));
    K::Iso_cuboid_3 bb(CGAL::bbox_3(pts, pts + 6));
    // K::Point_3 mid = CGAL::midpoint(bb.min(), bb.max());
    // K::Vector_3 diff = bb.max() - mid;
    // K::Point_3 expanded_min = mid - diff * 1.1;
    // K::Point_3 expanded_max = mid + diff * 1.1;
    // K::Iso_cuboid_3 expanded_bb(expanded_min, expanded_max);

    thread_local std::vector<Poly::Facet_handle> faceHandles;

    faceHandles.clear();
    cylinderMeshBVTreeExact.all_intersected_primitives(bb, std::back_inserter(faceHandles));

    K::Point_2 uv0 = K::Point_2(IK_to_EK(coord0[0]), IK_to_EK(coord0[1]));
    K::Point_2 uv1 = K::Point_2(IK_to_EK(coord1[0]), IK_to_EK(coord1[1]));
    seg[0] = K::Segment_2(uv0, uv1);

    targetMeshVertexLocks[vids[0]].lock();
    targetMeshVertexRadius[vids[0]] = coord0[2];
    targetMeshVertexLocks[vids[0]].unlock();

    targetMeshVertexLocks[vids[1]].lock();
    targetMeshVertexRadius[vids[1]] = coord1[2];
    targetMeshVertexLocks[vids[1]].unlock();

    seg[1] = K::Segment_2(K::Point_2(rotateAngle(uv0[0]), uv0[1]), K::Point_2(rotateAngle(uv1[0]), uv1[1]));

    cuttingSegments[0][ei] = seg[0];
    cuttingSegments[1][ei] = seg[1];

    K::FT diff[2] = {
      CGAL::abs(seg[0][0][0] - seg[0][1][0]),
      CGAL::abs(seg[1][0][0] - seg[1][1][0]),
    };

    bool valid[2] = {
      diff[0] < K_pi_threshold,
      diff[1] < K_pi_threshold
    };

    if (valid[0] || valid[1]) {
    }
    else {
      std::cout << "s0: " << tod(seg[0][0][0]) << ',' << tod(seg[0][1][0]) << std::endl;
      std::cout << "s1: " << tod(seg[1][0][0]) << ',' << tod(seg[1][1][0]) << std::endl;
      std::cout << "diff: " << tod(diff[0]) << ',' << tod(diff[1]) << std::endl;
      abort();
    }

    for (auto fh : faceHandles) {
      int ti = (int)fh->id();
      int projID = cylinderMeshProjectedTriangles[ti].projID;
      std::array<int, 2> isValid{ 0, 0 };
      for (int ci = 0; ci < 2; ci++) {
        if (cylinderMeshProjectedTriangles[ti].isProjectionValid[ci] && valid[ci]) {
          isValid[ci] = 1;
          projID = ci;
        }
      }

      if (isValid[0] || isValid[1]) {
        bool ret = CGAL::do_intersect(cylinderMeshProjectedTriangles[ti].triangle2D[projID], seg[projID]);
        if (ret) {
          edgeIntersectedTriangles[ei].emplace_back(std::tuple(ti, isValid));
        }
      }
      // else {
      //   std::cout << "t0: " << tod(cylinderMeshProjectedTriangles[ti].triangle2D[0][0][0]) << ','
      //             << tod(cylinderMeshProjectedTriangles[ti].triangle2D[0][1][0]) << ','
      //             << tod(cylinderMeshProjectedTriangles[ti].triangle2D[0][2][0]) << std::endl;

      //   std::cout << "t1: " << tod(cylinderMeshProjectedTriangles[ti].triangle2D[1][0][0]) << ','
      //             << tod(cylinderMeshProjectedTriangles[ti].triangle2D[1][1][0]) << ','
      //             << tod(cylinderMeshProjectedTriangles[ti].triangle2D[1][2][0]) << std::endl;
      // }

      // if (ei == 520 && (ti == 324 || ti == 325)) {
      //   std::cout << ti << ',' << projID << std::endl;
      //   pgo::asVec3d ps[3];
      //   for (int j = 0; j < 3; j++) {
      //     ps[j][0] = tod(cylinderMeshProjectedTriangles[ti].triangle2D[projID][j][0]);
      //     ps[j][1] = tod(cylinderMeshProjectedTriangles[ti].triangle2D[projID][j][1]);
      //     ps[j][2] = 0.0;
      //   }

      //   pgo::Mesh::TriMeshGeo m;
      //   m.addMesh(pgo::Mesh::.createSingleTriangleMesh(ps[0], ps[1], ps[2]));
      //   m.save("mss.obj");

      //   pgo::asVec3d ss[2];
      //   for (int j = 0; j < 2; j++) {
      //     ss[j][0] = tod(seg[projID][j][0]);
      //     ss[j][1] = tod(seg[projID][j][1]);
      //     ss[j][2] = 0.0;
      //   }

      //   std::cout << ss[0] << ',' << ss[1] << std::endl;

      //   m.addMesh(pgo::Mesh::.createSingleTriangleMesh(ss[0], ss[1], ss[0] + pgo::asVec3d(1e-6)));
      //   m.save("mss.obj");
      // }
    }
  });

  // for each intersected triangle, gather intersections
  std::vector<std::vector<int>> triangleIntersectedEdges(cylinderMesh.numTriangles());
  std::vector<std::array<int, 2>> triangleValidProjection(cylinderMesh.numTriangles(), std::array<int, 2>{ 1, 1 });
  for (int ei = 0; ei < (int)selectedEdges.size(); ei++) {
    for (int i = 0; i < (int)edgeIntersectedTriangles[ei].size(); i++) {
      const auto &tp = edgeIntersectedTriangles[ei][i];
      int ti = std::get<0>(tp);

      triangleValidProjection[ti][0] = triangleValidProjection[ti][0] & std::get<1>(tp)[0];
      triangleValidProjection[ti][1] = triangleValidProjection[ti][1] & std::get<1>(tp)[1];
      triangleIntersectedEdges[ti].emplace_back(ei);
    }
  }

  for (int ti = 0; ti < cylinderMesh.numTriangles(); ti++) {
    if (triangleValidProjection[ti][0]) {
      cylinderMeshProjectedTriangles[ti].projID = 0;
    }
    else if (triangleValidProjection[ti][1]) {
      cylinderMeshProjectedTriangles[ti].projID = 1;
    }
    else {
      std::cout << ti << std::endl;
      throw std::runtime_error("Invalid projection");
    }
  }

  [[maybe_unused]] auto dumpTri = [](const K::Triangle_2 &tri) -> pgo::Mesh::TriMeshGeo {
    ES::V3d p[3];
    for (int j = 0; j < 3; j++) {
      p[j][0] = tod(tri[j][0]);
      p[j][1] = tod(tri[j][1]);
      p[j][2] = 0;
    }

    return pgo::Mesh::createSingleTriangleMesh(p[0], p[1], p[2]);
  };

  [[maybe_unused]] auto dumpSeg = [](const K::Segment_2 &seg) -> pgo::Mesh::TriMeshGeo {
    ES::V3d p[2];
    for (int j = 0; j < 2; j++) {
      p[j][0] = tod(seg[j][0]);
      p[j][1] = tod(seg[j][1]);
      p[j][2] = 0;
    }

    return pgo::Mesh::createSingleTriangleMesh(p[0], p[1], p[0] + pgo::asVec3d(1e-4));
  };

  // =====================================
  // for each intersected triangle, compute corefined triangles
  // we use cutting mesh to cut each sphere triangle
  std::vector<Patch2D> corefinedPatches(cylinderMesh.numTriangles());
  tbb::parallel_for(0, cylinderMesh.numTriangles(), [&](int ti) {
    int projID = cylinderMeshProjectedTriangles[ti].projID;
    // for (int ti = 0; ti < cylinderMesh.numTriangles(); ti++) {
    const K::Point_2 &uv_v0 = cylinderMeshProjectedTriangles[ti].triangle2D[projID][0];
    const K::Point_2 &uv_v1 = cylinderMeshProjectedTriangles[ti].triangle2D[projID][1];
    const K::Point_2 &uv_v2 = cylinderMeshProjectedTriangles[ti].triangle2D[projID][2];

    corefinedPatches[ti].rawTriangles.emplace_back(uv_v0, uv_v1, uv_v2);

    if (triangleIntersectedEdges[ti].size() == 0ull) {
      corefinedPatches[ti].vertices = { uv_v0, uv_v1, uv_v2 };
      corefinedPatches[ti].triangles.emplace_back(std::array<int, 3>{ 0, 1, 2 });

      corefinedPatches[ti].edgeVertices[0] = { 0, 1 };
      corefinedPatches[ti].edgeVertices[1] = { 1, 2 };
      corefinedPatches[ti].edgeVertices[2] = { 2, 0 };

      corefinedPatches[ti].cornerVertices[0] = 0;
      corefinedPatches[ti].cornerVertices[1] = 1;
      corefinedPatches[ti].cornerVertices[2] = 2;

      corefinedPatches[ti].vertexIsOnTargetMesh.assign(corefinedPatches[ti].vertices.size(), std::array<int, 2>{ -1, -1 });
    }
    else {
      // bfs
      for (int ei = 0; ei < (int)triangleIntersectedEdges[ti].size(); ei++) {
        int edgeID = triangleIntersectedEdges[ti][ei];
        const K::Segment_2 &seg = cuttingSegments[projID][edgeID];

        auto iter = corefinedPatches[ti].rawTriangles.begin();
        size_t count = corefinedPatches[ti].rawTriangles.size();

        for (size_t ci = 0; ci < count; ci++) {
          const auto result = CGAL::intersection(seg, *iter);
          if (result) {
            if (const K::Segment_2 *s = boost::get<K::Segment_2>(&*result)) {
              // create 3 triangles for the first point,
              // the segment must be in one of the triangle
              K::Triangle_2 remainingTri;
              bool found = false;
              for (int i = 0; i < 3; i++) {
                K::Triangle_2 tri((*s)[0], (*iter)[i], (*iter)[(i + 1) % 3]);
                if (tri.is_degenerate()) {
                  continue;
                }
                if (!tri.has_on_unbounded_side((*s)[1]) && found == false) {
                  remainingTri = tri;
                  found = true;
                }
                else {
                  corefinedPatches[ti].rawTriangles.emplace_back(tri);
                }
              }

              // pgo::Mesh::TriMeshGeo m1 = dumpTri(*iter);
              // pgo::Mesh::TriMeshGeo m2 = dumpSeg(*s);

              // pgo::Mesh::TriMeshGeo mall;
              // mall.addMesh(m1);
              // mall.addMesh(m2);
              // mall.save("rzzz.obj");

              PGO_ALOG(CGAL::abs(remainingTri.area()) > 0);

              // create 3 triangles for the second point in the remaining triangle
              // the segment must be on the edge of some triangles
              for (int i = 0; i < 3; i++) {
                K::Triangle_2 tri((*s)[1], remainingTri[i], remainingTri[(i + 1) % 3]);
                if (!tri.is_degenerate()) {
                  corefinedPatches[ti].rawTriangles.emplace_back(tri);
                }
              }
            }
            else if (const K::Point_2 *p = boost::get<K::Point_2>(&*result)) {
              // create 3 triangles
              for (int i = 0; i < 3; i++) {
                K::Triangle_2 tri(*p, (*iter)[i], (*iter)[(i + 1) % 3]);

                if (!tri.is_degenerate()) {
                  corefinedPatches[ti].rawTriangles.emplace_back(tri);
                }
              }
            }
            else {
              throw std::runtime_error("Impossible");
            }

            corefinedPatches[ti].rawTriangles.pop_front();
            iter = corefinedPatches[ti].rawTriangles.begin();
          }
          else {
            corefinedPatches[ti].rawTriangles.splice(corefinedPatches[ti].rawTriangles.end(), corefinedPatches[ti].rawTriangles, corefinedPatches[ti].rawTriangles.begin());
            iter = corefinedPatches[ti].rawTriangles.begin();
          }
        }  // end result
      }

      corefinedPatches[ti].vertices.reserve(corefinedPatches[ti].rawTriangles.size() * 3);
      corefinedPatches[ti].triangles.reserve(corefinedPatches[ti].rawTriangles.size());
      for (const auto &tri : corefinedPatches[ti].rawTriangles) {
        std::array<int, 3> triVtx = { -1, -1, -1 };
        for (int j = 0; j < 3; j++) {
          // find if the vertex exist
          for (int vi = 0; vi < (int)corefinedPatches[ti].vertices.size(); vi++) {
            const auto &vtx = corefinedPatches[ti].vertices[vi];
            if (tri[j] == vtx) {
              triVtx[j] = vi;
              break;
            }
          }

          // if not found
          if (triVtx[j] < 0) {
            corefinedPatches[ti].vertices.emplace_back(tri[j]);
            triVtx[j] = (int)corefinedPatches[ti].vertices.size() - 1;
          }
        }

        corefinedPatches[ti].triangles.emplace_back(triVtx);
      }

      corefinedPatches[ti].edgeVertices[0].reserve(corefinedPatches[ti].vertices.size());
      corefinedPatches[ti].edgeVertices[1].reserve(corefinedPatches[ti].vertices.size());
      corefinedPatches[ti].edgeVertices[2].reserve(corefinedPatches[ti].vertices.size());

      K::Segment_2 e01(uv_v0, uv_v1), e12(uv_v1, uv_v2), e20(uv_v2, uv_v0);
      for (int vi = 0; vi < (int)corefinedPatches[ti].vertices.size(); vi++) {
        if (e01.has_on(corefinedPatches[ti].vertices[vi])) {
          corefinedPatches[ti].edgeVertices[0].emplace_back(vi);
        }

        if (e12.has_on(corefinedPatches[ti].vertices[vi])) {
          corefinedPatches[ti].edgeVertices[1].emplace_back(vi);
        }

        if (e20.has_on(corefinedPatches[ti].vertices[vi])) {
          corefinedPatches[ti].edgeVertices[2].emplace_back(vi);
        }

        if (corefinedPatches[ti].vertices[vi] == uv_v0) {
          corefinedPatches[ti].cornerVertices[0] = vi;
        }
        else if (corefinedPatches[ti].vertices[vi] == uv_v1) {
          corefinedPatches[ti].cornerVertices[1] = vi;
        }
        else if (corefinedPatches[ti].vertices[vi] == uv_v2) {
          corefinedPatches[ti].cornerVertices[2] = vi;
        }
      }

      // record the target edge ID for each newly generated triangle edges
      corefinedPatches[ti].triangleEdgeIDs.reserve((int)corefinedPatches[ti].triangles.size());
      for (int tii = 0; tii < (int)corefinedPatches[ti].triangles.size(); tii++) {
        for (int j = 0; j < 3; j++) {
          const K::Point_2 &p0 = corefinedPatches[ti].vertices[corefinedPatches[ti].triangles[tii][j]];
          const K::Point_2 &p1 = corefinedPatches[ti].vertices[corefinedPatches[ti].triangles[tii][(j + 1) % 3]];

          for (int ei = 0; ei < (int)triangleIntersectedEdges[ti].size(); ei++) {
            int edgeID = triangleIntersectedEdges[ti][ei];

            if (cuttingSegments[projID][edgeID].has_on(p0) &&
              cuttingSegments[projID][edgeID].has_on(p1)) {
              corefinedPatches[ti].triangleEdgeIDs.emplace_back(tii, j, triangleIntersectedEdges[ti][ei]);
              break;
            }
          }
        }
      }

      corefinedPatches[ti].vertexIsOnTargetMesh.assign(corefinedPatches[ti].vertices.size(), std::array<int, 2>{ -1, -1 });
      for (int vi = 0; vi < (int)corefinedPatches[ti].vertices.size(); vi++) {
        for (int ei = 0; ei < (int)triangleIntersectedEdges[ti].size(); ei++) {
          int edgeID = triangleIntersectedEdges[ti][ei];
          const K::Segment_2 &seg = cuttingSegments[projID][edgeID];

          if (seg[0] == corefinedPatches[ti].vertices[vi]) {
            corefinedPatches[ti].vertexIsOnTargetMesh[vi][0] = 0;
            corefinedPatches[ti].vertexIsOnTargetMesh[vi][1] = selectedEdges[edgeID].first;
          }
          else if (seg[1] == corefinedPatches[ti].vertices[vi]) {
            corefinedPatches[ti].vertexIsOnTargetMesh[vi][0] = 0;
            corefinedPatches[ti].vertexIsOnTargetMesh[vi][1] = selectedEdges[edgeID].second;
          }
          else if (seg.has_on(corefinedPatches[ti].vertices[vi])) {
            corefinedPatches[ti].vertexIsOnTargetMesh[vi][0] = 1;
            corefinedPatches[ti].vertexIsOnTargetMesh[vi][1] = edgeID;
          }
        }
      }
    }

    if (ti % 100 == 0)
      std::cout << ti << "  " << std::flush;
  });

  std::cout << std::endl;

  t2 = HClock::now();
  SPDLOG_LOGGER_INFO(pgo::Logging::lgr(), "S1 Done. Time cost: {:.4f}s", duraSecond(t1, t2));

  // ================================================================
  // next we merge all cut triangles into a closed manifold mesh
  pgo::Mesh::TriMeshNeighbor cylinderMeshNeighbor(cylinderMesh);

  std::vector<int> triIDQ;
  triIDQ.reserve(cylinderMesh.numTriangles() * 2);

  std::vector<int> isTriVisited(cylinderMesh.numTriangles(), 0);

  std::vector<K::Point_2> finalPoints2D;
  finalPoints2D.reserve(cylinderMesh.numTriangles() * 10);

  std::vector<IK::Point_3> finalPoints3D;
  finalPoints3D.reserve(cylinderMesh.numTriangles() * 10);

  std::vector<ES::V3i> finalTriangles;
  finalTriangles.reserve(cylinderMesh.numTriangles() * 10);

  std::map<std::pair<int, int>, std::vector<int>> edgeVertexIDs;
  std::vector<int> newVertexIDs(cylinderMesh.numVertices(), -1);

  triIDQ.emplace_back(0);
  isTriVisited[0] = 1;

  size_t start = 0;
  std::vector<int> vertexMasks;
  vertexMasks.reserve(10000);

  std::vector<int> triangleNewIDs;
  triangleNewIDs.reserve(10000);

  std::vector<std::vector<std::array<int, 2>>> targetEdgeOverlappingTriangleIDs(selectedEdges.size());
  std::map<int, std::array<int, 2>> vertexIsOnTargetMesh;

  while (triIDQ.size() > start) {
    int triIDCur = triIDQ[start++];

    // this one store the vertex global IDs
    // if it is -1, it means it is not addressed or found
    vertexMasks.assign(corefinedPatches[triIDCur].vertices.size(), -1);
    int edgeMask[3] = { 0, 0, 0 };

    for (int j = 0; j < 3; j++) {
      std::pair<int, int> edgeID{ cylinderMesh.triVtxID(triIDCur, j), cylinderMesh.triVtxID(triIDCur, (j + 1) % 3) };
      if (edgeID.first > edgeID.second)
        std::swap(edgeID.first, edgeID.second);

      // first check if this global sphere edge has been visited
      auto eitt = edgeVertexIDs.find(edgeID);
      // if no
      if (eitt == edgeVertexIDs.end()) {
        edgeMask[j] = 1;
      }
      // if yes, reuse existing vertex id
      else {
        if (eitt->second.size() != corefinedPatches[triIDCur].edgeVertices[j].size()) {
          std::vector<IK::Point_3> vtx;
          for (int i = 0; i < (int)corefinedPatches[triIDCur].vertices.size(); i++) {
            vtx.emplace_back(pt2D_to_pt3D(corefinedPatches[triIDCur].vertices[i], triIDCur));
          }

          pgo::CGALInterface::Polyhedron<IK> finalMesh;
          pgo::CGALInterface::PolyhedronBuilderNonTriangle<IK, std::vector<IK::Point_3>::iterator, std::vector<std::array<int, 3>>::iterator> builder(
            vtx.begin(), corefinedPatches[triIDCur].triangles.begin(), (int)vtx.size(), (int)corefinedPatches[triIDCur].triangles.size());
          finalMesh.delegate(builder);

          std::ofstream("rzr1.off") << finalMesh;

          if (1) {
            pgo::CGALInterface::Polyhedron<IK> finalMesh;
            pgo::CGALInterface::PolyhedronBuilderNonTriangle<IK, std::vector<IK::Point_3>::iterator, std::vector<ES::V3i>::iterator> builder(
              finalPoints3D.begin(), finalTriangles.begin(), (int)finalPoints3D.size(), (int)finalTriangles.size());
            finalMesh.delegate(builder);

            std::ofstream("rzr0.off") << finalMesh;
          }

          abort();
        }

        // find global vertex ID
        for (int k = 0; k < (int)corefinedPatches[triIDCur].edgeVertices[j].size(); k++) {
          if (vertexMasks[corefinedPatches[triIDCur].edgeVertices[j][k]] >= 0)
            continue;

          bool found = false;
          for (int r = 0; r < (int)eitt->second.size(); r++) {
            if (pt2DEqual(corefinedPatches[triIDCur].vertices[corefinedPatches[triIDCur].edgeVertices[j][k]], finalPoints2D[eitt->second[r]])) {
              vertexMasks[corefinedPatches[triIDCur].edgeVertices[j][k]] = eitt->second[r];
              found = true;
              break;
            }
          }

          if (found != true) {
            std::vector<IK::Point_3> vtx;
            for (int i = 0; i < (int)corefinedPatches[triIDCur].vertices.size(); i++) {
              vtx.emplace_back(pt2D_to_pt3D(corefinedPatches[triIDCur].vertices[i], triIDCur));
            }

            pgo::CGALInterface::Polyhedron<IK> finalMesh;
            pgo::CGALInterface::PolyhedronBuilderNonTriangle<IK, std::vector<IK::Point_3>::iterator, std::vector<std::array<int, 3>>::iterator> builder(
              vtx.begin(), corefinedPatches[triIDCur].triangles.begin(), (int)vtx.size(), (int)corefinedPatches[triIDCur].triangles.size());
            finalMesh.delegate(builder);
            // std::ofstream("rzr1.off") << finalMesh;

            abort();
          }
        }
      }
    }  // end for j \in [0, 1, 2]

    // address original triangle corners
    for (int j = 0; j < 3; j++) {
      int vi = corefinedPatches[triIDCur].cornerVertices[j];
      if (newVertexIDs[cylinderMesh.triVtxID(triIDCur, j)] >= 0) {
        if (vertexMasks[vi] >= 0) {
          PGO_ALOG(vertexMasks[vi] == newVertexIDs[cylinderMesh.triVtxID(triIDCur, j)]);
        }
        else {
          vertexMasks[vi] = newVertexIDs[cylinderMesh.triVtxID(triIDCur, j)];
        }
      }
    }

    // compute all vertex ID,
    // especially for the ones that has not been visited
    for (int j = 0; j < (int)corefinedPatches[triIDCur].vertices.size(); j++) {
      if (vertexMasks[j] < 0) {
        finalPoints2D.emplace_back(corefinedPatches[triIDCur].vertices[j]);
        finalPoints3D.emplace_back(pt2D_to_pt3D(corefinedPatches[triIDCur].vertices[j], triIDCur));

        // corefinedPatches[triIDCur].vertices[j] to 3D
        vertexMasks[j] = (int)finalPoints2D.size() - 1;
      }
    }

    // address original triangle corners again
    for (int j = 0; j < 3; j++) {
      int vi = corefinedPatches[triIDCur].cornerVertices[j];
      if (newVertexIDs[cylinderMesh.triVtxID(triIDCur, j)] < 0) {
        newVertexIDs[cylinderMesh.triVtxID(triIDCur, j)] = vertexMasks[vi];
      }
    }

    // adding delete boundary info and edge info
    for (int j = 0; j < 3; j++) {
      std::pair<int, int> edgeID{ cylinderMesh.triVtxID(triIDCur, j), cylinderMesh.triVtxID(triIDCur, (j + 1) % 3) };
      if (edgeID.first > edgeID.second)
        std::swap(edgeID.first, edgeID.second);

      // if it is an edge that is visited for the first time
      if (edgeMask[j]) {
        auto ret = edgeVertexIDs.emplace(edgeID, corefinedPatches[triIDCur].edgeVertices[j]);
        PGO_ALOG(ret.second == true);
        auto iter = ret.first;
        for (int k = 0; k < (int)iter->second.size(); k++) {
          iter->second[k] = vertexMasks[iter->second[k]];
        }
      }
      // otherwise, this edge will not be visited any more (becaue the mesh is a manifold)
      else {
        edgeVertexIDs.erase(edgeID);
      }
    }

    // adding cut triangles to the global array
    triangleNewIDs.assign(corefinedPatches[triIDCur].triangles.size(), -1);
    for (int j = 0; j < (int)corefinedPatches[triIDCur].triangles.size(); j++) {
      ES::V3i triID = ES::Mp<const ES::V3i>(corefinedPatches[triIDCur].triangles[j].data());
      for (int k = 0; k < 3; k++) {
        triID[k] = vertexMasks[triID[k]];
      }
      finalTriangles.emplace_back(triID);
      triangleNewIDs[j] = (int)finalTriangles.size() - 1;
    }

    // store edge info
    for (int j = 0; j < (int)corefinedPatches[triIDCur].triangleEdgeIDs.size(); j++) {
      int newTriID = std::get<0>(corefinedPatches[triIDCur].triangleEdgeIDs[j]);
      newTriID = triangleNewIDs[newTriID];

      int ithEdge = std::get<1>(corefinedPatches[triIDCur].triangleEdgeIDs[j]);
      int edgeID = std::get<2>(corefinedPatches[triIDCur].triangleEdgeIDs[j]);
      targetEdgeOverlappingTriangleIDs[edgeID].emplace_back(std::array<int, 2>{ newTriID, ithEdge });
    }

    // store vtx info
    for (int j = 0; j < (int)corefinedPatches[triIDCur].vertices.size(); j++) {
      if (corefinedPatches[triIDCur].vertexIsOnTargetMesh[j][0] >= 0) {
        auto it = vertexIsOnTargetMesh.lower_bound(vertexMasks[j]);
        if (it != vertexIsOnTargetMesh.end() && it->first == vertexMasks[j]) {
        }
        else {
          vertexIsOnTargetMesh.emplace_hint(it, vertexMasks[j], corefinedPatches[triIDCur].vertexIsOnTargetMesh[j]);
        }
      }
    }

    // expand bfs search
    ES::V3i triN = cylinderMeshNeighbor.getTriangleNeighbors(triIDCur);
    for (int j = 0; j < 3; j++) {
      if (triN[j] < 0) {
        continue;
      }

      if (isTriVisited[triN[j]])
        continue;

      triIDQ.emplace_back(triN[j]);
      isTriVisited[triN[j]] = 1;
    }
  }

  t2 = std::chrono::high_resolution_clock::now();
  SPDLOG_LOGGER_INFO(pgo::Logging::lgr(), "S2 Done. Time cost: {:.4f}s", double(std::chrono::duration_cast<std::chrono::microseconds>(t2 - t1).count()) * 1e-6);

  if (dumpMesh) {
    pgo::CGALInterface::Polyhedron<IK> finalMesh;
    pgo::CGALInterface::PolyhedronBuilderNonTriangle<IK, std::vector<IK::Point_3>::iterator, std::vector<ES::V3i>::iterator> builder(
      finalPoints3D.begin(), finalTriangles.begin(), (int)finalPoints3D.size(), (int)finalTriangles.size());
    finalMesh.delegate(builder);
    // std::ofstream("rzr0.off") << finalMesh;

    pgo::Mesh::TriMeshGeo mesh;
    int idx = 0;
    std::map<pgo::CGALInterface::Polyhedron<IK>::Vertex_const_handle, int> vertexIndices;
    for (pgo::CGALInterface::Polyhedron<IK>::Vertex_const_iterator itt = finalMesh.vertices_begin(); itt != finalMesh.vertices_end(); ++itt) {
      const IK::Point_3 &pt = itt->point();
      mesh.addPos(ES::V3d{ tod(pt[0]), tod(pt[1]), tod(pt[2]) });

      vertexIndices.insert(make_pair(itt, idx));
      idx++;
    }

    for (pgo::CGALInterface::Polyhedron<IK>::Face_const_iterator itt = finalMesh.facets_begin(); itt != finalMesh.facets_end(); ++itt) {
      ES::V3i face;
      int inc = 0;

      pgo::CGALInterface::Polyhedron<IK>::Halfedge_around_facet_const_circulator ctr = itt->facet_begin();
      do {
        auto idxItt = vertexIndices.find(ctr->vertex());
        if (idxItt == vertexIndices.end()) {
          std::cerr << "Cannot find the vertex" << std::endl;
          abort();
        }
        if (inc == 3) {
          std::cerr << "Not a triangle mesh" << std::endl;
          abort();
        }

        face[inc++] = idxItt->second;
        ctr++;
      } while (ctr != itt->facet_begin());
      mesh.addTri(face);
    }

    mesh.save("rzr0.obj");
  }

  // remove triangles that are not necessary
  using I2 = std::array<int, 2>;
  using I2Hash = ES::IntArrayHash<2>;
  using I2Equal = ES::IntArrayEqual<2>;
  using I2Less = ES::IntArrayLess<2>;
  auto makeKey = [](int tri, int j) -> I2 {
    return I2{ tri, j };
  };

  // compute a map, where each triangle edge is mapped to the edge ID of the target mesh
  std::unordered_map<I2, int, I2Hash, I2Equal> triangleTargetEdgeIDs;
  int potentialTotalCount = 0;
  for (int ei = 0; ei < (int)targetEdgeOverlappingTriangleIDs.size(); ei++) {
    potentialTotalCount += (int)targetEdgeOverlappingTriangleIDs[ei].size();
  }
  triangleTargetEdgeIDs.reserve(potentialTotalCount);
  for (int ei = 0; ei < (int)targetEdgeOverlappingTriangleIDs.size(); ei++) {
    for (const auto &triID : targetEdgeOverlappingTriangleIDs[ei]) {
      triangleTargetEdgeIDs.emplace(triID, ei);
    }
  }

  // build a connectivity structure
  pgo::Mesh::TriMeshGeo dummyCorefinedMesh;
  dummyCorefinedMesh.positions().assign(finalPoints3D.size(), pgo::asVec3d(0.0));
  dummyCorefinedMesh.triangles() = finalTriangles;
  pgo::Mesh::TriMeshNeighbor corefinedMeshNeighbor(dummyCorefinedMesh);

  // flag that markers whether each triangle is visited
  // <0: not visit; >= 0; the set ID for each triangle
  std::vector<int> isFinalTriVisited(finalTriangles.size(), -1);

  // the bfs Q
  std::vector<int> Q;
  Q.reserve(10000);

  // the edges the triangle has hit
  std::vector<int> hitEdgeIDs;
  hitEdgeIDs.reserve(10000);

  // initial valid set ID
  std::vector<int> validSetIDs;
  validSetIDs.reserve(targetMesh.numTriangles());

  int setID = 0;
  while (1) {
    Q.clear();
    hitEdgeIDs.clear();

    // find the first triangle with edge overlapping with the target mesh
    bool found = false;
    for (int ti = 0; ti < (int)finalTriangles.size(); ti++) {
      // is visited
      if (isFinalTriVisited[ti] >= 0)
        continue;

      // is a triangle that overlaps the target edge
      if (triangleTargetEdgeIDs.find(makeKey(ti, 0)) == triangleTargetEdgeIDs.end() &&
        triangleTargetEdgeIDs.find(makeKey(ti, 1)) == triangleTargetEdgeIDs.end() &&
        triangleTargetEdgeIDs.find(makeKey(ti, 2)) == triangleTargetEdgeIDs.end())
        continue;

      Q.emplace_back(ti);
      isFinalTriVisited[ti] = setID;
      found = true;
      break;
    }

    if (found == false)
      break;

    size_t start = 0;
    while (start < Q.size()) {
      int triCurID = Q[start++];

      // for each valid set, maximal 3 different edges can be hit,
      // if exceed, we early quit
      if (hitEdgeIDs.size() < 4ull) {
        for (int j = 0; j < 3; j++) {
          // if it is a boundary triangle
          auto it = triangleTargetEdgeIDs.find(makeKey(triCurID, j));
          if (it != triangleTargetEdgeIDs.end()) {
            auto hitEdge_it = std::find(hitEdgeIDs.begin(), hitEdgeIDs.end(), it->second);
            // if it is a new edge ID
            if (hitEdge_it == hitEdgeIDs.end()) {
              hitEdgeIDs.emplace_back(it->second);
            }
          }
        }
      }

      for (int j = 0; j < 3; j++) {
        int ntriID = corefinedMeshNeighbor.getTriangleNeighbors(triCurID)[j];
        if (ntriID < 0)
          continue;

        // is visited
        if (isFinalTriVisited[ntriID] >= 0)
          continue;

        // if it is a boundary
        if (triangleTargetEdgeIDs.find(makeKey(triCurID, j)) != triangleTargetEdgeIDs.end())
          continue;

        Q.emplace_back(ntriID);
        isFinalTriVisited[ntriID] = setID;
      }
    }

    if (hitEdgeIDs.size() == 3ull) {
      std::array<int, 3> edgeIDSet{ hitEdgeIDs[0], hitEdgeIDs[1], hitEdgeIDs[2] };
      std::sort(edgeIDSet.begin(), edgeIDSet.end());
      auto it = selectedTrianglesQueryByEdgeIDs.find(edgeIDSet);
      if (it != selectedTrianglesQueryByEdgeIDs.end()) {
        validSetIDs.emplace_back(setID);
      }
      else {
        SPDLOG_LOGGER_DEBUG(pgo::Logging::lgr(), "Find a false positive edge set. ({}, {}, {})", edgeIDSet[0], edgeIDSet[1], edgeIDSet[2]);
      }
    }

    setID++;
  }

  if (dumpMesh) {
    std::vector<ES::V3d> kk;
    for (int i = 0; i < (int)finalPoints3D.size(); i++) {
      kk.emplace_back(tod(finalPoints3D[i].x()), tod(finalPoints3D[i].y()), tod(finalPoints3D[i].z()));
    }

    pgo::Mesh::TriMeshGeo asas;
    asas.positions() = kk;

    for (int si : validSetIDs) {
      for (int ti = 0; ti < (int)finalTriangles.size(); ti++) {
        if (isFinalTriVisited[ti] == si) {
          asas.addTri(finalTriangles[ti]);
        }
      }
    }

    asas.save(fmt::format("rzr1.obj"));
  }

  std::vector<int> setIDAll(setID, 0);
  std::iota(setIDAll.begin(), setIDAll.end(), 0);

  // compute the set ID of triangles that cannot be edited
  std::vector<int> invalidSetIDs;
  std::set_difference(setIDAll.begin(), setIDAll.end(),
    validSetIDs.begin(), validSetIDs.end(),
    std::back_inserter(invalidSetIDs));

  // compute the vertices that cannot be edited
  std::vector<int> unremovableVertex(finalPoints3D.size(), 0);
  for (int ti = 0; ti < (int)finalTriangles.size(); ti++) {
    if (std::binary_search(validSetIDs.begin(), validSetIDs.end(), isFinalTriVisited[ti]))
      continue;

    for (int j = 0; j < 3; j++) {
      unremovableVertex[finalTriangles[ti][j]] = 1;
    }
  }

  std::vector<std::vector<int>> finalFaces;
  finalFaces.reserve(finalTriangles.size());

  std::set<I2, I2Less> edgeHashTable;
  std::vector<int> edgeLoops;
  std::vector<int> edgeLoops1;

  std::vector<int> isFaceFlexible;
  isFaceFlexible.reserve(finalTriangles.size());

  std::cout << finalTriangles.size() << ',' << validSetIDs.size() << ',' << invalidSetIDs.size() << std::endl;

  for (int si : validSetIDs) {
    edgeHashTable.clear();
    for (int ti = 0; ti < (int)finalTriangles.size(); ti++) {
      if (isFinalTriVisited[ti] != si) {
        continue;
      }

      for (int j = 0; j < 3; j++) {
        I2 edgeID{ finalTriangles[ti][j], finalTriangles[ti][(j + 1) % 3] };
        if (edgeID[0] > edgeID[1])
          std::swap(edgeID[0], edgeID[1]);

        auto it = edgeHashTable.lower_bound(edgeID);
        if (it != edgeHashTable.end() && *it == edgeID) {
          edgeHashTable.erase(it);
        }
        else {
          edgeHashTable.emplace_hint(it, edgeID);
        }
      }
    }

    // compute edge loop
    edgeLoops.clear();
    edgeLoops.emplace_back(edgeHashTable.begin()->at(0));
    while (edgeHashTable.size()) {
      for (auto it = edgeHashTable.begin(); it != edgeHashTable.end(); ++it) {
        if (it->at(0) == edgeLoops.back()) {
          edgeLoops.emplace_back(it->at(1));
          edgeHashTable.erase(it);
          break;
        }
        else if (it->at(1) == edgeLoops.back()) {
          edgeLoops.emplace_back(it->at(0));
          edgeHashTable.erase(it);
          break;
        }
      }
    }

    PGO_ALOG(edgeLoops.front() == edgeLoops.back());
    edgeLoops.pop_back();

    // compute corner vertex
    int ctr = 0;
    std::array<int, 3> cornerVertexIDs = { -1, -1, -1 };
    for (int vid : edgeLoops) {
      auto it = vertexIsOnTargetMesh.find(vid);
      if (it != vertexIsOnTargetMesh.end() &&
        it->second[0] == 0) {
        if (ctr >= 3) {
          pgo::Mesh::TriMeshGeo mm;
          for (int vi : edgeLoops) {
            IK::Point_3 p = finalPoints3D[vi];
            mm.addPos(ES::V3d(tod(p[0]), tod(p[1]), tod(p[2])));
          }
          mm.save("bug.obj");

          abort();
        }
        cornerVertexIDs[ctr] = vid;
        ctr++;
      }
    }

    // std::cout << ctr << std::endl;
    // std::cout << cornerVertexIDs[0] << ',' << cornerVertexIDs[0] << ',' << cornerVertexIDs[0] << std::endl;
    if (ctr == 1) {
      edgeLoops1.assign(edgeLoops.begin(), edgeLoops.end());
    }
    else if (ctr == 2) {
      // for each edge
      int ith = 0;
      edgeLoops1.clear();

      // find first point in the loop
      while (edgeLoops[ith] != cornerVertexIDs[0]) {
        ith = (ith + 1) % (int)edgeLoops.size();
      }

      // find next point in the loop
      int ithNext = ith;
      while (edgeLoops[ithNext] != cornerVertexIDs[1]) {
        ithNext = (ithNext + 1) % (int)edgeLoops.size();
      }

      int ith1 = (ith + 1) % (int)edgeLoops.size();
      int count0 = 0, count1 = 0;
      int lastEdgeID = -1;
      while (ith1 != ithNext) {
        auto it = vertexIsOnTargetMesh.find(edgeLoops[ith1]);
        if (it != vertexIsOnTargetMesh.end() && it->second[0]) {
          if (lastEdgeID < 0) {
            lastEdgeID = it->second[1];
            count1++;
          }
          else if (lastEdgeID == it->second[1]) {
            count1++;
          }
        }
        count0++;
        ith1 = (ith1 + 1) % (int)edgeLoops.size();
      }
      std::cout << count0 << ',' << count1 << std::endl;

      if (count0 != count1) {
        std::cout << "mismatch" << count1 << '/' << count0 << std::endl;

        std::reverse(edgeLoops.begin(), edgeLoops.end());

        ith = 0;
        // find first point in the loop
        while (edgeLoops[ith] != cornerVertexIDs[0]) {
          ith = (ith + 1) % (int)edgeLoops.size();
        }

        // find next point in the loop
        ithNext = ith;
        while (edgeLoops[ithNext] != cornerVertexIDs[1]) {
          ithNext = (ithNext + 1) % (int)edgeLoops.size();
        }

        ith1 = (ith + 1) % (int)edgeLoops.size();
        count0 = 0, count1 = 0;
        lastEdgeID = -1;
        while (ith1 != ithNext) {
          auto it = vertexIsOnTargetMesh.find(edgeLoops[ith1]);
          if (it != vertexIsOnTargetMesh.end() && it->second[0]) {
            if (lastEdgeID < 0) {
              lastEdgeID = it->second[1];
              count1++;
            }
            else if (lastEdgeID == it->second[1]) {
              count1++;
            }
          }
          count0++;
          ith1 = (ith1 + 1) % (int)edgeLoops.size();
        }
        PGO_ALOG(count0 == count1);
      }

      edgeLoops1.emplace_back(edgeLoops[ith]);
      ith1 = (ith + 1) % (int)edgeLoops.size();
      while (ith1 != ithNext) {
        if (unremovableVertex[edgeLoops[ith1]]) {
          edgeLoops1.emplace_back(edgeLoops[ith1]);
        }
        ith1 = (ith1 + 1) % (int)edgeLoops.size();
      }

      while (ith1 != ith) {
        edgeLoops1.emplace_back(edgeLoops[ith1]);
        ith1 = (ith1 + 1) % (int)edgeLoops.size();
      }
    }
    else if (ctr == 3) {
      // for each edge
      int ith = 0;
      edgeLoops1.clear();
      for (int j = 0; j < 3; j++) {
        // find first point in the loop
        while (edgeLoops[ith] != cornerVertexIDs[j]) {
          ith = (ith + 1) % (int)edgeLoops.size();
        }

        // find next point in the loop
        int ithNext = ith;
        while (edgeLoops[ithNext] != cornerVertexIDs[(j + 1) % 3]) {
          ithNext = (ithNext + 1) % (int)edgeLoops.size();
        }

        edgeLoops1.emplace_back(edgeLoops[ith]);
        int ith1 = (ith + 1) % (int)edgeLoops.size();
        while (ith1 != ithNext) {
          if (unremovableVertex[edgeLoops[ith1]]) {
            edgeLoops1.emplace_back(edgeLoops[ith1]);
          }
          ith1 = (ith1 + 1) % (int)edgeLoops.size();
        }

        ith = ithNext;
      }
    }

    finalFaces.emplace_back(edgeLoops1);
    if (edgeLoops1.size() == 3ull && ctr == 3) {
      isFaceFlexible.emplace_back(0);
    }
    else {
      isFaceFlexible.emplace_back(1);
    }
  }

  invalidSetIDs.emplace_back(-1);
  for (int si : invalidSetIDs) {
    for (int ti = 0; ti < (int)finalTriangles.size(); ti++) {
      if (isFinalTriVisited[ti] != si) {
        continue;
      }

      finalFaces.emplace_back(finalTriangles[ti].begin(), finalTriangles[ti].end());
      isFaceFlexible.emplace_back(1);
    }
  }

  std::unordered_map<int, int> usedVertices;
  for (const auto &f : finalFaces) {
    for (int vid : f) {
      usedVertices.emplace(vid, 0);
    }
  }

  std::vector<K::Point_3> finalPoints1;
  finalPoints1.reserve(finalPoints3D.size());
  for (int inc = 0; auto &v : usedVertices) {
    v.second = inc++;

    const IK::Point_3 &p = finalPoints3D[v.first];
    finalPoints1.emplace_back(toP3_EK(p));
  }

  std::unordered_map<int, int> vtxIDNew2Old;
  vtxIDNew2Old.reserve(usedVertices.size());
  for (auto &pr : usedVertices) {
    vtxIDNew2Old.emplace(pr.second, pr.first);
  }
  std::cout << usedVertices.size() << ',' << finalPoints1.size() << ',' << vtxIDNew2Old.size() << std::endl;

  for (auto &f : finalFaces) {
    for (int &vid : f) {
      auto it = usedVertices.find(vid);
      PGO_ALOG(it != usedVertices.end());
      vid = it->second;
    }
  }

  Poly finalMesh;
  orientMesh(finalPoints1, finalFaces, finalMesh);

  int orientationCount[2] = { 0, 0 };
  for (Poly::Facet_handle fi = finalMesh.facets_begin(); fi != finalMesh.facets_end(); ++fi) {
    if (fi->is_triangle()) {
      Poly::Halfedge_handle h0 = fi->halfedge();
      Poly::Halfedge_handle h1 = h0->next();
      Poly::Halfedge_handle h2 = h1->next();

      K::Vector_3 e01 = h1->vertex()->point() - h0->vertex()->point();
      K::Vector_3 e02 = h2->vertex()->point() - h0->vertex()->point();
      K::Vector_3 n = CGAL::cross_product(e01, e02);

      if (n.squared_length() < 1e-10) {
        continue;
      }

      K::Point_3 center = h0->vertex()->point() + e01 / 3.0 + e02 / 3.0;
      auto projPt = computeProjectionPt(center, primitiveCenterEK);
      K::Vector_3 dir = center - std::get<0>(projPt);

      if (CGAL::scalar_product(dir, n) < 0) {
        orientationCount[1] += 1;
      }
      else {
        orientationCount[0] += 1;
      }
    }
  }

  if (orientationCount[0] < orientationCount[1]) {
    CGAL::Polygon_mesh_processing::reverse_face_orientations(finalMesh);
  }

  if constexpr (dumpMesh) {
    std::ofstream("rzr1.off") << finalMesh;
    std::ofstream outfile("rzr.obj");
    for (auto vit = finalMesh.vertices_begin(); vit != finalMesh.vertices_end(); ++vit) {
      ES::V3d p = toVec3(vit->point());
      outfile << "v " << p[0] << ' ' << p[1] << ' ' << p[2] << std::endl;
    }

    for (auto fit = finalMesh.facets_begin(); fit != finalMesh.facets_end(); ++fit) {
      Poly::Halfedge_handle h0 = fit->halfedge();
      Poly::Halfedge_handle h1 = h0;
      outfile << "f";
      do {
        outfile << ' ' << (h1->vertex()->id() + 1);
        h1 = h1->next();
      } while (h1 != h0);
      outfile << std::endl;
    }
    outfile.close();
  }

  // bool triangulateRet = CGAL::Polygon_mesh_processing::triangulate_faces(finalMesh);
  // if (triangulateRet == false) {
  if (triangulateMesh(finalMesh) == false) {
    std::cout << "First iteration cannot clear the non triangles." << std::endl;

    for (auto fi = finalMesh.facets_begin(); fi != finalMesh.facets_end(); ++fi) {
      if (fi->is_triangle())
        continue;

      Poly::Halfedge_handle h0 = fi->halfedge();
      Poly::Halfedge_handle h1 = h0;
      std::vector<std::array<V3, 2>> uvws;
      uvws.reserve(100);
      do {
        const K::Point_3 &pt = h1->vertex()->point();
        auto projPtRet = computeProjectionPt(pt, primitiveCenterEK);
        V3 uvw0 = computeCylinderCoordinate(pt, projPtRet, ax_x, ax_y, primitiveCenterIK);
        V3 uvw1 = uvw0;
        uvw1[0] = EK_to_IK<real>(rotateAngle(IK_to_EK(uvw0[0])));

        uvws.emplace_back(std::array<V3, 2>{ uvw0, uvw1 });
        h1 = h1->next();
      } while (h0 != h1);

      real maxDiff[2] = { 0, 0 };
      for (int ci = 0; ci < 2; ci++) {
        for (int i = 0; i < (int)uvws.size(); i++) {
          for (int j = i + 1; j < (int)uvws.size(); j++) {
            maxDiff[ci] = std::max(maxDiff[ci], abs(uvws[i][ci][0] - uvws[j][ci][0]));
          }
        }
      }

      int projID = 0;
      if (maxDiff[0] < maxDiff[1]) {
        projID = 0;
      }
      else {
        projID = 1;
      }

      V3 center(0, 0, 0);
      for (int i = 0; i < (int)uvws.size(); i++) {
        center += uvws[i][projID];
      }
      center /= real((int)uvws.size());

      V3 x;
      if (projID == 0) {
        real angle = center[0];
        x = (ax_x * cos(angle) + ax_y * sin(angle)) * center[2];
      }
      else {
        real angle = center[0];
        x = -(ax_x * cos(angle) + ax_y * sin(angle)) * center[2];
      }

      V3 vert = axisDiff * center[1] + primitiveCenterIK[0];
      x += vert;

      int vidCur = (int)finalMesh.size_of_vertices();
      Poly::Halfedge_handle h = finalMesh.create_center_vertex(fi->halfedge());
      h->vertex()->point() = toP3_EK(x);
      h->vertex()->id() = vidCur;
    }
    int nv = finalMesh.size_of_vertices();
    int max_vid = 0;
    for (auto vi = finalMesh.vertices_begin(); vi != finalMesh.vertices_end(); ++vi) {
      max_vid = std::max(max_vid, (int)vi->id());
    }

    std::cout << "ID cmp:" << max_vid << ',' << nv << std::endl;
    std::cout << "Is triangle mesh: " << finalMesh.is_pure_triangle() << std::endl;
    std::ofstream("tri1.off") << finalMesh;
  }
  else {
    pgo::Mesh::TriMeshGeo meshOut;
    pgo::CGALInterface::polyhedron2TriangleMesh(finalMesh, meshOut);
    // meshOut.save("rzr2.obj");

    pgo::Mesh::TriMeshNeighbor m(meshOut);
    std::cout << "Is triangle mesh: " << finalMesh.is_pure_triangle() << std::endl;
  }
  // }

  cylinderVertexProjection(finalMesh, targetMesh, targetMeshPoly, targetMeshBVTreeExact,
    primitiveCenterEK, maxAllowedRadiusDist, vtxIDNew2Old, vertexIsOnTargetMesh, selectedEdges);

  if constexpr (dumpMesh) {
    std::ofstream("rzr2.off") << finalMesh;
  }
  // exit(1);

  std::unordered_set<int> borderVertexIDs;
  int numVertices = (int)finalMesh.size_of_vertices();

  // create a closed mesh
  if constexpr (1) {
    // find loop 0
    Poly::Halfedge_handle border_h[2];
    int cids[2];
    std::set<Poly::Halfedge_handle> border_h_set[2];
    for (auto hit = finalMesh.halfedges_begin(); hit != finalMesh.halfedges_end(); ++hit) {
      if (hit->is_border()) {
        border_h[0] = hit;
        break;
      }
    }

    Poly::Halfedge_handle hh = border_h[0];
    do {
      border_h_set[0].emplace(hh);
      hh = hh->next();
    } while (hh != border_h[0]);

    // find loop 1
    for (auto hit = finalMesh.halfedges_begin(); hit != finalMesh.halfedges_end(); ++hit) {
      if (border_h_set[0].find(hit) != border_h_set[0].end())
        continue;

      if (hit->is_border()) {
        border_h[1] = hit;
        break;
      }
    }

    hh = border_h[1];
    do {
      border_h_set[1].emplace(hh);
      hh = hh->next();
    } while (hh != border_h[1]);

    for (int i = 0; i < 2; i++) {
      for (const auto &pr : border_h_set[i]) {
        borderVertexIDs.emplace(pr->vertex()->id());
      }
    }

    auto proj0 = computeProjectionPt(border_h[0]->vertex()->point(), primitiveCenterEK);
    if ((std::get<0>(proj0) - primitiveCenterEK[0]).squared_length() <
      (std::get<0>(proj0) - primitiveCenterEK[1]).squared_length()) {
      cids[0] = 0;
      cids[1] = 1;
    }
    else {
      cids[0] = 1;
      cids[1] = 0;
    }

    int vidCur = (int)finalMesh.size_of_vertices();
    std::cout << vidCur << std::endl;

    for (int j = 0; j < 2; j++) {
      Poly::Halfedge_handle h_next = border_h[j]->next();
      Poly::Halfedge_handle h = finalMesh.add_vertex_and_facet_to_border(border_h[j]->prev(), border_h[j]);
      h->vertex()->point() = primitiveCenterEK[cids[j]];
      h->vertex()->id() = vidCur++;

      Poly::Halfedge_handle h0 = h->next()->opposite();

      do {
        Poly::Halfedge_handle h1 = h_next->next();
        finalMesh.add_facet_to_border(h0, h_next);

        h_next = h1;
      } while (h_next->next() != h0);
      finalMesh.fill_hole(h_next);
    }
    std::cout << vidCur << std::endl;
  }

  if constexpr (dumpMesh) {
    std::ofstream("rzr3.off") << finalMesh;
  }

  EdgeMap targetEdges;
  for (int ti = 0; ti < targetMesh.numTriangles(); ti++) {
    for (int j = 0; j < 3; j++) {
      std::pair edgeID{ targetMesh.triVtxID(ti, j), targetMesh.triVtxID(ti, (j + 1) % 3) };
      if (edgeID.first > edgeID.second) {
        std::swap(edgeID.first, edgeID.second);
      }
      targetEdges.emplace(edgeID);
    }
  }

  std::vector<VertexProperty> finalMeshVertexProperties(finalMesh.size_of_vertices());
  for (auto vit = finalMesh.vertices_begin(); vit != finalMesh.vertices_end(); ++vit) {
    int vid = (int)vit->id();
    if (vid >= numVertices) {
      finalMeshVertexProperties[vid].vaType = VertexAttachmentType::SKELETON_VERTEX;
      continue;
    }

    auto it = vtxIDNew2Old.find(vid);
    if (it != vtxIDNew2Old.end()) {
      int oldVID = it->second;

      auto itType = vertexIsOnTargetMesh.find(oldVID);
      if (itType != vertexIsOnTargetMesh.end()) {
        if (itType->second[0] == 0) {
          finalMeshVertexProperties[vid].vaType = VertexAttachmentType::TARGET_VERTEX;
          finalMeshVertexProperties[vid].attachedTargets[0] = vit->point();
          finalMeshVertexProperties[vid].targetIDs[0] = itType->second[1];
        }
        else if (itType->second[0] == 1) {
          finalMeshVertexProperties[vid].vaType = VertexAttachmentType::TARGET_EDGE;
          finalMeshVertexProperties[vid].attachedTargets[0] = toP3_EK(targetMesh.pos(selectedEdges[itType->second[1]].first));
          finalMeshVertexProperties[vid].attachedTargets[1] = toP3_EK(targetMesh.pos(selectedEdges[itType->second[1]].second));
          finalMeshVertexProperties[vid].targetIDs[0] = selectedEdges[itType->second[1]].first;
          finalMeshVertexProperties[vid].targetIDs[1] = selectedEdges[itType->second[1]].second;
        }
        else {
          throw std::runtime_error("impossible");
        }
      }
    }

    if (borderVertexIDs.find(vid) != borderVertexIDs.end()) {
      finalMeshVertexProperties[vid].isBorder = 1;
    }
  }

  cleanUpMesh(targetMeshPoly, targetMeshBVTreeExact, finalMesh, 5e-3, finalMeshVertexProperties, targetEdges);

  double step[2] = { 0.5, -0.6 };
  for (int iter = 0; iter < 0; iter++) {
    for (auto vit = finalMesh.vertices_begin(); vit != finalMesh.vertices_end(); ++vit) {
      int vid = (int)vit->id();

      // pgo::asVec3d pp(-0.141025, 0.0707617, -0.0211217);
      // pgo::asVec3d p = toVec3(vit->point());
      // if (len(p - pp) < 1e-3) {
      //   std::cout << p << std::endl;
      //   std::cout << "z";
      // }

      if (!(finalMeshVertexProperties[vid].vaType == VertexAttachmentType::ATTACHED) &&
        !(finalMeshVertexProperties[vid].vaType == VertexAttachmentType::DETACHED) &&
        !(finalMeshVertexProperties[vid].vaType == VertexAttachmentType::TARGET_EDGE))
        continue;

      if (finalMeshVertexProperties[vid].isBorder)
        continue;

      if (finalMeshVertexProperties[vid].vaType == VertexAttachmentType::TARGET_EDGE) {
        Poly::Halfedge_handle h = vit->halfedge();
        Poly::Halfedge_handle h1 = h;

        K::Vector_3 nprev;
        bool firstTime = true;
        bool found = false;
        do {
          const K::Point_3 &p0 = h1->vertex()->point();
          const K::Point_3 &p1 = h1->next()->vertex()->point();
          const K::Point_3 &p2 = h1->next()->next()->vertex()->point();

          K::Vector_3 e01 = p1 - p0;
          K::Vector_3 e02 = p2 - p0;
          K::Vector_3 n = CGAL::cross_product(e01, e02);

          if (firstTime) {
            nprev = n;
            firstTime = false;
          }
          else {
            if (CGAL::scalar_product(nprev, n) < 0) {
              found = true;
              break;
            }
          }

          h1 = h1->next()->opposite();
        } while (h != h1);

        if (found == false)
          continue;
      }

      Poly::Halfedge_handle h = vit->halfedge();
      Poly::Halfedge_handle h1 = h;
      std::vector<std::array<V3, 2>> uvws;
      uvws.reserve(100);
      do {
        const K::Point_3 &pt = h1->opposite()->vertex()->point();
        auto projPtRet = computeProjectionPt(pt, primitiveCenterEK);
        V3 uvw0 = computeCylinderCoordinate(pt, projPtRet, ax_x, ax_y, primitiveCenterIK);
        V3 uvw1 = uvw0;
        uvw1[0] = EK_to_IK<real>(rotateAngle(IK_to_EK(uvw0[0])));
        uvws.emplace_back(std::array<V3, 2>{ uvw0, uvw1 });

        h1 = h1->next()->opposite();
      } while (h != h1);

      real maxDiff[2] = { 0, 0 };
      for (int ci = 0; ci < 2; ci++) {
        for (int i = 0; i < (int)uvws.size(); i++) {
          for (int j = i + 1; j < (int)uvws.size(); j++) {
            maxDiff[ci] = std::max(maxDiff[ci], abs(uvws[i][ci][0] - uvws[j][ci][0]));
          }
        }
      }

      int projID = 0;
      if (maxDiff[0] < maxDiff[1]) {
        projID = 0;
      }
      else {
        projID = 1;
      }

      V3 center(0, 0, 0);
      for (int i = 0; i < (int)uvws.size(); i++) {
        center += uvws[i][projID];
      }
      center /= real((int)(uvws.size()));

      V3 x;
      if (projID == 0) {
        real angle = center[0];
        x = (ax_x * cos(angle) + ax_y * sin(angle)) * center[2];
      }
      else {
        real angle = center[0];
        x = -(ax_x * cos(angle) + ax_y * sin(angle)) * center[2];
      }

      V3 vert = axisDiff * center[1] + primitiveCenterIK[0];
      x += vert;

      K::Vector_3 diff = toP3_EK(x) - vit->point();
      vit->point() = vit->point() + diff * 1;  // step[iter % 2];
    }
  }

  if constexpr (dumpMesh) {
    std::ofstream("rzr4.off") << finalMesh;
  }

  pgo::CGALInterface::polyhedron2TriangleMesh(finalMesh, meshOut);

  if (filename) {
    saveExact(finalMesh, filename);
  }

  // exit(1);
  // orientMesh(finalPoints1, finalFaces, cer, 1);

  /*
  std::ofstream aa("qq.obj");
  for (int vi = 0; vi < finalPoints3D.size(); vi++) {
    aa << "v " << tod(finalPoints3D[vi][0]) << ' '
       << tod(finalPoints3D[vi][1]) << ' '
       << tod(finalPoints3D[vi][2]) << '\n';
  }

  for (int fi = 0; fi < (int)finalFaces.size(); fi++) {
    if (finalFaces[fi].size() < 3) {
      std::cout << fi << ' ' << finalFaces[fi].size() << std::endl;
    }
    aa << 'f';
    for (int vi : finalFaces[fi]) {
      aa << ' ' << vi + 1;
    }
    aa << '\n';
  }
  aa.close();

  aa.clear();
  aa.open("qq1.obj");
  for (int vi = 0; vi < finalPoints3D.size(); vi++) {
    aa << "v " << tod(finalPoints3D[vi][0]) << ' '
       << tod(finalPoints3D[vi][1]) << ' '
       << tod(finalPoints3D[vi][2]) << '\n';
  }

  for (int fi = 0; fi < (int)finalTriangles.size(); fi++) {
    aa << 'f';
    for (int vi : finalTriangles[fi]) {
      aa << ' ' << vi + 1;
    }
    aa << '\n';
  }
  aa.close();
  */

  return 0;
}

int MedialAxisRepresentation::fillCylinder(const pgo::Mesh::TriMeshGeo &sphereMesh, const ES::V3d centers[2], pgo::Mesh::TriMeshGeo &meshOut)
{
  K::Point_3 primitiveCenterEK[2] = {
    K::Point_3(centers[0][0], centers[0][1], centers[0][2]),
    K::Point_3(centers[1][0], centers[1][1], centers[1][2]),
  };

  Poly finalMesh;
  pgo::CGALInterface::triangleMesh2Polyhedron(sphereMesh, finalMesh);

  // create a closed mesh
  if constexpr (1) {
    // find loop 0
    Poly::Halfedge_handle border_h[2];
    int cids[2];
    std::set<Poly::Halfedge_handle> border_h_set[2];
    for (auto hit = finalMesh.halfedges_begin(); hit != finalMesh.halfedges_end(); ++hit) {
      if (hit->is_border()) {
        border_h[0] = hit;
        break;
      }
    }

    Poly::Halfedge_handle hh = border_h[0];
    do {
      border_h_set[0].emplace(hh);
      hh = hh->next();
    } while (hh != border_h[0]);

    // find loop 1
    for (auto hit = finalMesh.halfedges_begin(); hit != finalMesh.halfedges_end(); ++hit) {
      if (border_h_set[0].find(hit) != border_h_set[0].end())
        continue;

      if (hit->is_border()) {
        border_h[1] = hit;
        break;
      }
    }

    hh = border_h[1];
    do {
      border_h_set[1].emplace(hh);
      hh = hh->next();
    } while (hh != border_h[1]);

    auto proj0 = computeProjectionPt(border_h[0]->vertex()->point(), primitiveCenterEK);
    if ((std::get<0>(proj0) - primitiveCenterEK[0]).squared_length() <
      (std::get<0>(proj0) - primitiveCenterEK[1]).squared_length()) {
      cids[0] = 0;
      cids[1] = 1;
    }
    else {
      cids[0] = 1;
      cids[1] = 0;
    }

    int vidCur = (int)finalMesh.size_of_vertices();
    std::cout << vidCur << std::endl;

    for (int j = 0; j < 2; j++) {
      Poly::Halfedge_handle h_next = border_h[j]->next();
      Poly::Halfedge_handle h = finalMesh.add_vertex_and_facet_to_border(border_h[j]->prev(), border_h[j]);
      h->vertex()->point() = primitiveCenterEK[cids[j]];
      h->vertex()->id() = vidCur++;

      Poly::Halfedge_handle h0 = h->next()->opposite();

      do {
        Poly::Halfedge_handle h1 = h_next->next();
        finalMesh.add_facet_to_border(h0, h_next);

        h_next = h1;
      } while (h_next->next() != h0);
      finalMesh.fill_hole(h_next);
    }
    std::cout << vidCur << std::endl;
  }

  pgo::CGALInterface::polyhedron2TriangleMesh(finalMesh, meshOut);

  return 0;
}