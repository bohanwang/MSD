#include "prismCorefinement.h"
#include "corefinementUtilities.h"

#include "pgoLogging.h"
#include "libiglInterface.h"

#include "basicAlgorithms.h"
#include "createTriMesh.h"
#include "triMeshNeighbor.h"
#include "EigenSupport.h"

#include <CGAL/convex_hull_3.h>
#include <CGAL/Side_of_triangle_mesh.h>

#include <tbb/enumerable_thread_specific.h>
#include <tbb/parallel_for.h>
#include <tbb/concurrent_hash_map.h>
#include <tbb/spin_mutex.h>

#include <chrono>
#include <unordered_map>
#include <unordered_set>

namespace MedialAxisRepresentation
{
K::Point_3 dummyFarPt(10, 10, 10);

std::tuple<K::Point_3, int> computePrismProjectionPt(const K::Point_3 &pt, K::Triangle_3 &skeletonTri)
{
  using Primitive = CGAL::AABB_triangle_primitive<K, K::Triangle_3 *>;

  Primitive prim(&skeletonTri);
  CGAL::AABB_traits<K, Primitive> aabbTraits;
  K::Point_3 closestPt = aabbTraits.closest_point_object()(pt, prim, dummyFarPt);

  if (closestPt == skeletonTri[0]) {
    K::Vector_3 e01 = skeletonTri[0] - skeletonTri[1];
    K::Vector_3 e02 = skeletonTri[0] - skeletonTri[2];
    K::Vector_3 dir = pt - closestPt;

    if (CGAL::scalar_product(e01, dir) == 0) {
      return std::tuple(closestPt, 0);
    }
    else if (CGAL::scalar_product(e02, dir) == 0) {
      return std::tuple(closestPt, 2);
    }
    else {
      return std::tuple(closestPt, -1);
    }
  }
  else if (closestPt == skeletonTri[1]) {
    K::Vector_3 e10 = skeletonTri[1] - skeletonTri[0];
    K::Vector_3 e12 = skeletonTri[1] - skeletonTri[2];
    K::Vector_3 dir = pt - closestPt;

    if (CGAL::scalar_product(e10, dir) == 0) {
      return std::tuple(closestPt, 0);
    }
    else if (CGAL::scalar_product(e12, dir) == 0) {
      return std::tuple(closestPt, 1);
    }
    else {
      return std::tuple(closestPt, -1);
    }
  }
  else if (closestPt == skeletonTri[2]) {
    K::Vector_3 e20 = skeletonTri[2] - skeletonTri[0];
    K::Vector_3 e21 = skeletonTri[2] - skeletonTri[1];
    K::Vector_3 dir = pt - closestPt;

    if (CGAL::scalar_product(e20, dir) == 0) {
      return std::tuple(closestPt, 2);
    }
    else if (CGAL::scalar_product(e21, dir) == 0) {
      return std::tuple(closestPt, 1);
    }
    else {
      return std::tuple(closestPt, -1);
    }
  }
  else {
    K::Segment_3 s01(skeletonTri[0], skeletonTri[1]);
    K::Segment_3 s12(skeletonTri[1], skeletonTri[2]);
    K::Segment_3 s20(skeletonTri[2], skeletonTri[0]);

    if (s01.has_on(closestPt)) {
      return std::tuple(closestPt, 0);
    }
    else if (s12.has_on(closestPt)) {
      return std::tuple(closestPt, 1);
    }
    else if (s20.has_on(closestPt)) {
      return std::tuple(closestPt, 2);
    }
    else {
      return std::tuple(closestPt, 3);
    }
  }
}

K::Point_3 computePlaneProjectionPt(const K::Point_3 &pt, K::Triangle_3 &skeletonTri, const K::Vector_3 &nEK)
{
  // N(p - p0) = 0
  //
  K::Vector_3 diff = pt - skeletonTri[0];
  K::Vector_3 nDiff = CGAL::scalar_product(diff, nEK) / nEK.squared_length() * nEK;
  K::Point_3 projPt = pt - nDiff;

  return projPt;
}

std::tuple<K::Point_3, int> computeCylinderProjectionPt(const K::Point_3 &pt, const K::Point_3 primitiveCenter[3], int j)
{
  K::Vector_3 vec = primitiveCenter[(j + 1) % 3] - primitiveCenter[j];
  K::Vector_3 dir = pt - primitiveCenter[j];

  // vecT(t vec + v0 - x)
  K::FT b = vec.squared_length();
  K::FT a = CGAL::scalar_product(vec, dir);
  K::FT t = a / b;
  K::Point_3 closestPt = primitiveCenter[j] + vec * t;

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

V3 computeCylinderCoordinate(const K::Point_3 &pt, const std::tuple<K::Point_3, int> &projPt,
  const V3 &ax_x, const V3 &ax_y, const V3 primitiveCenter[3], int eid)
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

  // PGO_ALOG(sinAngle >= 0);
  if (sinAngle < 0) {
    // std::cout << "sin < 0; angle=" << tod(angle) << std::endl;
    real angle1 = real(2 * M_PI) - angle;
    if (angle1 > real(M_PI)) {
    }
    else {
      angle = angle1;
    }

    // std::cout << "sin < 0; corrected angle=" << tod(angle) << std::endl;
  }

  V3 diff = primitiveCenter[(eid + 1) % 3] - primitiveCenter[eid];
  V3 diffAbs(abs(diff[0]), abs(diff[1]), abs(diff[2]));
  int maxIDX;
  diffAbs.maxCoeff(&maxIDX);

  V3 diff1 = origin - primitiveCenter[eid];
  real t = diff1[maxIDX] / diff[maxIDX];

  return V3(angle, t, r);
};

void prismVertexProjection(Poly &finalMesh, const pgo::Mesh::TriMeshGeo &targetMesh, const Poly &targetMeshPoly, const FaceGraphTree &targetMeshBVTreeExact,
  K::Triangle_3 &skeletonTri, double maxAllowedDist,
  const std::unordered_map<int, int> &vtxIDNew2Old, const std::map<int, std::array<int, 2>> &vertexIsOnTargetMesh,
  const std::vector<std::pair<int, int>> &selectedEdges)
{
  pgo::Mesh::TriMeshNeighbor targetMeshNeighbor(targetMesh);
  for (auto it = finalMesh.vertices_begin(); it != finalMesh.vertices_end(); ++it) {
    K::Point_3 pEK = it->point();
    auto projPtRet = computePrismProjectionPt(pEK, skeletonTri);
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
            maxAngle = CGAL::max(a, maxAngle);
          }

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
            int edgeID = std::get<1>(projPtRet);
            if (edgeID >= 0 && edgeID < 3) {
              K::Vector_3 n = skeletonTri[(edgeID + 1) % 3] - skeletonTri[edgeID];
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
                maxAngle = CGAL::max(a, maxAngle);
              }
            }
          }

          if (maxAngle < 20) {
            std::cout << "tgt vtx: Max angle too small a=" << tod(maxAngle) << std::endl;
            if (curLength > len) {
              it->point() = CGAL::midpoint(cpt, pEK);
            }
          }
          else {
            std::cout << "tgt edge: Max angle ok a=" << tod(maxAngle) << std::endl;
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

    // ES::V3d zz(-0.181634, 0.151783, 0.00689945);
    // ES::V3d zz1 = toVec3(it->point());
    // if (std::sqrt(dot(zz1 - zz, zz1 - zz)) < 1e-4) {
    //   std::cout << "catched" << std::endl;
    //   std::cout << it->id() << std::endl;
    // }
  }
}
}  // namespace MedialAxisRepresentation

int MedialAxisRepresentation::corefinePrismMeshWithTarget(const pgo::Mesh::TriMeshGeo &prismMesh, const ES::V3d centers[2],
  const pgo::Mesh::TriMeshGeo &targetMesh, const pgo::Mesh::TriMeshBVTree &targetMeshBVTree, const pgo::Mesh::TriMeshPseudoNormal &targetMeshNormals,
  double maxConfidentDist, double maxAllowedRadiusDist, pgo::Mesh::TriMeshGeo &meshOut, const char *filename)
{
  // double maxConfidentDist = 1e-3;
  // double maxAllowedRadiusDist = 0.1;
  constexpr int dumpMesh = 0;
  constexpr int removeSphereMeshTopology = 1;
  constexpr int triangulationOnly = 0;

  K::Point_3 primitiveCenterEK[3] = { toP3_EK(centers[0]), toP3_EK(centers[1]), toP3_EK(centers[2]) };
  V3 primitiveCenterIK[3] = { toV3(centers[0]), toV3(centers[1]), toV3(centers[2]) };
  K::Triangle_3 skeletonTri(primitiveCenterEK[0], primitiveCenterEK[1], primitiveCenterEK[2]);
  K::Segment_3 skeletonSegs[] = {
    K::Segment_3(primitiveCenterEK[0], primitiveCenterEK[1]),
    K::Segment_3(primitiveCenterEK[1], primitiveCenterEK[2]),
    K::Segment_3(primitiveCenterEK[2], primitiveCenterEK[0]),
  };
  K::Vector_3 skeletonEdgeDiff[] = {
    primitiveCenterEK[1] - primitiveCenterEK[0],
    primitiveCenterEK[2] - primitiveCenterEK[1],
    primitiveCenterEK[0] - primitiveCenterEK[2],
  };

  HClockPt t1 = HClock::now();

  pgo::Mesh::TriMeshBVTree prismMeshBVTree;
  prismMeshBVTree.buildByInertiaPartition(prismMesh);
  pgo::Mesh::TriMeshPseudoNormal prismMeshNormal;
  prismMeshNormal.buildPseudoNormals(prismMesh);
  std::vector<int> targetMeshChoosenVertex(targetMesh.numVertices(), 0);

  // compute the selected vertices for the target mesh to project
  tbb::parallel_for(0, targetMesh.numVertices(), [&](int vi) {
    // for (int vi = 0; vi < targetMesh.numVertices(); vi++) {
    // project target vertex to the skeleton
    // if the projection is successful, then keep going
    // otherwise, it is done
    ES::V3d p = targetMesh.pos(vi);
    auto projPt = computePrismProjectionPt(toP3_EK(p), skeletonTri);
    if (std::get<1>(projPt) < 0) {
      targetMeshChoosenVertex[vi] = 0;
      // continue;
      return;
    }

    thread_local std::vector<std::tuple<double, double, int>> closestDistDueryStack;
    thread_local std::stack<int> lineSegmentQueryStack;

    // compute distance to the cylinder mesh
    // close to the sphere
    closestDistDueryStack.clear();
    auto ret = prismMeshBVTree.closestTriangleQuery(prismMesh, p, closestDistDueryStack);
    if (ret.dist2 < maxConfidentDist * maxConfidentDist) {
      targetMeshChoosenVertex[vi] = 1;
    }

    // in contact with the sphere
    ES::V3d dir = (ret.closestPosition - p).normalized();
    ES::V3d n = prismMeshNormal.getPseudoNormal(prismMesh.triangles().data(), ret.triID, ret.feature);
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
      return;
      // continue;
    }

    while (!lineSegmentQueryStack.empty())
      lineSegmentQueryStack.pop();

    // intersect with sphere mesh
    retID = prismMeshBVTree.lineSegmentFirstIntersectionPointSemiExact(prismMesh, segStart, segEnd, segW, nullptr, &lineSegmentQueryStack);

    // if there is a intersection point
    if (retID >= 0) {
      ES::V3d pt = segStart * segW[0] + segEnd * segW[1];
      double dist = (pt - p).norm();

      // if dist is small enough
      if (dist <= maxAllowedRadiusDist * 1.2) {
        targetMeshChoosenVertex[vi] = 1;
      }
    }

    while (!lineSegmentQueryStack.empty())
      lineSegmentQueryStack.pop();

    segEnd = (p - segStart) * 5.0 + segStart;
    retID = prismMeshBVTree.lineSegmentFirstIntersectionPointSemiExact(prismMesh, segStart, segEnd, segW, nullptr, &lineSegmentQueryStack);

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
  });

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

        auto projPt0 = computePrismProjectionPt(toP3_EK(p0), skeletonTri);
        auto projPt1 = computePrismProjectionPt(toP3_EK(p1), skeletonTri);
        if (std::get<1>(projPt0) >= 0 && std::get<1>(projPt1) < 0) {
          targetMeshChoosenVertex[targetMesh.triVtxID(ti, (j + 1) % 3)] = 1;
        }
      }
      else if (targetMeshChoosenVertex[targetMesh.triVtxID(ti, j)] == 0 &&
        targetMeshChoosenVertex[targetMesh.triVtxID(ti, (j + 1) % 3)]) {
        ES::V3d p0 = targetMesh.pos(ti, j);
        ES::V3d p1 = targetMesh.pos(ti, (j + 1) % 3);

        auto projPt0 = computePrismProjectionPt(toP3_EK(p0), skeletonTri);
        auto projPt1 = computePrismProjectionPt(toP3_EK(p1), skeletonTri);
        if (std::get<1>(projPt0) < 0 && std::get<1>(projPt1) >= 0) {
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
  // x is always the triangle normal
  // y points outwards of the triangle
  // z is the edge dir
  K::Vector_3 nEK = CGAL::cross_product(primitiveCenterEK[1] - primitiveCenterEK[0], primitiveCenterEK[2] - primitiveCenterEK[0]);
  double nLength = std::sqrt(CGAL::to_double(nEK.squared_length()));
  V3 skeletonTriangleN = (primitiveCenterIK[1] - primitiveCenterIK[0]).cross(primitiveCenterIK[2] - primitiveCenterIK[0]);
  skeletonTriangleN /= sqrt(skeletonTriangleN.squaredNorm());

  V3 ax_x[3] = { skeletonTriangleN, skeletonTriangleN, skeletonTriangleN };
  V3 ax_z[3] = {
    primitiveCenterIK[1] - primitiveCenterIK[0],
    primitiveCenterIK[2] - primitiveCenterIK[1],
    primitiveCenterIK[0] - primitiveCenterIK[2],
  };

  ax_z[0] /= sqrt(ax_z[0].squaredNorm());
  ax_z[1] /= sqrt(ax_z[1].squaredNorm());
  ax_z[2] /= sqrt(ax_z[2].squaredNorm());

  V3 ax_y[3] = {
    ax_z[0].cross(ax_x[0]),
    ax_z[1].cross(ax_x[1]),
    ax_z[2].cross(ax_x[2]),
  };

  ax_y[0] /= sqrt(ax_y[0].squaredNorm());
  ax_y[1] /= sqrt(ax_y[1].squaredNorm());
  ax_y[2] /= sqrt(ax_y[2].squaredNorm());

  std::vector<K::Point_3> prismMeshVerticesEK(prismMesh.numVertices());
  for (int vi = 0; vi < prismMesh.numVertices(); vi++) {
    prismMeshVerticesEK[vi] = K::Point_3(prismMesh.pos(vi)[0], prismMesh.pos(vi)[1], prismMesh.pos(vi)[2]);
  }

  // analyze prism mesh ray type
  enum class PrismVertexRegionType : int
  {
    VERTEX0 = 0,
    VERTEX1,
    VERTEX2,
    CYLINDER0,
    CYLINDER1,
    CYLINDER2,
    TRIANGLE_POS,
    TRIANGLE_NEG,
    BORDER_CT0,
    BORDER_CT1,
    BORDER_CT2,
  };

  double ik_eps = 1e-5;

  std::vector<int> prismMeshVertexType(prismMesh.numVertices(), -1);
  for (int vi = 0; vi < (int)prismMeshVertexType.size(); vi++) {
    auto projPt = computePrismProjectionPt(prismMeshVerticesEK[vi], skeletonTri);
    if (std::get<1>(projPt) < 0) {
      int minID = -1;
      for (int j = 0; j < 3; j++) {
        if (std::get<0>(projPt) == primitiveCenterEK[j]) {
          minID = j;
          break;
        }
      }
      PGO_ALOG(minID >= 0);

      K::Vector_3 diff = prismMeshVerticesEK[vi] - std::get<0>(projPt);
      K::Vector_3 e0 = primitiveCenterEK[(minID + 1) % 3] - primitiveCenterEK[minID];
      K::Vector_3 e1 = primitiveCenterEK[(minID + 2) % 3] - primitiveCenterEK[minID];

      double dotValue[2] = {
        tod(CGAL::scalar_product(diff, e0)) / (std::sqrt(tod(diff.squared_length()) * tod(e0.squared_length()))),
        tod(CGAL::scalar_product(diff, e1)) / (std::sqrt(tod(diff.squared_length()) * tod(e1.squared_length()))),
      };

      if (std::abs(dotValue[0]) < ik_eps && std::abs(std::abs(dotValue[1])) < ik_eps) {
        // keep components that are in the direction of the normal
        K::FT c = CGAL::scalar_product(nEK, diff) / nEK.squared_length();
        diff = nEK * c;
        prismMeshVerticesEK[vi] = std::get<0>(projPt) + diff;

        projPt = computePrismProjectionPt(prismMeshVerticesEK[vi], skeletonTri);
        PGO_ALOG(std::get<1>(projPt) >= 0);
        K::Vector_3 newDiff = prismMeshVerticesEK[vi] - std::get<0>(projPt);

        PGO_ALOG(CGAL::scalar_product(newDiff, skeletonEdgeDiff[0]) == 0 &&
          CGAL::scalar_product(newDiff, skeletonEdgeDiff[1]) == 0);

        if (std::get<0>(projPt) == primitiveCenterEK[0]) {
          prismMeshVertexType[vi] = (int)PrismVertexRegionType::VERTEX0;
        }
        else if (std::get<0>(projPt) == primitiveCenterEK[1]) {
          prismMeshVertexType[vi] = (int)PrismVertexRegionType::VERTEX1;
        }
        else if (std::get<0>(projPt) == primitiveCenterEK[2]) {
          prismMeshVertexType[vi] = (int)PrismVertexRegionType::VERTEX2;
        }
        else {
          throw std::runtime_error("Impossible");
        }
      }
      else if (std::abs(dotValue[0]) < std::abs(std::abs(dotValue[1]))) {
        K::Vector_3 diff2 = prismMeshVerticesEK[vi] - primitiveCenterEK[minID];
        K::FT t = CGAL::scalar_product(diff2, e0) / e0.squared_length();
        K::Point_3 edgeProjPos = primitiveCenterEK[minID] + e0 * t;
        K::Vector_3 dir = prismMeshVerticesEK[vi] - edgeProjPos;
        PGO_ALOG(CGAL::scalar_product(dir, e0) == 0);

        K::Vector_3 err = edgeProjPos - std::get<0>(projPt);
        double errDist = std::sqrt(tod(err.squared_length()));
        PGO_ALOG(errDist < ik_eps);

        prismMeshVerticesEK[vi] = std::get<0>(projPt) + dir;

        int edgeID = minID;
        prismMeshVertexType[vi] = (int)PrismVertexRegionType::CYLINDER0 + edgeID;
      }
      else if (std::abs(dotValue[1]) < std::abs(std::abs(dotValue[0]))) {
        K::Vector_3 diff2 = prismMeshVerticesEK[vi] - primitiveCenterEK[minID];
        K::FT t = CGAL::scalar_product(diff2, e1) / e1.squared_length();
        K::Point_3 edgeProjPos = primitiveCenterEK[minID] + e1 * t;
        K::Vector_3 dir = prismMeshVerticesEK[vi] - edgeProjPos;
        PGO_ALOG(CGAL::scalar_product(dir, e1) == 0);

        K::Vector_3 err = edgeProjPos - std::get<0>(projPt);
        double errDist = std::sqrt(tod(err.squared_length()));
        PGO_ALOG(errDist < ik_eps);

        prismMeshVerticesEK[vi] = std::get<0>(projPt) + dir;

        int edgeID = (minID + 2) % 3;
        prismMeshVertexType[vi] = (int)PrismVertexRegionType::CYLINDER0 + edgeID;
      }
      else {
        throw std::runtime_error("impossible");
      }
    }
    else if (std::get<1>(projPt) >= 0 && std::get<1>(projPt) < 3) {
      int eid = std::get<1>(projPt);
      V3 cc = computeCylinderCoordinate(prismMeshVerticesEK[vi], projPt,
        ax_x[eid], ax_y[eid], primitiveCenterIK, eid);
      double angle_d = tod(cc[0]);

      K::FT d2[3] = {
        CGAL::squared_distance(primitiveCenterEK[0], std::get<0>(projPt)),
        CGAL::squared_distance(primitiveCenterEK[1], std::get<0>(projPt)),
        CGAL::squared_distance(primitiveCenterEK[2], std::get<0>(projPt))
      };

      int minID = -1;
      K::FT minDist = 1e100;
      for (int i = 0; i < 3; i++) {
        if (minDist > d2[i]) {
          minID = i;
          minDist = d2[i];
        }
      }

      PGO_ALOG(minID >= 0 && minID < 3);
      double d = std::sqrt(CGAL::to_double(minDist));
      K::Vector_3 diff = prismMeshVerticesEK[vi] - std::get<0>(projPt);
      ES::V3d ray0 = toVec3(diff);
      ES::V3d ray1 = toVec3(nEK);

      ray0 /= (ray0).norm();
      ray1 /= (ray1).norm();

      double cosAngle = ray0.dot(ray1);

      if (d < ik_eps && std::abs(std::abs(cosAngle) - 1) < ik_eps * 10) {
        K::Point_3 origin;
        origin = primitiveCenterEK[minID];
        prismMeshVertexType[vi] = (int)PrismVertexRegionType::VERTEX0 + minID;

        K::FT c = CGAL::scalar_product(nEK, diff) / nEK.squared_length();
        diff = nEK * c;
        prismMeshVerticesEK[vi] = origin + diff;
        projPt = computePrismProjectionPt(prismMeshVerticesEK[vi], skeletonTri);
        PGO_ALOG(std::get<1>(projPt) >= 0 && std::get<1>(projPt) < 3);
      }
      else if (std::abs(angle_d) < ik_eps * 10 || std::abs(angle_d - M_PI) < ik_eps * 10 ||
        std::abs(angle_d - 2 * M_PI) < ik_eps * 10) {
        K::FT c = CGAL::scalar_product(nEK, diff) / nEK.squared_length();
        diff = nEK * c;

        prismMeshVerticesEK[vi] = std::get<0>(projPt) + diff;
        projPt = computePrismProjectionPt(prismMeshVerticesEK[vi], skeletonTri);
        PGO_ALOG(std::get<1>(projPt) >= 0 && std::get<1>(projPt) < 3);

        K::Vector_3 newDiff = prismMeshVerticesEK[vi] - std::get<0>(projPt);
        PGO_ALOG(CGAL::scalar_product(newDiff, skeletonEdgeDiff[0]) == 0 &&
          CGAL::scalar_product(newDiff, skeletonEdgeDiff[1]) == 0);

        if (std::get<0>(projPt) == primitiveCenterEK[0]) {
          prismMeshVertexType[vi] = (int)PrismVertexRegionType::VERTEX0;
        }
        else if (std::get<0>(projPt) == primitiveCenterEK[1]) {
          prismMeshVertexType[vi] = (int)PrismVertexRegionType::VERTEX1;
        }
        else if (std::get<0>(projPt) == primitiveCenterEK[2]) {
          prismMeshVertexType[vi] = (int)PrismVertexRegionType::VERTEX2;
        }
        else {
          prismMeshVertexType[vi] = (int)PrismVertexRegionType::BORDER_CT0 + eid;
        }
      }
      else {
        prismMeshVertexType[vi] = (int)PrismVertexRegionType::CYLINDER0 + eid;
      }
    }
    else {
      K::FT d2[6] = {
        CGAL::squared_distance(primitiveCenterEK[0], std::get<0>(projPt)),
        CGAL::squared_distance(primitiveCenterEK[1], std::get<0>(projPt)),
        CGAL::squared_distance(primitiveCenterEK[2], std::get<0>(projPt)),
        CGAL::squared_distance(skeletonSegs[0], std::get<0>(projPt)),
        CGAL::squared_distance(skeletonSegs[1], std::get<0>(projPt)),
        CGAL::squared_distance(skeletonSegs[2], std::get<0>(projPt)),
      };

      double d_d[6];
      for (int j = 0; j < 6; j++) {
        d_d[j] = std::sqrt(CGAL::to_double(d2[j]));
      }

      // std::cout << CGAL::to_double(d2[0]) << ','
      //           << CGAL::to_double(d2[1]) << ','
      //           << CGAL::to_double(d2[2]) << ','
      //           << CGAL::to_double(d2[3]) << ','
      //           << CGAL::to_double(d2[4]) << ','
      //           << CGAL::to_double(d2[5]) << std::endl;

      int minID = -1;
      K::FT minDist = 1e100;
      for (int i = 0; i < 6; i++) {
        if (minDist > d2[i]) {
          minID = i;
          minDist = d2[i];
        }
      }

      PGO_ALOG(minID >= 0 && minID < 6);
      double d = std::sqrt(CGAL::to_double(minDist));
      if (d_d[0] < ik_eps && d_d[3] < ik_eps && d_d[5] < ik_eps) {
        K::Vector_3 diff = prismMeshVerticesEK[vi] - std::get<0>(projPt);
        prismMeshVertexType[vi] = (int)PrismVertexRegionType::VERTEX0;
        prismMeshVerticesEK[vi] = primitiveCenterEK[0] + diff;

        projPt = computePrismProjectionPt(prismMeshVerticesEK[vi], skeletonTri);
        PGO_ALOG(std::get<1>(projPt) >= 0 && std::get<1>(projPt) < 3);
      }
      else if (d_d[1] < ik_eps && d_d[3] < ik_eps && d_d[4] < ik_eps) {
        K::Vector_3 diff = prismMeshVerticesEK[vi] - std::get<0>(projPt);
        prismMeshVertexType[vi] = (int)PrismVertexRegionType::VERTEX1;
        prismMeshVerticesEK[vi] = primitiveCenterEK[1] + diff;

        projPt = computePrismProjectionPt(prismMeshVerticesEK[vi], skeletonTri);
        PGO_ALOG(std::get<1>(projPt) >= 0 && std::get<1>(projPt) < 3);
      }
      else if (d_d[2] < ik_eps && d_d[4] < ik_eps && d_d[5] < ik_eps) {
        K::Vector_3 diff = prismMeshVerticesEK[vi] - std::get<0>(projPt);
        prismMeshVertexType[vi] = (int)PrismVertexRegionType::VERTEX2;
        prismMeshVerticesEK[vi] = primitiveCenterEK[2] + diff;

        projPt = computePrismProjectionPt(prismMeshVerticesEK[vi], skeletonTri);
        PGO_ALOG(std::get<1>(projPt) >= 0 && std::get<1>(projPt) < 3);
      }
      else if (d < ik_eps) {
        K::Vector_3 diff = prismMeshVerticesEK[vi] - std::get<0>(projPt);
        K::Point_3 origin;
        if (minID >= 0 && minID < 3) {
          origin = primitiveCenterEK[minID];
          prismMeshVertexType[vi] = (int)PrismVertexRegionType::VERTEX0 + minID;
        }
        else {
          int vids[2] = { minID - 3, (minID - 3 + 1) % 3 };
          K::Vector_3 d1 = std::get<0>(projPt) - primitiveCenterEK[vids[0]];
          K::Vector_3 d0 = primitiveCenterEK[vids[1]] - primitiveCenterEK[vids[0]];
          K::FT t = CGAL::scalar_product(d0, d1) / d0.squared_length();
          origin = primitiveCenterEK[vids[0]] + t * d0;

          prismMeshVertexType[vi] = (int)PrismVertexRegionType::BORDER_CT0 + vids[0];
        }

        prismMeshVerticesEK[vi] = origin + diff;
        projPt = computePrismProjectionPt(prismMeshVerticesEK[vi], skeletonTri);
        PGO_ALOG(std::get<1>(projPt) >= 0 && std::get<1>(projPt) < 3);
      }
      else {
        K::Vector_3 diff = prismMeshVerticesEK[vi] - std::get<0>(projPt);
        if (CGAL::scalar_product(diff, nEK) > 0) {
          prismMeshVertexType[vi] = (int)PrismVertexRegionType::TRIANGLE_POS;
        }
        else {
          prismMeshVertexType[vi] = (int)PrismVertexRegionType::TRIANGLE_NEG;
        }
      }
    }
  }

  int flagCounters[3] = { 0, 0, 0 };
  for (int vi = 0; vi < (int)prismMeshVertexType.size(); vi++) {
    if (prismMeshVertexType[vi] == (int)PrismVertexRegionType::VERTEX0)
      flagCounters[0] += 1;
    else if (prismMeshVertexType[vi] == (int)PrismVertexRegionType::VERTEX1)
      flagCounters[1] += 1;
    else if (prismMeshVertexType[vi] == (int)PrismVertexRegionType::VERTEX2)
      flagCounters[2] += 1;
  }
  PGO_ALOG(flagCounters[0] == 2);
  PGO_ALOG(flagCounters[1] == 2);
  PGO_ALOG(flagCounters[2] == 2);

  auto isRegionBorderVertex = [&prismMeshVertexType](int vi) {
    if (prismMeshVertexType[vi] == (int)PrismVertexRegionType::VERTEX0 ||
      prismMeshVertexType[vi] == (int)PrismVertexRegionType::VERTEX1 ||
      prismMeshVertexType[vi] == (int)PrismVertexRegionType::VERTEX2 ||
      prismMeshVertexType[vi] == (int)PrismVertexRegionType::BORDER_CT0 ||
      prismMeshVertexType[vi] == (int)PrismVertexRegionType::BORDER_CT1 ||
      prismMeshVertexType[vi] == (int)PrismVertexRegionType::BORDER_CT2) {
      return true;
    }
    else {
      return false;
    }
  };

  if (dumpMesh) {
    pgo::Mesh::TriMeshGeo m;
    for (int i = 0; i < (int)prismMeshVertexType.size(); i++) {
      if (isRegionBorderVertex(i)) {
        m.addPos(toVec3(prismMeshVerticesEK[i]));
      }
    }
    m.save("zza.obj");
  }

  // build cgal poly
  Poly prismMeshPoly;
  pgo::CGALInterface::triangleMesh2Polyhedron(prismMesh, prismMeshPoly);
  for (auto vit = prismMeshPoly.vertices_begin(); vit != prismMeshPoly.vertices_end(); vit++) {
    vit->point() = prismMeshVerticesEK[vit->id()];
  }

  // build aabb tree
  FaceGraphTree prismMeshBVTreeExact(CGAL::faces(prismMeshPoly).first, CGAL::faces(prismMeshPoly).second, prismMeshPoly);

  struct PrismTriangleInfo
  {
    K::Triangle_2 triangle2D;
    std::array<real, 3> vertexRadii;
    PrismVertexRegionType regionType;
  };

  std::vector<PrismTriangleInfo> prismMeshProjectedTriangles(prismMesh.numTriangles());
  K::FT K_pi(M_PI);
  K::FT K_zero(0);

  // direct computation of the border vertex positions
  auto pt2D_is_on_border = [&](const K::Point_2 &pt, int triID, K::Point_3 &p3D_EK) -> bool {
    const ES::V3i &vertexID = prismMesh.tri(triID);
    PrismVertexRegionType rType = prismMeshProjectedTriangles[triID].regionType;
    const auto &tri = prismMeshProjectedTriangles[triID].triangle2D;
    K::FT w[3];
    // std::cout << CGAL::to_double(tri[0][0]) << ',' << CGAL::to_double(tri[0][1]) << std::endl;
    // std::cout << CGAL::to_double(tri[1][0]) << ',' << CGAL::to_double(tri[1][1]) << std::endl;
    // std::cout << CGAL::to_double(tri[2][0]) << ',' << CGAL::to_double(tri[2][1]) << std::endl;
    // std::cout << CGAL::to_double(pt[0]) << ',' << CGAL::to_double(pt[1]) << std::endl;
    CGAL::Barycentric_coordinates::triangle_coordinates_2(tri[0], tri[1], tri[2], pt, w);
    // std::cout << "w:" << CGAL::to_double(w[0]) << ',' << CGAL::to_double(w[1]) << ',' << CGAL::to_double(w[2]) << std::endl;

    // if this triangle has a border edge
    // and pt is on this border edge
    for (int i = 0; i < 3; i++) {
      if (isRegionBorderVertex(vertexID[i]) && pt == tri[i]) {
        p3D_EK = prismMeshVerticesEK[vertexID[i]];
        return true;
      }
    }

    for (int i = 0; i < 3; i++) {
      for (int j = i + 1; j < 3; j++) {
        if (isRegionBorderVertex(vertexID[i]) && isRegionBorderVertex(vertexID[j])) {
          int k = 3 - i - j;
          PGO_ALOG(isRegionBorderVertex(vertexID[k]) == false);
          if (w[k] == 0) {
            // when in this case, angle = 0 or pi, so it doesn't matter any more
            // we can directly linear interpolate the vertex

            p3D_EK = K::Point_3(
              prismMeshVerticesEK[vertexID[i]][0] * w[i] + prismMeshVerticesEK[vertexID[j]][0] * w[j],
              prismMeshVerticesEK[vertexID[i]][1] * w[i] + prismMeshVerticesEK[vertexID[j]][1] * w[j],
              prismMeshVerticesEK[vertexID[i]][2] * w[i] + prismMeshVerticesEK[vertexID[j]][2] * w[j]);
            // p3D_EK = primitiveCenterEK[cylinderID] + skeletonEdgeDiff[cylinderID] * pt[1];

            // std::cout << "3D:" << toVec3(p3D_EK) << std::endl;

            return true;
          }
        }
      }
    }

    return false;
  };

  // compute 3D positions
  auto pt2D_to_pt3D = [&](const K::Point_2 &pt, int triID) -> IK::Point_3 {
    PrismVertexRegionType rType = prismMeshProjectedTriangles[triID].regionType;
    int edgeID = 0;
    if (rType == PrismVertexRegionType::CYLINDER0 ||
      rType == PrismVertexRegionType::CYLINDER1 ||
      rType == PrismVertexRegionType::CYLINDER2) {
      edgeID = ((int)rType - (int)PrismVertexRegionType::CYLINDER0);
    }
    else {
      throw std::invalid_argument("Impossible");
    }

    const auto &tri = prismMeshProjectedTriangles[triID].triangle2D;
    K::FT w[3];
    CGAL::Barycentric_coordinates::triangle_coordinates_2(tri[0], tri[1], tri[2], pt, w);
    /*
    std::cout << "angle: " << tod(pt[0]) << ',' << tod(pt[1]) << std::endl;
    std::cout << tod(w[0]) << ',' << tod(w[1]) << ',' << tod(w[2]) << std::endl;
    */

    IK::FT r = 0;
    for (int j = 0; j < 3; j++) {
      r += prismMeshProjectedTriangles[triID].vertexRadii[j] * EK_to_IK<real>(w[j]);
    }

    V3 x;
    real angle = EK_to_IK<real>(pt[0]);
    x = (ax_x[edgeID] * cos(angle) + ax_y[edgeID] * sin(angle)) * r;

    V3 vert = (primitiveCenterIK[(edgeID + 1) % 3] - primitiveCenterIK[edgeID]) * EK_to_IK<real>(pt[1]) + primitiveCenterIK[edgeID];
    x += vert;

    return IK::Point_3(x[0], x[1], x[2]);
  };

  for (int ti = 0; ti < prismMesh.numTriangles(); ti++) {
    // tbb::parallel_for(0, prismMesh.numTriangles(), [&](int ti) {
    int edgeID = -1;
    for (int j = 0; j < 3; j++) {
      if (prismMeshVertexType[prismMesh.tri(ti)[j]] == (int)PrismVertexRegionType::CYLINDER0 ||
        prismMeshVertexType[prismMesh.tri(ti)[j]] == (int)PrismVertexRegionType::CYLINDER1 ||
        prismMeshVertexType[prismMesh.tri(ti)[j]] == (int)PrismVertexRegionType::CYLINDER2) {
        edgeID = (int)prismMeshVertexType[prismMesh.tri(ti)[j]] - (int)PrismVertexRegionType::CYLINDER0;
        break;
      }
    }

    if (edgeID >= 0) {
      for (int j = 0; j < 3; j++) {
        PGO_ALOG(prismMeshVertexType[prismMesh.tri(ti)[j]] == (int)PrismVertexRegionType::CYLINDER0 + edgeID ||
          isRegionBorderVertex(prismMesh.tri(ti)[j]));
      }

      K::Point_2 p2d[3];
      for (int j = 0; j < 3; j++) {
        int vid = prismMesh.tri(ti)[j];
        const K::Point_3 &pt = prismMeshVerticesEK[vid];
        auto projPt = computePrismProjectionPt(pt, skeletonTri);
        V3 uvw = computeCylinderCoordinate(pt, projPt, ax_x[edgeID], ax_y[edgeID], primitiveCenterIK, edgeID);

        // std::cout << "vtype: " << prismMeshVertexType[vid] << "; is border: " << isRegionBorderVertex(vid) << std::endl;
        if (isRegionBorderVertex(vid)) {
          PGO_ALOG((edgeID == 0 &&
                     (prismMeshVertexType[vid] == (int)PrismVertexRegionType::BORDER_CT0 ||
                       prismMeshVertexType[vid] == (int)PrismVertexRegionType::VERTEX0 ||
                       prismMeshVertexType[vid] == (int)PrismVertexRegionType::VERTEX1)) ||
            (edgeID == 1 &&
              (prismMeshVertexType[vid] == (int)PrismVertexRegionType::BORDER_CT1 ||
                prismMeshVertexType[vid] == (int)PrismVertexRegionType::VERTEX1 ||
                prismMeshVertexType[vid] == (int)PrismVertexRegionType::VERTEX2)) ||
            (edgeID == 2 &&
              (prismMeshVertexType[vid] == (int)PrismVertexRegionType::BORDER_CT2 ||
                prismMeshVertexType[vid] == (int)PrismVertexRegionType::VERTEX2 ||
                prismMeshVertexType[vid] == (int)PrismVertexRegionType::VERTEX0)));

          K::Vector_3 diff = pt - primitiveCenterEK[edgeID];
          K::FT t = CGAL::scalar_product(skeletonEdgeDiff[edgeID], diff) / skeletonEdgeDiff[edgeID].squared_length();
          // std::cout << "t: " << tod(t) << ',' << tod(uvw[1]) << std::endl;

          if (std::abs(tod(uvw[0])) < 1e-5) {
            p2d[j] = K::Point_2(0, t);
          }
          else if (std::abs(tod(uvw[0]) - M_PI) < 1e-5) {
            p2d[j] = K::Point_2(M_PI, t);
          }
          else {
            throw std::runtime_error("impossible");
          }
        }
        else {
          p2d[j] = K::Point_2(IK_to_EK(uvw[0]), IK_to_EK(uvw[1]));
        }

        prismMeshProjectedTriangles[ti].vertexRadii[j] = uvw[2];
      }

      prismMeshProjectedTriangles[ti].triangle2D = K::Triangle_2(p2d[0], p2d[1], p2d[2]);
      prismMeshProjectedTriangles[ti].regionType = (PrismVertexRegionType)((int)PrismVertexRegionType::CYLINDER0 + edgeID);
    }
    else {
      K::Vector_3 e01 = prismMeshVerticesEK[prismMesh.tri(ti)[1]] - prismMeshVerticesEK[prismMesh.tri(ti)[0]];
      K::Vector_3 e02 = prismMeshVerticesEK[prismMesh.tri(ti)[2]] - prismMeshVerticesEK[prismMesh.tri(ti)[0]];
      K::Vector_3 n = CGAL::cross_product(e01, e02);

      if (CGAL::scalar_product(n, nEK) > 0) {
        prismMeshProjectedTriangles[ti].regionType = PrismVertexRegionType::TRIANGLE_POS;
      }
      else {
        prismMeshProjectedTriangles[ti].regionType = PrismVertexRegionType::TRIANGLE_NEG;
      }
    }
  }  //);

  if constexpr (dumpMesh) {
    pgo::Mesh::TriMeshGeo mm;
    for (int ti = 0; ti < (int)prismMeshProjectedTriangles.size(); ti++) {
      if (prismMeshProjectedTriangles[ti].regionType == PrismVertexRegionType::CYLINDER0 ||
        prismMeshProjectedTriangles[ti].regionType == PrismVertexRegionType::CYLINDER1 ||
        prismMeshProjectedTriangles[ti].regionType == PrismVertexRegionType::CYLINDER2) {
        IK::Point_3 p[3];
        for (int j = 0; j < 3; j++) {
          K::Point_3 pp;
          if (pt2D_is_on_border(prismMeshProjectedTriangles[ti].triangle2D[j], ti, pp)) {
            p[j] = toP3_IK(pp);
          }
          else {
            p[j] = pt2D_to_pt3D(prismMeshProjectedTriangles[ti].triangle2D[j], ti);
          }
        }

        ES::V3d pV3d[3] = {
          ES::V3d(tod(p[0][0]), tod(p[0][1]), tod(p[0][2])),
          ES::V3d(tod(p[1][0]), tod(p[1][1]), tod(p[1][2])),
          ES::V3d(tod(p[2][0]), tod(p[2][1]), tod(p[2][2])),
        };

        mm.addMesh(pgo::Mesh::createSingleTriangleMesh(pV3d[0], pV3d[1], pV3d[2]));
        // std::cout << ti << ' ';
      }
      else if (prismMeshProjectedTriangles[ti].regionType == PrismVertexRegionType::TRIANGLE_POS) {
      }
      else if (prismMeshProjectedTriangles[ti].regionType == PrismVertexRegionType::TRIANGLE_NEG) {
      }
      else {
        throw std::runtime_error("aaaaa");
      }
    }
    std::cout << std::endl;

    mm.save("rarrr.obj");
  }

  // construct prism for the center
  std::vector<K::Point_3> prismCenterPrismVertex(6);
  prismCenterPrismVertex[0] = skeletonTri[0] + nEK / nLength;
  prismCenterPrismVertex[1] = skeletonTri[1] + nEK / nLength;
  prismCenterPrismVertex[2] = skeletonTri[2] + nEK / nLength;
  prismCenterPrismVertex[3] = skeletonTri[0] - nEK / nLength;
  prismCenterPrismVertex[4] = skeletonTri[1] - nEK / nLength;
  prismCenterPrismVertex[5] = skeletonTri[2] - nEK / nLength;

  Poly prismCenterPoly;
  CGAL::convex_hull_3(prismCenterPrismVertex.begin(), prismCenterPrismVertex.end(), prismCenterPoly);
  FaceGraphTree prismCenterPrismBVTree(CGAL::faces(prismCenterPoly).first, CGAL::faces(prismCenterPoly).second, prismCenterPoly);
  CGAL::Side_of_triangle_mesh<Poly, K> prismCenterInsideChecker(prismCenterPoly);

  struct EdgeSegments
  {
    // for cylinder region
    std::array<K::Segment_3, 3> cylinderSegPieces;
    std::array<std::array<int, 2>, 3> segPieceEndIsTargetVertex;
    std::array<int, 3> cylinderSegPieceDoesExist;
    std::array<K::Segment_2, 3> cylinderProjections;
    // for triangle region
    std::array<K::Triangle_3, 2> cuttingTriangles;
    Poly cuttingMesh;
    std::array<K::Segment_3, 2> cuttingMeshTargetVertexSegments;
    std::array<int, 2> cuttingMeshIsTargetVertexSegments;
  };

  // build cutting mesh
  std::vector<EdgeSegments> cuttingSegments(selectedEdges.size());
  std::vector<std::vector<int>> edgeIntersectedTriangles(selectedEdges.size());

  auto triangleRegionCutting = [&](const K::Segment_3 &seg, int ei,
                                 const int segEndIsTargetVertex[2],
                                 std::vector<Poly::Facet_handle> &faceHandles) {
    // the cutting will only happens on the triangle
    // it forms a rectangle
    auto proj0 = computePrismProjectionPt(seg[0], skeletonTri);
    auto proj1 = computePrismProjectionPt(seg[1], skeletonTri);

    K::Vector_3 diff0 = seg[0] - std::get<0>(proj0);
    K::Vector_3 diff1 = seg[1] - std::get<0>(proj1);
    PGO_ALOG(CGAL::scalar_product(diff0, skeletonEdgeDiff[0]) == 0 &&
      CGAL::scalar_product(diff0, skeletonEdgeDiff[1]) == 0);

    PGO_ALOG(CGAL::scalar_product(diff1, skeletonEdgeDiff[0]) == 0 &&
      CGAL::scalar_product(diff1, skeletonEdgeDiff[1]) == 0);

    if (CGAL::scalar_product(diff0, diff1) < 0) {
      std::cout << "???? " << ei << std::endl;
      throw std::runtime_error("penetrate skeleton");
    }

    K::Point_3 v0, v1, v2, v3;
    v0 = std::get<0>(proj0);
    v1 = std::get<0>(proj1);

    if (CGAL::scalar_product(diff0, nEK) > 0) {
      v2 = std::get<0>(proj1) + nEK / nLength * 2.0;
      v3 = std::get<0>(proj0) + nEK / nLength * 2.0;
    }
    else {
      v2 = std::get<0>(proj1) - nEK / nLength * 2.0;
      v3 = std::get<0>(proj0) - nEK / nLength * 2.0;
    }

    // std::cout << "len: " << std::sqrt(tod(diff0.squared_length())) << ',' << std::sqrt(tod(diff1.squared_length())) << std::endl;
    // std::cout << "dot: " << tod(CGAL::scalar_product(diff0, nEK)) << ',' << tod(CGAL::scalar_product(diff1, nEK)) << std::endl;

    K::Triangle_3 tri0(v0, v1, v2);
    K::Triangle_3 tri1(v0, v2, v3);

    faceHandles.clear();
    prismMeshBVTreeExact.all_intersected_primitives(tri0, std::back_inserter(faceHandles));
    prismMeshBVTreeExact.all_intersected_primitives(tri1, std::back_inserter(faceHandles));

    for (auto fh : faceHandles) {
      int ti = (int)fh->id();
      edgeIntersectedTriangles[ei].emplace_back(ti);
    }

    cuttingSegments[ei].cuttingTriangles[0] = tri0;
    cuttingSegments[ei].cuttingTriangles[1] = tri1;

    cuttingSegments[ei].cuttingMesh.make_triangle(v0, v1, v2);

    Poly::Halfedge_handle v0h;
    for (auto vit = cuttingSegments[ei].cuttingMesh.vertices_begin(); vit != cuttingSegments[ei].cuttingMesh.vertices_end(); ++vit) {
      if (vit->point() == v0) {
        v0h = vit->halfedge();
        PGO_ALOG(v0h->opposite()->vertex()->point() == v2);
        break;
      }
    }

    Poly::Halfedge_handle g = v0h->opposite();
    Poly::Halfedge_handle h = g->prev();
    h = cuttingSegments[ei].cuttingMesh.add_vertex_and_facet_to_border(h, g);
    h->vertex()->point() = v3;

    cuttingSegments[ei].cuttingMeshIsTargetVertexSegments[0] = segEndIsTargetVertex[0];
    cuttingSegments[ei].cuttingMeshIsTargetVertexSegments[1] = segEndIsTargetVertex[1];

    if (segEndIsTargetVertex[0]) {
      cuttingSegments[ei].cuttingMeshTargetVertexSegments[0] = K::Segment_3(v0, v3);

      // std::cout << ei << ',' << 0 << ": seg on: ";
      // std::cout << cuttingSegments[ei].cuttingMeshTargetVertexSegments[0].has_on(targetMeshVerticesER[selectedEdges[ei].first]) << ','
      //           << cuttingSegments[ei].cuttingMeshTargetVertexSegments[0].has_on(targetMeshVerticesER[selectedEdges[ei].second]) << std::endl;

      // std::cout << ei << ',' << 0 << ": line on: ";
      // std::cout << cuttingSegments[ei].cuttingMeshTargetVertexSegments[0].supporting_line().has_on(targetMeshVerticesER[selectedEdges[ei].first]) << ','
      //           << cuttingSegments[ei].cuttingMeshTargetVertexSegments[0].supporting_line().has_on(targetMeshVerticesER[selectedEdges[ei].second]) << std::endl;
    }

    if (segEndIsTargetVertex[1]) {
      cuttingSegments[ei].cuttingMeshTargetVertexSegments[1] = K::Segment_3(v1, v2);

      // std::cout << ei << ',' << 1 << ": seg on: ";
      // std::cout << cuttingSegments[ei].cuttingMeshTargetVertexSegments[1].has_on(targetMeshVerticesER[selectedEdges[ei].first]) << ','
      //           << cuttingSegments[ei].cuttingMeshTargetVertexSegments[1].has_on(targetMeshVerticesER[selectedEdges[ei].second]) << std::endl;

      // std::cout << ei << ',' << 1 << ": line on: ";
      // std::cout << cuttingSegments[ei].cuttingMeshTargetVertexSegments[1].supporting_line().has_on(targetMeshVerticesER[selectedEdges[ei].first]) << ','
      //           << cuttingSegments[ei].cuttingMeshTargetVertexSegments[1].supporting_line().has_on(targetMeshVerticesER[selectedEdges[ei].second]) << std::endl;
    }
  };

  K::Plane_3 cylinderCuttingPlane[3];
  K::Vector_3 ax_y_nonunitEK[3];
  for (int i = 0; i < 3; i++) {
    ax_y_nonunitEK[i] = CGAL::cross_product(skeletonEdgeDiff[i], nEK);
    // n (x - x0) = 0
    cylinderCuttingPlane[i] = K::Plane_3(primitiveCenterEK[i], ax_y_nonunitEK[i]);
    // std::cout << cylinderCuttingPlane[i].has_on_positive_side(primitiveCenterEK[i] + ax_y_nonunitEK[i]) << std::endl;
  }

  auto edgeRegionCutting = [&](const K::Segment_3 &seg, int ei) {
    // std::cout << ei << std::endl;
    // compute cylinder coordinate
    int hasBorder[3] = { 0, 0, 0 };
    K::Point_3 borderPoint3D[3];
    K::Point_2 borderPoint2D[3];
    for (int j = 0; j < 3; j++) {
      cuttingSegments[ei].cylinderSegPieceDoesExist[j] = 0;
      auto ret = CGAL::intersection(seg, cylinderCuttingPlane[j]);
      int hasBorderVtx = -1;
      if (ret) {
        if (const K::Segment_3 *s = boost::get<K::Segment_3>(&*ret)) {
          throw std::runtime_error("impossible");
        }
        else if (const K::Point_3 *p = boost::get<K::Point_3>(&*ret)) {
          if (cylinderCuttingPlane[j].has_on_positive_side(seg[0])) {
            cuttingSegments[ei].cylinderSegPieces[j] = K::Segment_3(seg[0], *p);
            cuttingSegments[ei].segPieceEndIsTargetVertex[j][0] = 1;
            cuttingSegments[ei].segPieceEndIsTargetVertex[j][1] = 0;

            cuttingSegments[ei].cylinderSegPieceDoesExist[j] = 1;

            hasBorderVtx = 1;
          }
          else if (cylinderCuttingPlane[j].has_on_positive_side(seg[1])) {
            cuttingSegments[ei].cylinderSegPieces[j] = K::Segment_3(*p, seg[1]);
            cuttingSegments[ei].segPieceEndIsTargetVertex[j][0] = 0;
            cuttingSegments[ei].segPieceEndIsTargetVertex[j][1] = 1;

            cuttingSegments[ei].cylinderSegPieceDoesExist[j] = 1;
            hasBorderVtx = 0;
          }
          else {
            throw std::runtime_error("impossible");
          }
        }
      }
      else {
        if (cylinderCuttingPlane[j].has_on_positive_side(seg[0]) &&
          cylinderCuttingPlane[j].has_on_positive_side(seg[1])) {
          cuttingSegments[ei].cylinderSegPieces[j] = seg;
          cuttingSegments[ei].cylinderSegPieceDoesExist[j] = 1;
          cuttingSegments[ei].segPieceEndIsTargetVertex[j][0] = 1;
          cuttingSegments[ei].segPieceEndIsTargetVertex[j][1] = 1;
        }
      }

      if (cuttingSegments[ei].cylinderSegPieceDoesExist[j]) {
        auto outSegProj0 = computeCylinderProjectionPt(cuttingSegments[ei].cylinderSegPieces[j][0], primitiveCenterEK, j);
        auto outSegProj1 = computeCylinderProjectionPt(cuttingSegments[ei].cylinderSegPieces[j][1], primitiveCenterEK, j);
        V3 uvw0 = computeCylinderCoordinate(cuttingSegments[ei].cylinderSegPieces[j][0], outSegProj0, ax_x[j], ax_y[j], primitiveCenterIK, j);
        V3 uvw1 = computeCylinderCoordinate(cuttingSegments[ei].cylinderSegPieces[j][1], outSegProj1, ax_x[j], ax_y[j], primitiveCenterIK, j);
        K::Point_2 segEnds[2] = {
          K::Point_2(IK_to_EK(uvw0[0]), IK_to_EK(uvw0[1])),
          K::Point_2(IK_to_EK(uvw1[0]), IK_to_EK(uvw1[1]))
        };

        if (hasBorderVtx >= 0) {
          K::Vector_3 diff = cuttingSegments[ei].cylinderSegPieces[j][hasBorderVtx] - primitiveCenterEK[j];
          K::FT t = CGAL::scalar_product(diff, skeletonEdgeDiff[j]) / skeletonEdgeDiff[j].squared_length();
          // std::cout << "t: " << tod(t) << ',' << tod(segEnds[hasBorderVtx][1]) << std::endl;
          // std::cout << "a: " << tod(segEnds[hasBorderVtx][0]) << std::endl;

          if (std::abs(tod(segEnds[hasBorderVtx][0])) < 1e-5) {
            segEnds[hasBorderVtx] = K::Point_2(0, t);
          }
          else if (std::abs(tod(segEnds[hasBorderVtx][0]) - M_PI) < 1e-5) {
            segEnds[hasBorderVtx] = K::Point_2(M_PI, t);
          }
          else {
            throw std::runtime_error("impossible");
          }

          borderPoint3D[j] = cuttingSegments[ei].cylinderSegPieces[j][hasBorderVtx];
          borderPoint2D[j] = segEnds[hasBorderVtx];
          hasBorder[j] = 1;
        }
        else {
          hasBorder[j] = 0;
        }

        cuttingSegments[ei].cylinderProjections[j] = K::Segment_2(segEnds[0], segEnds[1]);
        // std::cout << j << ":(" << tod(segEnds[0][0]) << ',' << tod(segEnds[0][1]) << "),(" << tod(segEnds[1][0]) << ',' << tod(segEnds[1][1]) << ")" << std::endl;
      }
    }

    for (int ti = 0; ti < (int)prismMeshProjectedTriangles.size(); ti++) {
      // if (ti == 1123) {
      //   std::cout << "aa\n";

      //   const auto &tri = prismMeshProjectedTriangles[ti].triangle2D;
      //   std::cout << 0 << ":" << tod(tri[0][0]) << ',' << tod(tri[0][1]) << std::endl;
      //   std::cout << 1 << ":" << tod(tri[1][0]) << ',' << tod(tri[1][1]) << std::endl;
      //   std::cout << 2 << ":" << tod(tri[2][0]) << ',' << tod(tri[2][1]) << std::endl;
      // }

      if (prismMeshProjectedTriangles[ti].regionType == PrismVertexRegionType::CYLINDER0) {
        if (cuttingSegments[ei].cylinderSegPieceDoesExist[0] &&
          CGAL::do_intersect(prismMeshProjectedTriangles[ti].triangle2D, cuttingSegments[ei].cylinderProjections[0])) {
          /*
        int eid = 0;
        if (hasBorder[eid] && borderPoint2D[eid][1] >= 0 && borderPoint2D[eid][1] <= 1) {
          std::cout << tod(cuttingSegments[ei].cylinderProjections[eid][0][0]) << ',' << tod(cuttingSegments[ei].cylinderProjections[eid][0][1]) << std::endl;
          std::cout << tod(cuttingSegments[ei].cylinderProjections[eid][1][0]) << ',' << tod(cuttingSegments[ei].cylinderProjections[eid][1][1]) << std::endl;
          std::cout << "seg0:" << toVec3(borderPoint3D[eid]) << std::endl;

          const auto &tri = prismMeshProjectedTriangles[ti].triangle2D;
          K::FT w[3];
          CGAL::Barycentric_coordinates::triangle_coordinates_2(tri[0], tri[1], tri[2], borderPoint2D[eid], w);

          if (w[0] >= 0 && w[0] <= 1 &&
            w[1] >= 0 && w[1] <= 1 &&
            w[2] >= 0 && w[2] <= 1) {
            int vids[3] = {
              prismMesh.triVtxID(ti, 0),
              prismMesh.triVtxID(ti, 1),
              prismMesh.triVtxID(ti, 2),
            };

            K::FT finalPt[3] = {
              0, 0, 0
            };

            std::cout << std::setprecision(17);
            std::cout << tod(w[0]) << ',' << tod(w[1]) << ',' << tod(w[2]) << std::endl;

            K::FT tri1Vtx[3];
            for (int i = 0; i < 3; i++) {
              K::Vector_3 diff = prismMeshVerticesEK[vids[i]] - primitiveCenterEK[eid];
              K::FT t = CGAL::scalar_product(skeletonEdgeDiff[eid], diff) / skeletonEdgeDiff[eid].squared_length();
              tri1Vtx[i] = tri[i][1];
            }

            K::Segment_3 borderSeg;

            if (w[0] == 0) {
              // t * tri[1] + (1-t) * tri[2] = s
              // t (tri[1] - tri[2]) + tri[2] = s
              K::FT a = borderPoint2D[eid][1] - tri1Vtx[2];
              K::FT b = tri1Vtx[1] - tri1Vtx[2];
              w[1] = a / b;
              w[2] = 1 - w[1];

              borderSeg = K::Segment_3(prismMeshVerticesEK[vids[1]], prismMeshVerticesEK[vids[2]]);
            }
            else if (w[1] == 0) {
              K::FT a = borderPoint2D[eid][1] - tri1Vtx[2];
              K::FT b = tri1Vtx[0] - tri1Vtx[2];
              w[0] = a / b;
              w[2] = 1 - w[0];

              borderSeg = K::Segment_3(prismMeshVerticesEK[vids[0]], prismMeshVerticesEK[vids[2]]);
            }
            else if (w[2] == 0) {
              K::FT a = borderPoint2D[eid][1] - tri1Vtx[1];
              K::FT b = tri1Vtx[0] - tri1Vtx[1];
              w[0] = a / b;
              w[1] = 1 - w[0];

              borderSeg = K::Segment_3(prismMeshVerticesEK[vids[0]], prismMeshVerticesEK[vids[1]]);
            }

            for (int r = 0; r < 3; r++) {
              finalPt[0] += prismMeshVerticesEK[vids[r]][0] * w[r];
              finalPt[1] += prismMeshVerticesEK[vids[r]][1] * w[r];
              finalPt[2] += prismMeshVerticesEK[vids[r]][2] * w[r];
            }

            auto projPt = computePrismProjectionPt(borderPoint3D[eid], skeletonTri);
            K::Vector_3 dir = borderPoint3D[eid] - std::get<0>(projPt);
            K::Segment_3 verticalSeg(borderPoint3D[eid] + dir, std::get<0>(projPt));
            if (borderSeg.is_degenerate()) {
              throw std::runtime_error("rr");
            }

            if (verticalSeg.is_degenerate()) {
              throw std::runtime_error("rr");
            }

            auto ret = CGAL::intersection(borderSeg, verticalSeg);
            K::Point_3 tgtp;
            if (ret) {
              if (const K::Point_3 *p = boost::get<K::Point_3>(&*ret)) {
                tgtp = *p;
              }
              else {
                throw std::runtime_error("1");
              }
            }
            else {
              std::cout << std::sqrt(tod(CGAL::squared_distance(borderSeg, verticalSeg))) << std::endl;
              throw std::runtime_error("1");
            }

            std::cout << "tri:\n";
            std::cout << tod(tri[0][0]) << ',' << tod(tri[0][1]) << std::endl;
            std::cout << tod(tri[1][0]) << ',' << tod(tri[1][1]) << std::endl;
            std::cout << tod(tri[2][0]) << ',' << tod(tri[2][1]) << std::endl;

            std::cout << tod(tri1Vtx[0]) << std::endl;
            std::cout << tod(tri1Vtx[1]) << std::endl;
            std::cout << tod(tri1Vtx[2]) << std::endl;

            K::Point_3 pt1(finalPt[0], finalPt[1], finalPt[2]);

            std::cout << "edge:\n";
            std::cout << tod(borderPoint2D[eid][0]) << ',' << tod(borderPoint2D[eid][1]) << std::endl;
            std::cout << tod(w[0]) << ',' << tod(w[1]) << ',' << tod(w[2]) << std::endl;
            if (pt1 == tgtp) {
            }
            else {
              throw std::runtime_error("1");
            }
            std::cout << "p0:" << toVec3(tgtp) << std::endl;
            std::cout << "p1:" << toVec3(pt1) << std::endl;
            std::cout << tod((pt1 - tgtp).squared_length()) << std::endl;
          }
        }*/
          edgeIntersectedTriangles[ei].emplace_back(ti);
        }
      }
      else if (prismMeshProjectedTriangles[ti].regionType == PrismVertexRegionType::CYLINDER1) {
        if (cuttingSegments[ei].cylinderSegPieceDoesExist[1] &&
          CGAL::do_intersect(prismMeshProjectedTriangles[ti].triangle2D, cuttingSegments[ei].cylinderProjections[1])) {
          /*
        int eid = 1;
        if (hasBorder[eid] && borderPoint2D[eid][1] >= 0 && borderPoint2D[eid][1] <= 1) {
          std::cout << tod(cuttingSegments[ei].cylinderProjections[eid][0][0]) << ',' << tod(cuttingSegments[ei].cylinderProjections[eid][0][1]) << std::endl;
          std::cout << tod(cuttingSegments[ei].cylinderProjections[eid][1][0]) << ',' << tod(cuttingSegments[ei].cylinderProjections[eid][1][1]) << std::endl;
          std::cout << "seg0:" << toVec3(borderPoint3D[eid]) << std::endl;

          const auto &tri = prismMeshProjectedTriangles[ti].triangle2D;
          K::FT w[3];
          CGAL::Barycentric_coordinates::triangle_coordinates_2(tri[0], tri[1], tri[2], borderPoint2D[eid], w);

          if (w[0] >= 0 && w[0] <= 1 &&
            w[1] >= 0 && w[1] <= 1 &&
            w[2] >= 0 && w[2] <= 1) {
            int vids[3] = {
              prismMesh.triVtxID(ti, 0),
              prismMesh.triVtxID(ti, 1),
              prismMesh.triVtxID(ti, 2),
            };

            K::FT finalPt[3] = {
              0, 0, 0
            };

            std::cout << std::setprecision(17);
            std::cout << tod(w[0]) << ',' << tod(w[1]) << ',' << tod(w[2]) << std::endl;

            K::FT tri1Vtx[3];
            for (int i = 0; i < 3; i++) {
              K::Vector_3 diff = prismMeshVerticesEK[vids[i]] - primitiveCenterEK[eid];
              K::FT t = CGAL::scalar_product(skeletonEdgeDiff[eid], diff) / skeletonEdgeDiff[eid].squared_length();
              tri1Vtx[i] = tri[i][1];
            }

            K::Segment_3 borderSeg;

            if (w[0] == 0) {
              // t * tri[1] + (1-t) * tri[2] = s
              // t (tri[1] - tri[2]) + tri[2] = s
              K::FT a = borderPoint2D[eid][1] - tri1Vtx[2];
              K::FT b = tri1Vtx[1] - tri1Vtx[2];
              w[1] = a / b;
              w[2] = 1 - w[1];

              borderSeg = K::Segment_3(prismMeshVerticesEK[vids[1]], prismMeshVerticesEK[vids[2]]);
            }
            else if (w[1] == 0) {
              K::FT a = borderPoint2D[eid][1] - tri1Vtx[2];
              K::FT b = tri1Vtx[0] - tri1Vtx[2];
              w[0] = a / b;
              w[2] = 1 - w[0];

              borderSeg = K::Segment_3(prismMeshVerticesEK[vids[0]], prismMeshVerticesEK[vids[2]]);
            }
            else if (w[2] == 0) {
              K::FT a = borderPoint2D[eid][1] - tri1Vtx[1];
              K::FT b = tri1Vtx[0] - tri1Vtx[1];
              w[0] = a / b;
              w[1] = 1 - w[0];

              borderSeg = K::Segment_3(prismMeshVerticesEK[vids[0]], prismMeshVerticesEK[vids[1]]);
            }

            for (int r = 0; r < 3; r++) {
              finalPt[0] += prismMeshVerticesEK[vids[r]][0] * w[r];
              finalPt[1] += prismMeshVerticesEK[vids[r]][1] * w[r];
              finalPt[2] += prismMeshVerticesEK[vids[r]][2] * w[r];
            }

            auto projPt = computePrismProjectionPt(borderPoint3D[eid], skeletonTri);
            K::Vector_3 dir = borderPoint3D[eid] - std::get<0>(projPt);
            K::Segment_3 verticalSeg(borderPoint3D[eid] + dir, std::get<0>(projPt));
            if (borderSeg.is_degenerate()) {
              throw std::runtime_error("rr");
            }

            if (verticalSeg.is_degenerate()) {
              throw std::runtime_error("rr");
            }

            auto ret = CGAL::intersection(borderSeg, verticalSeg);
            K::Point_3 tgtp;
            if (ret) {
              if (const K::Point_3 *p = boost::get<K::Point_3>(&*ret)) {
                tgtp = *p;
              }
              else {
                throw std::runtime_error("1");
              }
            }
            else {
              throw std::runtime_error("1");
            }

            std::cout << "tri:\n";
            std::cout << tod(tri[0][0]) << ',' << tod(tri[0][1]) << std::endl;
            std::cout << tod(tri[1][0]) << ',' << tod(tri[1][1]) << std::endl;
            std::cout << tod(tri[2][0]) << ',' << tod(tri[2][1]) << std::endl;

            std::cout << tod(tri1Vtx[0]) << std::endl;
            std::cout << tod(tri1Vtx[1]) << std::endl;
            std::cout << tod(tri1Vtx[2]) << std::endl;

            K::Point_3 pt1(finalPt[0], finalPt[1], finalPt[2]);

            std::cout << "edge:\n";
            std::cout << tod(borderPoint2D[eid][0]) << ',' << tod(borderPoint2D[eid][1]) << std::endl;
            std::cout << tod(w[0]) << ',' << tod(w[1]) << ',' << tod(w[2]) << std::endl;
            if (pt1 == tgtp) {
            }
            else {
              throw std::runtime_error("1");
            }
            std::cout << "p0:" << toVec3(tgtp) << std::endl;
            std::cout << "p1:" << toVec3(pt1) << std::endl;
            std::cout << tod((pt1 - tgtp).squared_length()) << std::endl;
          }
        }*/

          edgeIntersectedTriangles[ei].emplace_back(ti);
        }
      }
      else if (prismMeshProjectedTriangles[ti].regionType == PrismVertexRegionType::CYLINDER2) {
        if (cuttingSegments[ei].cylinderSegPieceDoesExist[2] &&
          CGAL::do_intersect(prismMeshProjectedTriangles[ti].triangle2D, cuttingSegments[ei].cylinderProjections[2])) {
          /*
        int eid = 2;
        if (hasBorder[eid] && borderPoint2D[eid][1] >= 0 && borderPoint2D[eid][1] <= 1) {
          std::cout << tod(cuttingSegments[ei].cylinderProjections[eid][0][0]) << ',' << tod(cuttingSegments[ei].cylinderProjections[eid][0][1]) << std::endl;
          std::cout << tod(cuttingSegments[ei].cylinderProjections[eid][1][0]) << ',' << tod(cuttingSegments[ei].cylinderProjections[eid][1][1]) << std::endl;
          std::cout << "seg0:" << toVec3(borderPoint3D[eid]) << std::endl;

          const auto &tri = prismMeshProjectedTriangles[ti].triangle2D;
          K::FT w[3];
          CGAL::Barycentric_coordinates::triangle_coordinates_2(tri[0], tri[1], tri[2], borderPoint2D[eid], w);

          if (w[0] >= 0 && w[0] <= 1 &&
            w[1] >= 0 && w[1] <= 1 &&
            w[2] >= 0 && w[2] <= 1) {
            int vids[3] = {
              prismMesh.triVtxID(ti, 0),
              prismMesh.triVtxID(ti, 1),
              prismMesh.triVtxID(ti, 2),
            };

            K::FT finalPt[3] = {
              0, 0, 0
            };

            std::cout << std::setprecision(17);
            std::cout << tod(w[0]) << ',' << tod(w[1]) << ',' << tod(w[2]) << std::endl;

            K::FT tri1Vtx[3];
            for (int i = 0; i < 3; i++) {
              K::Vector_3 diff = prismMeshVerticesEK[vids[i]] - primitiveCenterEK[eid];
              K::FT t = CGAL::scalar_product(skeletonEdgeDiff[eid], diff) / skeletonEdgeDiff[eid].squared_length();
              tri1Vtx[i] = tri[i][1];
            }

            K::Segment_3 borderSeg;

            if (w[0] == 0) {
              // t * tri[1] + (1-t) * tri[2] = s
              // t (tri[1] - tri[2]) + tri[2] = s
              K::FT a = borderPoint2D[eid][1] - tri1Vtx[2];
              K::FT b = tri1Vtx[1] - tri1Vtx[2];
              w[1] = a / b;
              w[2] = 1 - w[1];

              borderSeg = K::Segment_3(prismMeshVerticesEK[vids[1]], prismMeshVerticesEK[vids[2]]);
            }
            else if (w[1] == 0) {
              K::FT a = borderPoint2D[eid][1] - tri1Vtx[2];
              K::FT b = tri1Vtx[0] - tri1Vtx[2];
              w[0] = a / b;
              w[2] = 1 - w[0];

              borderSeg = K::Segment_3(prismMeshVerticesEK[vids[0]], prismMeshVerticesEK[vids[2]]);
            }
            else if (w[2] == 0) {
              K::FT a = borderPoint2D[eid][1] - tri1Vtx[1];
              K::FT b = tri1Vtx[0] - tri1Vtx[1];
              w[0] = a / b;
              w[1] = 1 - w[0];

              borderSeg = K::Segment_3(prismMeshVerticesEK[vids[0]], prismMeshVerticesEK[vids[1]]);
            }

            for (int r = 0; r < 3; r++) {
              finalPt[0] += prismMeshVerticesEK[vids[r]][0] * w[r];
              finalPt[1] += prismMeshVerticesEK[vids[r]][1] * w[r];
              finalPt[2] += prismMeshVerticesEK[vids[r]][2] * w[r];
            }

            auto projPt = computePrismProjectionPt(borderPoint3D[eid], skeletonTri);
            K::Vector_3 dir = borderPoint3D[eid] - std::get<0>(projPt);
            K::Segment_3 verticalSeg(borderPoint3D[eid] + dir, std::get<0>(projPt));
            if (borderSeg.is_degenerate()) {
              throw std::runtime_error("rr");
            }

            if (verticalSeg.is_degenerate()) {
              throw std::runtime_error("rr");
            }

            auto ret = CGAL::intersection(borderSeg, verticalSeg);
            K::Point_3 tgtp;
            if (ret) {
              if (const K::Point_3 *p = boost::get<K::Point_3>(&*ret)) {
                tgtp = *p;
              }
              else {
                throw std::runtime_error("1");
              }
            }
            else {
              throw std::runtime_error("1");
            }

            std::cout << "tri:\n";
            std::cout << tod(tri[0][0]) << ',' << tod(tri[0][1]) << std::endl;
            std::cout << tod(tri[1][0]) << ',' << tod(tri[1][1]) << std::endl;
            std::cout << tod(tri[2][0]) << ',' << tod(tri[2][1]) << std::endl;

            std::cout << tod(tri1Vtx[0]) << std::endl;
            std::cout << tod(tri1Vtx[1]) << std::endl;
            std::cout << tod(tri1Vtx[2]) << std::endl;

            K::Point_3 pt1(finalPt[0], finalPt[1], finalPt[2]);

            std::cout << "edge:\n";
            std::cout << tod(borderPoint2D[eid][0]) << ',' << tod(borderPoint2D[eid][1]) << std::endl;
            std::cout << tod(w[0]) << ',' << tod(w[1]) << ',' << tod(w[2]) << std::endl;
            if (pt1 == tgtp) {
            }
            else {
              throw std::runtime_error("1");
            }
            std::cout << "p0:" << toVec3(tgtp) << std::endl;
            std::cout << "p1:" << toVec3(pt1) << std::endl;
            std::cout << tod((pt1 - tgtp).squared_length()) << std::endl;
          }
        }*/

          edgeIntersectedTriangles[ei].emplace_back(ti);
        }
      }
    }
  };

  // for each edge find all intersected triangles
  tbb::parallel_for(0, (int)selectedEdges.size(), [&](int ei) {
    // for (int ei = 0; ei < (int)selectedEdges.size(); ei++) {
    int vids[2] = {
      selectedEdges[ei].first,
      selectedEdges[ei].second
    };

    // pgo::Mesh::TriMeshGeo rr;
    // rr.addMesh(createSingleTriangleMesh(targetMesh.pos(vids[0]), targetMesh.pos(vids[1]), targetMesh.pos(vids[0]) + Vec3d(1e-5)));
    // rr.save("a.obj");

    K::Segment_3 seg(targetMeshVerticesER[vids[0]], targetMeshVerticesER[vids[1]]);
    CGAL::Bounded_side ret0 = prismCenterInsideChecker(seg[0]);
    CGAL::Bounded_side ret1 = prismCenterInsideChecker(seg[1]);

    thread_local std::vector<Poly::Facet_handle> faceHandles;
    // if both are inside
    if (ret0 == CGAL::ON_BOUNDED_SIDE && ret1 == CGAL::ON_BOUNDED_SIDE) {
      int isTargetVertex[2] = { 1, 1 };
      triangleRegionCutting(seg, ei, isTargetVertex, faceHandles);
      pgo::BasicAlgorithms::sortAndDeduplicateWithErase(edgeIntersectedTriangles[ei]);
    }
    // if one point is inside the other is not
    else if ((ret0 == CGAL::ON_BOUNDED_SIDE && ret1 != CGAL::ON_BOUNDED_SIDE) ||
      (ret0 != CGAL::ON_BOUNDED_SIDE && ret1 == CGAL::ON_BOUNDED_SIDE)) {
      // the cutting will happens on
      faceHandles.clear();
      prismCenterPrismBVTree.all_intersected_primitives(seg, std::back_inserter(faceHandles));
      PGO_ALOG(faceHandles.size() == 1ull);

      K::Triangle_3 prismTri(
        faceHandles[0]->halfedge()->vertex()->point(),
        faceHandles[0]->halfedge()->next()->vertex()->point(),
        faceHandles[0]->halfedge()->next()->next()->vertex()->point());
      auto ret = CGAL::intersection(prismTri, seg);

      K::Segment_3 cutSegs[2];
      int isTargetVertex[2] = { 1, 1 };

      if (ret) {
        if (const K::Segment_3 *s = boost::get<K::Segment_3>(&*ret)) {
          throw std::runtime_error("impossible");
        }
        else if (const K::Point_3 *p = boost::get<K::Point_3>(&*ret)) {
          if (ret0 == CGAL::ON_BOUNDED_SIDE) {
            cutSegs[0] = K::Segment_3(seg[0], *p);
            isTargetVertex[0] = 1;
            isTargetVertex[1] = 0;

            cutSegs[1] = K::Segment_3(*p, seg[1]);
          }
          else if (ret1 == CGAL::ON_BOUNDED_SIDE) {
            cutSegs[0] = K::Segment_3(*p, seg[1]);
            isTargetVertex[0] = 0;
            isTargetVertex[1] = 1;

            cutSegs[1] = K::Segment_3(seg[0], *p);
          }
          else {
            throw std::runtime_error("impossible");
          }
        }
      }

      triangleRegionCutting(cutSegs[0], ei, isTargetVertex, faceHandles);
      edgeRegionCutting(seg, ei);

      pgo::BasicAlgorithms::sortAndDeduplicateWithErase(edgeIntersectedTriangles[ei]);
    }
    else {
      faceHandles.clear();
      prismCenterPrismBVTree.all_intersected_primitives(seg, std::back_inserter(faceHandles));
      if (faceHandles.size()) {
        // at most there are two intersected points or one intersected segments
        K::Segment_3 intersectedSeg;
        K::Point_3 intersectedPt[2];

        int numFoundIntersectedSeg = 0;
        int numFoundIntersectedPt = 0;

        for (auto fh : faceHandles) {
          K::Triangle_3 prismTri(
            fh->halfedge()->vertex()->point(),
            fh->halfedge()->next()->vertex()->point(),
            fh->halfedge()->next()->next()->vertex()->point());

          // std::cout << toVec3(prismTri[0]) << ',' << toVec3(prismTri[1]) << ',' << toVec3(prismTri[2]) << std::endl;
          // std::cout << toVec3(seg[0]) << ',' << toVec3(seg[1]) << std::endl;

          auto ret = CGAL::intersection(prismTri, seg);
          if (ret) {
            if (const K::Segment_3 *s = boost::get<K::Segment_3>(&*ret)) {
              if (numFoundIntersectedSeg == 0) {
                intersectedSeg = *s;
                numFoundIntersectedSeg++;
              }
              else {
                PGO_ALOG(intersectedSeg.supporting_line().has_on((*s)[0]));
                PGO_ALOG(intersectedSeg.supporting_line().has_on((*s)[1]));

                if ((*s)[0] == intersectedSeg[0]) {
                  intersectedSeg = K::Segment_3((*s)[1], intersectedSeg[1]);
                }
                else if ((*s)[0] == intersectedSeg[1]) {
                  intersectedSeg = K::Segment_3((*s)[1], intersectedSeg[0]);
                }
                else if ((*s)[1] == intersectedSeg[0]) {
                  intersectedSeg = K::Segment_3((*s)[0], intersectedSeg[1]);
                }
                else if ((*s)[1] == intersectedSeg[1]) {
                  intersectedSeg = K::Segment_3((*s)[0], intersectedSeg[0]);
                }
                else {
                  throw std::runtime_error("impossible");
                }
              }

              // this part tries to remove points on the segments
              int hasOn[2] = { 0, 0 };
              for (int i = 0; i < numFoundIntersectedPt; i++) {
                if (intersectedSeg.has_on(intersectedPt[i])) {
                  hasOn[i] = 1;
                }
              }

              int inc = 0;
              for (int i = 0; i < numFoundIntersectedPt; i++) {
                if (hasOn[i]) {
                  intersectedPt[inc++] = intersectedPt[i];
                }
              }
              numFoundIntersectedPt = inc;
            }
            else if (const K::Point_3 *p = boost::get<K::Point_3>(&*ret)) {
              if (numFoundIntersectedSeg) {
                if (intersectedSeg.has_on(*p) == false) {
                  bool found = false;
                  for (int i = 0; i < numFoundIntersectedPt; i++) {
                    if (intersectedPt[i] == *p) {
                      found = true;
                      break;
                    }
                  }

                  if (found == false) {
                    intersectedPt[numFoundIntersectedPt++] = *p;
                  }
                }
              }
              else {
                bool found = false;
                for (int i = 0; i < numFoundIntersectedPt; i++) {
                  if (intersectedPt[i] == *p) {
                    found = true;
                    break;
                  }
                }

                if (found == false) {
                  intersectedPt[numFoundIntersectedPt++] = *p;
                }
              }
            }
          }
        }

        if (numFoundIntersectedSeg) {
          PGO_ALOG(numFoundIntersectedPt == 0);
          // it means no triangle project
          edgeRegionCutting(seg, ei);
        }
        else {
          // this is the only case there is a triangle projection
          if (numFoundIntersectedPt == 2) {
            int isTargetVertex[2] = { 0, 0 };
            triangleRegionCutting(K::Segment_3(intersectedPt[0], intersectedPt[1]), ei, isTargetVertex, faceHandles);
          }

          edgeRegionCutting(seg, ei);
        }
        pgo::BasicAlgorithms::sortAndDeduplicateWithErase(edgeIntersectedTriangles[ei]);
      }
      else {
        edgeRegionCutting(seg, ei);
        pgo::BasicAlgorithms::sortAndDeduplicateWithErase(edgeIntersectedTriangles[ei]);
      }
    }
  });

  // for each intersected triangle, gather intersections
  std::vector<std::vector<int>> triangleIntersectedEdges(prismMesh.numTriangles());

  for (int ei = 0; ei < (int)selectedEdges.size(); ei++) {
    for (int i = 0; i < (int)edgeIntersectedTriangles[ei].size(); i++) {
      int ti = edgeIntersectedTriangles[ei][i];
      triangleIntersectedEdges[ti].emplace_back(ei);
    }
  }

  int maxEdgeCount = 0;
  int maxCountID = -1;
  for (int ti = 0; ti < prismMesh.numTriangles(); ti++) {
    if ((int)triangleIntersectedEdges[ti].size() > maxEdgeCount) {
      maxEdgeCount = (int)triangleIntersectedEdges[ti].size();
      maxCountID = ti;
    }
  }

  std::cout << "#edges max: " << maxEdgeCount << std::endl;
  std::cout << "id: " << maxCountID << std::endl;

  // =====================================
  // for each intersected triangle, compute corefined triangles
  // we use cutting mesh to cut each sphere triangle
  std::vector<Patch2D> corefinedPatches2D(prismMesh.numTriangles());
  std::vector<Patch3D> corefinedPatches3D(prismMesh.numTriangles());

  std::atomic<int> counter;

  tbb::parallel_for(0, prismMesh.numTriangles(), [&](int ti) {
    // for (int ti = 0; ti < prismMesh.numTriangles(); ti++) {
    PrismVertexRegionType rType = prismMeshProjectedTriangles[ti].regionType;
    if (rType == PrismVertexRegionType::TRIANGLE_POS ||
      rType == PrismVertexRegionType::TRIANGLE_NEG) {
      // for (int ti = 0; ti < sphereMesh.numTriangles(); ti++) {
      const K::Point_3 &v0 = prismMeshVerticesEK[prismMesh.triVtxID(ti, 0)];
      const K::Point_3 &v1 = prismMeshVerticesEK[prismMesh.triVtxID(ti, 1)];
      const K::Point_3 &v2 = prismMeshVerticesEK[prismMesh.triVtxID(ti, 2)];

      if (triangleIntersectedEdges[ti].size() == 0ull) {
        corefinedPatches3D[ti].vertices = { v0, v1, v2 };
        corefinedPatches3D[ti].triangles.emplace_back(std::array<int, 3>{ 0, 1, 2 });

        corefinedPatches3D[ti].edgeVertices[0] = { 0, 1 };
        corefinedPatches3D[ti].edgeVertices[1] = { 1, 2 };
        corefinedPatches3D[ti].edgeVertices[2] = { 2, 0 };

        corefinedPatches3D[ti].cornerVertices[0] = 0;
        corefinedPatches3D[ti].cornerVertices[1] = 1;
        corefinedPatches3D[ti].cornerVertices[2] = 2;

        corefinedPatches3D[ti].vertexIsOnTargetMesh.assign(3, std::array<int, 2>{ -1, -1 });
      }
      else {
        Poly tri;
        tri.make_triangle(v0, v1, v2);

        for (int ei = 0; ei < (int)triangleIntersectedEdges[ti].size(); ei++) {
          int edgeID = triangleIntersectedEdges[ti][ei];
          CGAL::Polygon_mesh_processing::corefine(tri, cuttingSegments[edgeID].cuttingMesh,
            CGAL::parameters::do_not_modify(false),
            CGAL::parameters::do_not_modify(true));
        }

        // std::ofstream("b.off") << tri;

        corefinedPatches3D[ti].vertices.resize(tri.size_of_vertices());
        corefinedPatches3D[ti].triangles.resize(tri.size_of_facets());
        int inc = 0;
        for (auto vit = tri.vertices_begin(); vit != tri.vertices_end(); ++vit) {
          vit->id() = (size_t)inc;
          corefinedPatches3D[ti].vertices[vit->id()] = vit->point();
          inc++;
        }

        inc = 0;
        for (auto fit = tri.facets_begin(); fit != tri.facets_end(); ++fit) {
          corefinedPatches3D[ti].triangles[inc][0] = fit->halfedge()->vertex()->id();
          corefinedPatches3D[ti].triangles[inc][1] = fit->halfedge()->next()->vertex()->id();
          corefinedPatches3D[ti].triangles[inc][2] = fit->halfedge()->next()->next()->vertex()->id();
          inc++;
        }

        corefinedPatches3D[ti].edgeVertices[0].reserve(corefinedPatches3D[ti].vertices.size());
        corefinedPatches3D[ti].edgeVertices[1].reserve(corefinedPatches3D[ti].vertices.size());
        corefinedPatches3D[ti].edgeVertices[2].reserve(corefinedPatches3D[ti].vertices.size());

        K::Segment_3 e01(v0, v1), e12(v1, v2), e20(v2, v0);
        for (int vi = 0; vi < (int)corefinedPatches3D[ti].vertices.size(); vi++) {
          if (e01.has_on(corefinedPatches3D[ti].vertices[vi])) {
            corefinedPatches3D[ti].edgeVertices[0].emplace_back(vi);
          }

          if (e12.has_on(corefinedPatches3D[ti].vertices[vi])) {
            corefinedPatches3D[ti].edgeVertices[1].emplace_back(vi);
          }

          if (e20.has_on(corefinedPatches3D[ti].vertices[vi])) {
            corefinedPatches3D[ti].edgeVertices[2].emplace_back(vi);
          }

          if (corefinedPatches3D[ti].vertices[vi] == v0) {
            corefinedPatches3D[ti].cornerVertices[0] = vi;
          }
          else if (corefinedPatches3D[ti].vertices[vi] == v1) {
            corefinedPatches3D[ti].cornerVertices[1] = vi;
          }
          else if (corefinedPatches3D[ti].vertices[vi] == v2) {
            corefinedPatches3D[ti].cornerVertices[2] = vi;
          }
        }

        corefinedPatches3D[ti].triangleEdgeIDs.reserve((int)corefinedPatches3D[ti].triangles.size());
        for (int tii = 0; tii < (int)corefinedPatches3D[ti].triangles.size(); tii++) {
          for (int j = 0; j < 3; j++) {
            const K::Point_3 &p0 = corefinedPatches3D[ti].vertices[corefinedPatches3D[ti].triangles[tii][j]];
            const K::Point_3 &p1 = corefinedPatches3D[ti].vertices[corefinedPatches3D[ti].triangles[tii][(j + 1) % 3]];

            for (int ei = 0; ei < (int)triangleIntersectedEdges[ti].size(); ei++) {
              int edgeID = triangleIntersectedEdges[ti][ei];
              const K::Triangle_3 &cuttingTri0 = cuttingSegments[edgeID].cuttingTriangles[0];
              const K::Triangle_3 &cuttingTri1 = cuttingSegments[edgeID].cuttingTriangles[1];

              if ((cuttingTri0.has_on(p0) || cuttingTri1.has_on(p0)) &&
                (cuttingTri0.has_on(p1) || cuttingTri1.has_on(p1))) {
                corefinedPatches3D[ti].triangleEdgeIDs.emplace_back(tii, j, triangleIntersectedEdges[ti][ei]);
                break;
              }
            }
          }
        }

        corefinedPatches3D[ti].vertexIsOnTargetMesh.assign(corefinedPatches3D[ti].vertices.size(), std::array<int, 2>{ -1, -1 });
        for (int vi = 0; vi < (int)corefinedPatches3D[ti].vertices.size(); vi++) {
          for (int ei = 0; ei < (int)triangleIntersectedEdges[ti].size(); ei++) {
            int edgeID = triangleIntersectedEdges[ti][ei];
            const K::Triangle_3 &cuttingTri0 = cuttingSegments[edgeID].cuttingTriangles[0];
            const K::Triangle_3 &cuttingTri1 = cuttingSegments[edgeID].cuttingTriangles[1];

            bool onSeg = false;
            if (cuttingSegments[edgeID].cuttingMeshIsTargetVertexSegments[0]) {
              if (cuttingSegments[edgeID].cuttingMeshTargetVertexSegments[0].has_on(corefinedPatches3D[ti].vertices[vi])) {
                int vid;
                if (cuttingSegments[edgeID].cuttingMeshTargetVertexSegments[0].has_on(targetMeshVerticesER[selectedEdges[edgeID].first])) {
                  vid = selectedEdges[edgeID].first;
                }
                else {
                  PGO_ALOG(cuttingSegments[edgeID].cuttingMeshTargetVertexSegments[0].has_on(targetMeshVerticesER[selectedEdges[edgeID].second]));
                  vid = selectedEdges[edgeID].second;
                }

                corefinedPatches3D[ti].vertexIsOnTargetMesh[vi][0] = 0;
                corefinedPatches3D[ti].vertexIsOnTargetMesh[vi][1] = vid;
                onSeg = true;
              }
            }

            if (cuttingSegments[edgeID].cuttingMeshIsTargetVertexSegments[1]) {
              if (cuttingSegments[edgeID].cuttingMeshTargetVertexSegments[1].has_on(corefinedPatches3D[ti].vertices[vi])) {
                int vid;
                if (cuttingSegments[edgeID].cuttingMeshTargetVertexSegments[1].has_on(targetMeshVerticesER[selectedEdges[edgeID].first])) {
                  vid = selectedEdges[edgeID].first;
                }
                else {
                  if (cuttingSegments[edgeID].cuttingMeshTargetVertexSegments[1].has_on(targetMeshVerticesER[selectedEdges[edgeID].second]) == false) {
                    std::cout << edgeID << ',' << ti << std::endl;
                    throw std::runtime_error("impossible");
                  }
                  vid = selectedEdges[edgeID].second;
                }

                corefinedPatches3D[ti].vertexIsOnTargetMesh[vi][0] = 0;
                corefinedPatches3D[ti].vertexIsOnTargetMesh[vi][1] = vid;
                onSeg = true;
              }
            }

            if (!onSeg && (cuttingTri0.has_on(corefinedPatches3D[ti].vertices[vi]) || cuttingTri1.has_on(corefinedPatches3D[ti].vertices[vi]))) {
              corefinedPatches3D[ti].vertexIsOnTargetMesh[vi][0] = 1;
              corefinedPatches3D[ti].vertexIsOnTargetMesh[vi][1] = edgeID;
            }
          }
        }
      }
    }
    else if (rType == PrismVertexRegionType::CYLINDER0 ||
      rType == PrismVertexRegionType::CYLINDER1 ||
      rType == PrismVertexRegionType::CYLINDER2) {
      int skeletonEID = (int)rType - (int)PrismVertexRegionType::CYLINDER0;

      const K::Point_2 &uv_v0 = prismMeshProjectedTriangles[ti].triangle2D[0];
      const K::Point_2 &uv_v1 = prismMeshProjectedTriangles[ti].triangle2D[1];
      const K::Point_2 &uv_v2 = prismMeshProjectedTriangles[ti].triangle2D[2];

      corefinedPatches2D[ti].rawTriangles.emplace_back(uv_v0, uv_v1, uv_v2);

      if (triangleIntersectedEdges[ti].size() == 0ull) {
        corefinedPatches2D[ti].vertices = { uv_v0, uv_v1, uv_v2 };
        corefinedPatches2D[ti].triangles.emplace_back(std::array<int, 3>{ 0, 1, 2 });

        corefinedPatches2D[ti].edgeVertices[0] = { 0, 1 };
        corefinedPatches2D[ti].edgeVertices[1] = { 1, 2 };
        corefinedPatches2D[ti].edgeVertices[2] = { 2, 0 };

        corefinedPatches2D[ti].cornerVertices[0] = 0;
        corefinedPatches2D[ti].cornerVertices[1] = 1;
        corefinedPatches2D[ti].cornerVertices[2] = 2;

        corefinedPatches2D[ti].vertexIsOnTargetMesh.assign(3, std::array<int, 2>{ -1, -1 });
      }
      else {
        // bfs
        for (int ei = 0; ei < (int)triangleIntersectedEdges[ti].size(); ei++) {
          int edgeID = triangleIntersectedEdges[ti][ei];
          PGO_ALOG(cuttingSegments[edgeID].cylinderSegPieceDoesExist[skeletonEID] == 1);
          const K::Segment_2 &seg = cuttingSegments[edgeID].cylinderProjections[skeletonEID];

          auto iter = corefinedPatches2D[ti].rawTriangles.begin();
          size_t count = corefinedPatches2D[ti].rawTriangles.size();

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
                    corefinedPatches2D[ti].rawTriangles.emplace_back(tri);
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
                    corefinedPatches2D[ti].rawTriangles.emplace_back(tri);
                  }
                }
              }
              else if (const K::Point_2 *p = boost::get<K::Point_2>(&*result)) {
                // create 3 triangles
                for (int i = 0; i < 3; i++) {
                  K::Triangle_2 tri(*p, (*iter)[i], (*iter)[(i + 1) % 3]);

                  if (!tri.is_degenerate()) {
                    corefinedPatches2D[ti].rawTriangles.emplace_back(tri);
                  }
                }
              }
              else {
                throw std::runtime_error("Impossible");
              }

              corefinedPatches2D[ti].rawTriangles.pop_front();
              iter = corefinedPatches2D[ti].rawTriangles.begin();
            }
            else {
              corefinedPatches2D[ti].rawTriangles.splice(corefinedPatches2D[ti].rawTriangles.end(), corefinedPatches2D[ti].rawTriangles, corefinedPatches2D[ti].rawTriangles.begin());
              iter = corefinedPatches2D[ti].rawTriangles.begin();
            }
          }  // end result
        }

        corefinedPatches2D[ti].vertices.reserve(corefinedPatches2D[ti].rawTriangles.size() * 3);
        corefinedPatches2D[ti].triangles.reserve(corefinedPatches2D[ti].rawTriangles.size());
        for (const auto &tri : corefinedPatches2D[ti].rawTriangles) {
          std::array<int, 3> triVtx = { -1, -1, -1 };
          for (int j = 0; j < 3; j++) {
            // find if the vertex exist
            for (int vi = 0; vi < (int)corefinedPatches2D[ti].vertices.size(); vi++) {
              const auto &vtx = corefinedPatches2D[ti].vertices[vi];
              if (tri[j] == vtx) {
                triVtx[j] = vi;
                break;
              }
            }

            // if not found
            if (triVtx[j] < 0) {
              corefinedPatches2D[ti].vertices.emplace_back(tri[j]);
              triVtx[j] = (int)corefinedPatches2D[ti].vertices.size() - 1;
            }
          }

          corefinedPatches2D[ti].triangles.emplace_back(triVtx);
        }

        corefinedPatches2D[ti].edgeVertices[0].reserve(corefinedPatches2D[ti].vertices.size());
        corefinedPatches2D[ti].edgeVertices[1].reserve(corefinedPatches2D[ti].vertices.size());
        corefinedPatches2D[ti].edgeVertices[2].reserve(corefinedPatches2D[ti].vertices.size());

        K::Segment_2 e01(uv_v0, uv_v1), e12(uv_v1, uv_v2), e20(uv_v2, uv_v0);
        for (int vi = 0; vi < (int)corefinedPatches2D[ti].vertices.size(); vi++) {
          if (e01.has_on(corefinedPatches2D[ti].vertices[vi])) {
            corefinedPatches2D[ti].edgeVertices[0].emplace_back(vi);
          }

          if (e12.has_on(corefinedPatches2D[ti].vertices[vi])) {
            corefinedPatches2D[ti].edgeVertices[1].emplace_back(vi);
          }

          if (e20.has_on(corefinedPatches2D[ti].vertices[vi])) {
            corefinedPatches2D[ti].edgeVertices[2].emplace_back(vi);
          }

          if (corefinedPatches2D[ti].vertices[vi] == uv_v0) {
            corefinedPatches2D[ti].cornerVertices[0] = vi;
          }
          else if (corefinedPatches2D[ti].vertices[vi] == uv_v1) {
            corefinedPatches2D[ti].cornerVertices[1] = vi;
          }
          else if (corefinedPatches2D[ti].vertices[vi] == uv_v2) {
            corefinedPatches2D[ti].cornerVertices[2] = vi;
          }
        }

        // record the target edge ID for each newly generated triangle edges
        corefinedPatches2D[ti].triangleEdgeIDs.reserve((int)corefinedPatches2D[ti].triangles.size());

        for (int tii = 0; tii < (int)corefinedPatches2D[ti].triangles.size(); tii++) {
          for (int j = 0; j < 3; j++) {
            const K::Point_2 &p0 = corefinedPatches2D[ti].vertices[corefinedPatches2D[ti].triangles[tii][j]];
            const K::Point_2 &p1 = corefinedPatches2D[ti].vertices[corefinedPatches2D[ti].triangles[tii][(j + 1) % 3]];

            for (int ei = 0; ei < (int)triangleIntersectedEdges[ti].size(); ei++) {
              int edgeID = triangleIntersectedEdges[ti][ei];

              if (cuttingSegments[edgeID].cylinderProjections[skeletonEID].has_on(p0) &&
                cuttingSegments[edgeID].cylinderProjections[skeletonEID].has_on(p1)) {
                corefinedPatches2D[ti].triangleEdgeIDs.emplace_back(tii, j, triangleIntersectedEdges[ti][ei]);
                break;
              }
            }
          }
        }

        corefinedPatches2D[ti].vertexIsOnTargetMesh.assign(corefinedPatches2D[ti].vertices.size(), std::array<int, 2>{ -1, -1 });
        for (int vi = 0; vi < (int)corefinedPatches2D[ti].vertices.size(); vi++) {
          for (int ei = 0; ei < (int)triangleIntersectedEdges[ti].size(); ei++) {
            int edgeID = triangleIntersectedEdges[ti][ei];
            const K::Segment_2 &seg = cuttingSegments[edgeID].cylinderProjections[skeletonEID];
            bool onSeg = false;
            if (cuttingSegments[edgeID].segPieceEndIsTargetVertex[skeletonEID][0]) {
              if (seg[0] == corefinedPatches2D[ti].vertices[vi]) {
                corefinedPatches2D[ti].vertexIsOnTargetMesh[vi][0] = 0;
                corefinedPatches2D[ti].vertexIsOnTargetMesh[vi][1] = selectedEdges[edgeID].first;
                onSeg = true;
              }
            }

            if (cuttingSegments[edgeID].segPieceEndIsTargetVertex[skeletonEID][1]) {
              if (seg[1] == corefinedPatches2D[ti].vertices[vi]) {
                corefinedPatches2D[ti].vertexIsOnTargetMesh[vi][0] = 0;
                corefinedPatches2D[ti].vertexIsOnTargetMesh[vi][1] = selectedEdges[edgeID].second;
                onSeg = true;
              }
            }

            if (!onSeg && seg.has_on(corefinedPatches2D[ti].vertices[vi])) {
              corefinedPatches2D[ti].vertexIsOnTargetMesh[vi][0] = 1;
              corefinedPatches2D[ti].vertexIsOnTargetMesh[vi][1] = edgeID;
            }
          }
        }
      }
    }
    else {
      throw std::runtime_error("impossible");
    }

    int cc = counter.fetch_add(1);
    if (cc % 100 == 0)
      std::cout << cc << '/' << prismMesh.numTriangles() << "  " << std::flush;
  });

  std::cout << std::endl;

  t2 = HClock::now();
  SPDLOG_LOGGER_INFO(pgo::Logging::lgr(), "S1 Done. Time cost: {:.4f}s", duraSecond(t1, t2));

  if constexpr (dumpMesh) {
    pgo::Mesh::TriMeshGeo rr;
    for (int ti = 0; ti < prismMesh.numTriangles(); ti++) {
      PrismVertexRegionType rType = prismMeshProjectedTriangles[ti].regionType;
      if (rType == PrismVertexRegionType::TRIANGLE_POS ||
        rType == PrismVertexRegionType::TRIANGLE_NEG) {
        pgo::Mesh::TriMeshGeo zz;
        for (int vi = 0; vi < corefinedPatches3D[ti].vertices.size(); vi++) {
          zz.addPos(toVec3(corefinedPatches3D[ti].vertices[vi]));
        }

        for (int tt = 0; tt < corefinedPatches3D[ti].triangles.size(); tt++) {
          zz.addTri(ES::V3i(corefinedPatches3D[ti].triangles[tt][0], corefinedPatches3D[ti].triangles[tt][1], corefinedPatches3D[ti].triangles[tt][2]));
        }
        rr.addMesh(zz);
      }
      else if (rType == PrismVertexRegionType::CYLINDER0 ||
        rType == PrismVertexRegionType::CYLINDER1 ||
        rType == PrismVertexRegionType::CYLINDER2) {
        pgo::Mesh::TriMeshGeo zz;
        for (int vi = 0; vi < corefinedPatches2D[ti].vertices.size(); vi++) {
          K::Point_3 pp;
          if (pt2D_is_on_border(corefinedPatches2D[ti].vertices[vi], ti, pp)) {
            zz.addPos(toVec3(pp));
          }
          else {
            IK::Point_3 p = pt2D_to_pt3D(corefinedPatches2D[ti].vertices[vi], ti);
            zz.addPos(toVec3(p));
          }
        }

        for (int tt = 0; tt < corefinedPatches2D[ti].triangles.size(); tt++) {
          zz.addTri(ES::V3i(corefinedPatches2D[ti].triangles[tt][0], corefinedPatches2D[ti].triangles[tt][1], corefinedPatches2D[ti].triangles[tt][2]));
        }
        rr.addMesh(zz);
      }
      else {
      }
    }
    rr.save("c.obj");
  }

  // ================================================================
  // next we merge all cut triangles into a closed manifold mesh
  pgo::Mesh::TriMeshNeighbor prismMeshNeighbor(prismMesh);

  std::vector<int> triIDQ;
  triIDQ.reserve(prismMesh.numTriangles() * 2);

  std::vector<int> isTriVisited(prismMesh.numTriangles(), 0);

  std::vector<K::Point_2> finalPoints2D;
  finalPoints2D.reserve(prismMesh.numTriangles() * 10);

  std::vector<K::Point_3> finalPoints3D;
  finalPoints3D.reserve(prismMesh.numTriangles() * 10);

  std::vector<ES::V3i> finalTriangles;
  finalTriangles.reserve(prismMesh.numTriangles() * 10);

  std::map<std::pair<int, int>, std::vector<int>> edgeVertexIDs;
  std::set<std::pair<int, int>> edgeIsBorder;
  std::vector<int> newVertexIDs(prismMesh.numVertices(), -1);

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
    std::cout << triIDCur << std::endl;
    PrismVertexRegionType rType = prismMeshProjectedTriangles[triIDCur].regionType;
    bool isCylinder = false;
    if (rType == PrismVertexRegionType::TRIANGLE_POS ||
      rType == PrismVertexRegionType::TRIANGLE_NEG) {
      // this one store the vertex global IDs
      // if it is -1, it means it is not addressed or found
      vertexMasks.assign(corefinedPatches3D[triIDCur].vertices.size(), -1);
      isCylinder = false;
    }
    else if (rType == PrismVertexRegionType::CYLINDER0 ||
      rType == PrismVertexRegionType::CYLINDER1 ||
      rType == PrismVertexRegionType::CYLINDER2) {
      // this one store the vertex global IDs
      // if it is -1, it means it is not addressed or found
      vertexMasks.assign(corefinedPatches2D[triIDCur].vertices.size(), -1);
      isCylinder = true;
    }

    int edgeMask[3] = { 0, 0, 0 };
    for (int j = 0; j < 3; j++) {
      std::pair<int, int> edgeID{ prismMesh.triVtxID(triIDCur, j), prismMesh.triVtxID(triIDCur, (j + 1) % 3) };
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
        if (isCylinder) {
          if (eitt->second.size() != corefinedPatches2D[triIDCur].edgeVertices[j].size()) {
            std::vector<IK::Point_3> vtx;
            for (int i = 0; i < (int)corefinedPatches2D[triIDCur].vertices.size(); i++) {
              vtx.emplace_back(pt2D_to_pt3D(corefinedPatches2D[triIDCur].vertices[i], triIDCur));
            }

            pgo::CGALInterface::Polyhedron<IK> finalMesh;
            pgo::CGALInterface::PolyhedronBuilderNonTriangle<IK, std::vector<IK::Point_3>::iterator, std::vector<std::array<int, 3>>::iterator> builder(
              vtx.begin(), corefinedPatches2D[triIDCur].triangles.begin(), (int)vtx.size(), (int)corefinedPatches2D[triIDCur].triangles.size());
            finalMesh.delegate(builder);
            // std::ofstream("rzr1.off") << finalMesh;

            if (1) {
              pgo::CGALInterface::Polyhedron<K> finalMesh;
              pgo::CGALInterface::PolyhedronBuilderNonTriangle<K, std::vector<K::Point_3>::iterator, std::vector<ES::V3i>::iterator> builder(
                finalPoints3D.begin(), finalTriangles.begin(), (int)finalPoints3D.size(), (int)finalTriangles.size());
              finalMesh.delegate(builder);
              // std::ofstream("rzr0.off") << finalMesh;
            }

            abort();
          }

          // find global vertex ID
          for (int k = 0; k < (int)corefinedPatches2D[triIDCur].edgeVertices[j].size(); k++) {
            if (vertexMasks[corefinedPatches2D[triIDCur].edgeVertices[j][k]] >= 0)
              continue;

            // check if it is region border vertices
            bool found = false;
            K::Point_3 pp3;
            if (pt2D_is_on_border(corefinedPatches2D[triIDCur].vertices[corefinedPatches2D[triIDCur].edgeVertices[j][k]], triIDCur, pp3)) {
              // std::cout << "C1\n";
              // std::cout << "pp3: " << toVec3(pp3) << std::endl;

              for (int r = 0; r < (int)eitt->second.size(); r++) {
                if (pp3 == finalPoints3D[eitt->second[r]]) {
                  vertexMasks[corefinedPatches2D[triIDCur].edgeVertices[j][k]] = eitt->second[r];
                  found = true;
                  break;
                }
              }
            }
            // if it is not
            else {
              // std::cout << "C2\n";
              const K::Point_2 &pt2D = corefinedPatches2D[triIDCur].vertices[corefinedPatches2D[triIDCur].edgeVertices[j][k]];
              // std::cout << tod(pt2D[0]) << ',' << tod(pt2D[1]) << std::endl;

              int prismMeshVertexID = -1;
              for (int c = 0; c < 3; c++) {
                if (corefinedPatches2D[triIDCur].edgeVertices[j][k] == corefinedPatches2D[triIDCur].cornerVertices[c]) {
                  int originVID = prismMesh.triVtxID(triIDCur, c);
                  if (isRegionBorderVertex(originVID)) {
                    prismMeshVertexID = originVID;
                    break;
                  }
                }
              }
              // std::cout << "is corner: " << prismMeshVertexID << std::endl;

              // if a triangle is not a border triangle,
              // it is possible that there is one vertex on the region border,
              // this 'if' handles this case.
              if (prismMeshVertexID >= 0) {
                // std::cout << "corner pos: " << toVec3(prismMeshVerticesEK[prismMeshVertexID]) << std::endl;

                for (int r = 0; r < (int)eitt->second.size(); r++) {
                  if (prismMeshVerticesEK[prismMeshVertexID] == finalPoints3D[eitt->second[r]]) {
                    vertexMasks[corefinedPatches2D[triIDCur].edgeVertices[j][k]] = eitt->second[r];
                    found = true;
                    break;
                  }
                }
              }
              else {
                for (int r = 0; r < (int)eitt->second.size(); r++) {
                  if (corefinedPatches2D[triIDCur].vertices[corefinedPatches2D[triIDCur].edgeVertices[j][k]] == finalPoints2D[eitt->second[r]]) {
                    vertexMasks[corefinedPatches2D[triIDCur].edgeVertices[j][k]] = eitt->second[r];
                    found = true;
                    break;
                  }
                }
              }
            }
            std::cout << std::flush;

            if (found != true) {
              std::cout << "3d pt:\n";
              for (int r = 0; r < (int)eitt->second.size(); r++) {
                std::cout << r << ':' << toVec3(finalPoints3D[eitt->second[r]]) << std::endl;
              }
              std::cout << std::endl;

              std::cout << "2d pt:\n";
              for (int r = 0; r < (int)eitt->second.size(); r++) {
                std::cout << r << ':' << tod(finalPoints2D[eitt->second[r]][0]) << ',' << tod(finalPoints2D[eitt->second[r]][1]) << std::endl;
              }
              std::cout << std::endl;

              std::vector<IK::Point_3> vtx;
              for (int i = 0; i < (int)corefinedPatches2D[triIDCur].vertices.size(); i++) {
                vtx.emplace_back(pt2D_to_pt3D(corefinedPatches2D[triIDCur].vertices[i], triIDCur));
              }

              pgo::CGALInterface::Polyhedron<IK> finalMesh;
              pgo::CGALInterface::PolyhedronBuilderNonTriangle<IK, std::vector<IK::Point_3>::iterator, std::vector<std::array<int, 3>>::iterator> builder(
                vtx.begin(), corefinedPatches2D[triIDCur].triangles.begin(), (int)vtx.size(), (int)corefinedPatches2D[triIDCur].triangles.size());
              finalMesh.delegate(builder);

              std::ofstream afile("aa.obj");
              for (auto vit = finalMesh.vertices_begin(); vit != finalMesh.vertices_end(); ++vit) {
                afile << "v " << tod(vit->point()[0]) << ' ' << tod(vit->point()[1]) << ' ' << tod(vit->point()[2]) << '\n';
              }

              for (auto fit = finalMesh.facets_begin(); fit != finalMesh.facets_end(); ++fit) {
                afile << "f " << (fit->halfedge()->vertex()->id() + 1);
                afile << " " << (fit->halfedge()->next()->vertex()->id() + 1);
                afile << " " << (fit->halfedge()->next()->next()->vertex()->id() + 1) << std::endl;
              }
              afile.close();

              if (1) {
                pgo::CGALInterface::Polyhedron<K> finalMesh;
                pgo::CGALInterface::PolyhedronBuilderNonTriangle<K, std::vector<K::Point_3>::iterator, std::vector<ES::V3i>::iterator> builder(
                  finalPoints3D.begin(), finalTriangles.begin(), (int)finalPoints3D.size(), (int)finalTriangles.size());
                finalMesh.delegate(builder);
                std::ofstream("aa1.off") << finalMesh;
              }

              // std::ofstream("rzr1.off") << finalMesh;
              abort();
            }
          }
        }
        else {
          if (eitt->second.size() != corefinedPatches3D[triIDCur].edgeVertices[j].size()) {
            pgo::CGALInterface::Polyhedron<K> finalMesh;
            pgo::CGALInterface::PolyhedronBuilderNonTriangle<K, std::vector<K::Point_3>::iterator, std::vector<std::array<int, 3>>::iterator> builder(
              corefinedPatches3D[triIDCur].vertices.begin(), corefinedPatches3D[triIDCur].triangles.begin(), (int)corefinedPatches3D[triIDCur].vertices.size(), (int)corefinedPatches3D[triIDCur].triangles.size());
            finalMesh.delegate(builder);
            std::ofstream("aa1.off") << finalMesh;

            if (1) {
              pgo::CGALInterface::Polyhedron<K> finalMesh;
              pgo::CGALInterface::PolyhedronBuilderNonTriangle<K, std::vector<K::Point_3>::iterator, std::vector<ES::V3i>::iterator> builder(
                finalPoints3D.begin(), finalTriangles.begin(), (int)finalPoints3D.size(), (int)finalTriangles.size());
              finalMesh.delegate(builder);
              std::ofstream("aa0.off") << finalMesh;
            }

            abort();
          }

          // find global vertex ID
          for (int k = 0; k < (int)corefinedPatches3D[triIDCur].edgeVertices[j].size(); k++) {
            if (vertexMasks[corefinedPatches3D[triIDCur].edgeVertices[j][k]] >= 0)
              continue;

            bool found = false;
            for (int r = 0; r < (int)eitt->second.size(); r++) {
              if (corefinedPatches3D[triIDCur].vertices[corefinedPatches3D[triIDCur].edgeVertices[j][k]] == finalPoints3D[eitt->second[r]]) {
                vertexMasks[corefinedPatches3D[triIDCur].edgeVertices[j][k]] = eitt->second[r];
                found = true;
                break;
              }
            }

            if (found != true) {
              std::cout << "p3d: " << std::endl;
              for (int r = 0; r < (int)eitt->second.size(); r++) {
                std::cout << toVec3(finalPoints3D[eitt->second[r]]) << std::endl;
              }
              std::cout << "p to found: " << toVec3(corefinedPatches3D[triIDCur].vertices[corefinedPatches3D[triIDCur].edgeVertices[j][k]]) << std::endl;

              pgo::CGALInterface::Polyhedron<K> finalMesh;
              pgo::CGALInterface::PolyhedronBuilderNonTriangle<K, std::vector<K::Point_3>::iterator, std::vector<std::array<int, 3>>::iterator> builder(
                corefinedPatches3D[triIDCur].vertices.begin(), corefinedPatches3D[triIDCur].triangles.begin(),
                (int)corefinedPatches3D[triIDCur].vertices.size(), (int)corefinedPatches3D[triIDCur].triangles.size());
              finalMesh.delegate(builder);
              std::ofstream("aa1.off") << finalMesh;

              if (1) {
                pgo::CGALInterface::Polyhedron<K> finalMesh;
                pgo::CGALInterface::PolyhedronBuilderNonTriangle<K, std::vector<K::Point_3>::iterator, std::vector<ES::V3i>::iterator> builder(
                  finalPoints3D.begin(), finalTriangles.begin(), (int)finalPoints3D.size(), (int)finalTriangles.size());
                finalMesh.delegate(builder);
                std::ofstream("aa0.off") << finalMesh;
              }

              abort();
            }
          }
        }
      }
    }  // end for j \in [0, 1, 2]

    if (isCylinder) {
      // address original triangle corners
      for (int j = 0; j < 3; j++) {
        int vi = corefinedPatches2D[triIDCur].cornerVertices[j];
        if (newVertexIDs[prismMesh.triVtxID(triIDCur, j)] >= 0) {
          if (vertexMasks[vi] >= 0) {
            PGO_ALOG(vertexMasks[vi] == newVertexIDs[prismMesh.triVtxID(triIDCur, j)]);
          }
          else {
            vertexMasks[vi] = newVertexIDs[prismMesh.triVtxID(triIDCur, j)];
          }
        }
      }

      // compute all vertex ID,
      // especially for the ones that has not been visited
      for (int j = 0; j < (int)corefinedPatches2D[triIDCur].vertices.size(); j++) {
        if (vertexMasks[j] < 0) {
          finalPoints2D.emplace_back(corefinedPatches2D[triIDCur].vertices[j]);

          K::Point_3 pp3;
          if (pt2D_is_on_border(corefinedPatches2D[triIDCur].vertices[j], triIDCur, pp3)) {
            finalPoints3D.emplace_back(pp3);
          }
          else {
            IK::Point_3 p = pt2D_to_pt3D(corefinedPatches2D[triIDCur].vertices[j], triIDCur);
            finalPoints3D.emplace_back(toP3_EK(p));
          }

          // corefinedPatches[triIDCur].vertices[j] to 3D
          vertexMasks[j] = (int)finalPoints2D.size() - 1;
        }
      }

      // address original triangle corners again
      for (int j = 0; j < 3; j++) {
        int vi = corefinedPatches2D[triIDCur].cornerVertices[j];
        if (newVertexIDs[prismMesh.triVtxID(triIDCur, j)] < 0) {
          newVertexIDs[prismMesh.triVtxID(triIDCur, j)] = vertexMasks[vi];
        }
      }

      // adding delete boundary info and edge info
      for (int j = 0; j < 3; j++) {
        std::pair<int, int> edgeID{ prismMesh.triVtxID(triIDCur, j), prismMesh.triVtxID(triIDCur, (j + 1) % 3) };
        if (edgeID.first > edgeID.second)
          std::swap(edgeID.first, edgeID.second);

        // if it is an edge that is visited for the first time
        if (edgeMask[j]) {
          auto ret = edgeVertexIDs.emplace(edgeID, corefinedPatches2D[triIDCur].edgeVertices[j]);
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
      triangleNewIDs.assign(corefinedPatches2D[triIDCur].triangles.size(), -1);
      for (int j = 0; j < (int)corefinedPatches2D[triIDCur].triangles.size(); j++) {
        ES::V3i triID = ES::Mp<const ES::V3i>(corefinedPatches2D[triIDCur].triangles[j].data());
        for (int k = 0; k < 3; k++) {
          triID[k] = vertexMasks[triID[k]];
        }
        finalTriangles.emplace_back(triID);
        triangleNewIDs[j] = (int)finalTriangles.size() - 1;
      }

      // store edge info
      for (int j = 0; j < (int)corefinedPatches2D[triIDCur].triangleEdgeIDs.size(); j++) {
        int newTriID = std::get<0>(corefinedPatches2D[triIDCur].triangleEdgeIDs[j]);
        newTriID = triangleNewIDs[newTriID];

        int ithEdge = std::get<1>(corefinedPatches2D[triIDCur].triangleEdgeIDs[j]);
        int edgeID = std::get<2>(corefinedPatches2D[triIDCur].triangleEdgeIDs[j]);
        targetEdgeOverlappingTriangleIDs[edgeID].emplace_back(std::array<int, 2>{ newTriID, ithEdge });
      }

      // store vtx info
      for (int j = 0; j < (int)corefinedPatches2D[triIDCur].vertices.size(); j++) {
        if (corefinedPatches2D[triIDCur].vertexIsOnTargetMesh[j][0] >= 0) {
          auto it = vertexIsOnTargetMesh.lower_bound(vertexMasks[j]);
          if (it != vertexIsOnTargetMesh.end() && it->first == vertexMasks[j]) {
          }
          else {
            vertexIsOnTargetMesh.emplace_hint(it, vertexMasks[j], corefinedPatches2D[triIDCur].vertexIsOnTargetMesh[j]);
          }
        }
      }
    }
    else {
      // address original triangle corners
      for (int j = 0; j < 3; j++) {
        int vi = corefinedPatches3D[triIDCur].cornerVertices[j];
        if (newVertexIDs[prismMesh.triVtxID(triIDCur, j)] >= 0) {
          if (vertexMasks[vi] >= 0) {
            PGO_ALOG(vertexMasks[vi] == newVertexIDs[prismMesh.triVtxID(triIDCur, j)]);
          }
          else {
            vertexMasks[vi] = newVertexIDs[prismMesh.triVtxID(triIDCur, j)];
          }
        }
      }

      // compute all vertex ID,
      // especially for the ones that has not been visited
      for (int j = 0; j < (int)corefinedPatches3D[triIDCur].vertices.size(); j++) {
        if (vertexMasks[j] < 0) {
          finalPoints2D.emplace_back(-1.0, -1.0);
          finalPoints3D.emplace_back(corefinedPatches3D[triIDCur].vertices[j]);

          // corefinedPatches[triIDCur].vertices[j] to 3D
          vertexMasks[j] = (int)finalPoints2D.size() - 1;
        }
      }

      // address original triangle corners again
      for (int j = 0; j < 3; j++) {
        int vi = corefinedPatches3D[triIDCur].cornerVertices[j];
        if (newVertexIDs[prismMesh.triVtxID(triIDCur, j)] < 0) {
          newVertexIDs[prismMesh.triVtxID(triIDCur, j)] = vertexMasks[vi];
        }
      }

      // adding delete boundary info and edge info
      for (int j = 0; j < 3; j++) {
        std::pair<int, int> edgeID{ prismMesh.triVtxID(triIDCur, j), prismMesh.triVtxID(triIDCur, (j + 1) % 3) };
        if (edgeID.first > edgeID.second)
          std::swap(edgeID.first, edgeID.second);

        // if it is an edge that is visited for the first time
        if (edgeMask[j]) {
          auto ret = edgeVertexIDs.emplace(edgeID, corefinedPatches3D[triIDCur].edgeVertices[j]);
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
      triangleNewIDs.assign(corefinedPatches3D[triIDCur].triangles.size(), -1);
      for (int j = 0; j < (int)corefinedPatches3D[triIDCur].triangles.size(); j++) {
        ES::V3i triID = ES::Mp<const ES::V3i>(corefinedPatches3D[triIDCur].triangles[j].data());
        for (int k = 0; k < 3; k++) {
          triID[k] = vertexMasks[triID[k]];
        }
        finalTriangles.emplace_back(triID);
        triangleNewIDs[j] = (int)finalTriangles.size() - 1;
      }

      // store edge info
      for (int j = 0; j < (int)corefinedPatches3D[triIDCur].triangleEdgeIDs.size(); j++) {
        int newTriID = std::get<0>(corefinedPatches3D[triIDCur].triangleEdgeIDs[j]);
        newTriID = triangleNewIDs[newTriID];

        int ithEdge = std::get<1>(corefinedPatches3D[triIDCur].triangleEdgeIDs[j]);
        int edgeID = std::get<2>(corefinedPatches3D[triIDCur].triangleEdgeIDs[j]);
        targetEdgeOverlappingTriangleIDs[edgeID].emplace_back(std::array<int, 2>{ newTriID, ithEdge });
      }

      // store vtx info
      for (int j = 0; j < (int)corefinedPatches3D[triIDCur].vertices.size(); j++) {
        if (corefinedPatches3D[triIDCur].vertexIsOnTargetMesh[j][0] >= 0) {
          auto it = vertexIsOnTargetMesh.lower_bound(vertexMasks[j]);
          if (it != vertexIsOnTargetMesh.end() && it->first == vertexMasks[j]) {
          }
          else {
            vertexIsOnTargetMesh.emplace_hint(it, vertexMasks[j], corefinedPatches3D[triIDCur].vertexIsOnTargetMesh[j]);
          }
        }
      }
    }

    // expand bfs search
    ES::V3i triN = prismMeshNeighbor.getTriangleNeighbors(triIDCur);
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
    pgo::CGALInterface::Polyhedron<K> finalMesh;
    pgo::CGALInterface::PolyhedronBuilderNonTriangle<K, std::vector<K::Point_3>::iterator, std::vector<ES::V3i>::iterator> builder(
      finalPoints3D.begin(), finalTriangles.begin(), (int)finalPoints3D.size(), (int)finalTriangles.size());
    finalMesh.delegate(builder);
    std::ofstream("rzr0.off") << finalMesh;
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
            K::Point_3 p = finalPoints3D[vi];
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

    const K::Point_3 &p = finalPoints3D[v.first];
    finalPoints1.emplace_back(p);
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
      auto projPt = computePrismProjectionPt(center, skeletonTri);
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
      K::FT center[3] = { 0, 0, 0 };
      int counter = 0;
      do {
        const K::Point_3 &pt = h1->vertex()->point();
        for (int j = 0; j < 3; j++) {
          center[j] += pt[j];
        }
        counter += 1;
        h1 = h1->next();
      } while (h0 != h1);

      int vidCur = (int)finalMesh.size_of_vertices();
      Poly::Halfedge_handle h = finalMesh.create_center_vertex(fi->halfedge());
      h->vertex()->point() = K::Point_3(center[0] / counter, center[1] / counter, center[2] / counter);
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

  prismVertexProjection(finalMesh, targetMesh, targetMeshPoly, targetMeshBVTreeExact,
    skeletonTri, maxAllowedRadiusDist, vtxIDNew2Old, vertexIsOnTargetMesh, selectedEdges);

  if constexpr (dumpMesh) {
    std::ofstream("rzr2.off") << finalMesh;
  }

  std::unordered_set<int> borderVertexIDs;
  int numVertices = (int)finalMesh.size_of_vertices();

  // create a closed mesh
  if constexpr (1) {
    // find loop 0
    Poly::Halfedge_handle border_h[3];
    int cids[3];
    std::set<Poly::Halfedge_handle> border_h_set[3];

    for (int j = 0; j < 3; j++) {
      for (auto hit = finalMesh.halfedges_begin(); hit != finalMesh.halfedges_end(); ++hit) {
        // check if the halfedge is a previously visited edge
        bool visited = false;
        for (int k = 0; k < j; k++) {
          if (border_h_set[k].find(hit) != border_h_set[k].end()) {
            visited = true;
            break;
          }
        }

        if (visited)
          continue;

        if (hit->is_border()) {
          border_h[j] = hit;
          break;
        }
      }

      Poly::Halfedge_handle hh = border_h[j];
      do {
        border_h_set[j].emplace(hh);
        hh = hh->next();
      } while (hh != border_h[j]);
    }

    for (int i = 0; i < 3; i++) {
      for (const auto &pr : border_h_set[i]) {
        borderVertexIDs.emplace(pr->vertex()->id());
      }
    }

    for (int i = 0; i < 3; i++) {
      auto proj = computePrismProjectionPt(border_h[i]->vertex()->point(), skeletonTri);
      K::FT d2[3] = {
        (std::get<0>(proj) - primitiveCenterEK[0]).squared_length(),
        (std::get<0>(proj) - primitiveCenterEK[1]).squared_length(),
        (std::get<0>(proj) - primitiveCenterEK[2]).squared_length()
      };

      cids[i] = std::ptrdiff_t(std::min_element(d2, d2 + 3) - d2);
    }
    PGO_ALOG(cids[0] != cids[1] && cids[1] != cids[2] && cids[2] != cids[0]);

    int vidCur = (int)finalMesh.size_of_vertices();
    std::cout << vidCur << std::endl;

    for (int j = 0; j < 3; j++) {
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
    else {
    }

    if (borderVertexIDs.find(vid) != borderVertexIDs.end()) {
      finalMeshVertexProperties[vid].isBorder = 1;
    }
  }

  cleanUpMesh(targetMeshPoly, targetMeshBVTreeExact, finalMesh, 5e-3, finalMeshVertexProperties, targetEdges);

  if constexpr (dumpMesh) {
    std::ofstream("rzr4.off") << finalMesh;
  }

  double step[2] = { 0.5, -0.6 };
  for (int iter = 0; iter < 0; iter++) {
    for (auto vit = finalMesh.vertices_begin(); vit != finalMesh.vertices_end(); ++vit) {
      int vid = (int)vit->id();

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
      K::FT center[3] = { 0, 0, 0 };
      int counter = 0;
      do {
        const K::Point_3 &pt = h1->opposite()->vertex()->point();
        center[0] += pt[0];
        center[1] += pt[1];
        center[2] += pt[2];
        counter += 1;

        h1 = h1->next()->opposite();
      } while (h != h1);

      K::Point_3 tgt(center[0] / counter, center[1] / counter, center[2] / counter);
      K::Vector_3 diff = tgt - vit->point();
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

  return 0;
}

void MedialAxisRepresentation::fillPrism(const pgo::Mesh::TriMeshGeo &prismMesh, const ES::V3d centers[3], pgo::Mesh::TriMeshGeo &meshOut)
{
  K::Point_3 primitiveCenterEK[3] = { toP3_EK(centers[0]), toP3_EK(centers[1]), toP3_EK(centers[2]) };
  K::Triangle_3 skeletonTri(primitiveCenterEK[0], primitiveCenterEK[1], primitiveCenterEK[2]);

  Poly finalMesh;
  pgo::CGALInterface::triangleMesh2Polyhedron(prismMesh, finalMesh);

  // create a closed mesh
  if constexpr (1) {
    // find loop 0
    Poly::Halfedge_handle border_h[3];
    int cids[3];
    std::set<Poly::Halfedge_handle> border_h_set[3];

    for (int j = 0; j < 3; j++) {
      for (auto hit = finalMesh.halfedges_begin(); hit != finalMesh.halfedges_end(); ++hit) {
        // check if the halfedge is a previously visited edge
        bool visited = false;
        for (int k = 0; k < j; k++) {
          if (border_h_set[k].find(hit) != border_h_set[k].end()) {
            visited = true;
            break;
          }
        }

        if (visited)
          continue;

        if (hit->is_border()) {
          border_h[j] = hit;
          break;
        }
      }

      Poly::Halfedge_handle hh = border_h[j];
      do {
        border_h_set[j].emplace(hh);
        hh = hh->next();
      } while (hh != border_h[j]);
    }

    for (int i = 0; i < 3; i++) {
      auto proj = computePrismProjectionPt(border_h[i]->vertex()->point(), skeletonTri);
      K::FT d2[3] = {
        (std::get<0>(proj) - primitiveCenterEK[0]).squared_length(),
        (std::get<0>(proj) - primitiveCenterEK[1]).squared_length(),
        (std::get<0>(proj) - primitiveCenterEK[2]).squared_length()
      };

      cids[i] = std::ptrdiff_t(std::min_element(d2, d2 + 3) - d2);
    }
    PGO_ALOG(cids[0] != cids[1] && cids[1] != cids[2] && cids[2] != cids[0]);

    int vidCur = (int)finalMesh.size_of_vertices();
    std::cout << vidCur << std::endl;

    for (int j = 0; j < 3; j++) {
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
}