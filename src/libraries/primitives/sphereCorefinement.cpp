#include "sphereCorefinement.h"
#include "corefinementUtilities.h"

#include "pgoLogging.h"
#include "basicAlgorithms.h"
#include "createTriMesh.h"
#include "triMeshNeighbor.h"
#include "EigenSupport.h"

#include <tbb/enumerable_thread_specific.h>
#include <tbb/parallel_for.h>
#include <tbb/concurrent_hash_map.h>

#include <chrono>
#include <unordered_map>
#include <unordered_set>

namespace MedialAxisRepresentation
{
void vertexProjection(Poly &finalMesh, const pgo::Mesh::TriMeshGeo &targetMesh, const FaceGraphTree &targetMeshBVTreeExact,
  const K::Point_3 &centerER, double maxAllowedRadiusDist,
  const std::unordered_map<int, int> &vtxIDNew2Old, const std::map<int, std::array<int, 2>> &vertexIsOnTargetMesh,
  const std::vector<std::pair<int, int>> &selectedEdges)
{
  pgo::Mesh::TriMeshNeighbor targetMeshNeighbor(targetMesh);

  for (auto it = finalMesh.vertices_begin(); it != finalMesh.vertices_end(); ++it) {
    K::Point_3 pEK = it->point();
    const K::Point_3 &projPtEK = centerER;
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

    if (curLength <= len + maxAllowedRadiusDist) {
      it->point() = closestPtEK;
    }
    else {
      it->point() = projPtEK + dir * (len + maxAllowedRadiusDist) / len;
    }

#if 1
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
          else if (lenNew > len + maxAllowedRadiusDist) {
            it->point() = projPtEK + dir * (len + maxAllowedRadiusDist) / lenNew;
          }
        }
        else if (vit->second[0] == 1) {
          int vids[2] = { selectedEdges[vit->second[1]].first, selectedEdges[vit->second[1]].second };
          const K::Point_3 cpt = closestPtEK;
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
          else if (lenNew > len + maxAllowedRadiusDist) {
            it->point() = projPtEK + dir * (len + maxAllowedRadiusDist) / lenNew;
          }
        }
      }
    }
#endif
  }
}

}  // namespace MedialAxisRepresentation

int MedialAxisRepresentation::corefineSphereMeshWithTarget(const pgo::Mesh::TriMeshGeo &sphereMesh, const double center[3],
  const pgo::Mesh::TriMeshGeo &targetMesh, const pgo::Mesh::TriMeshBVTree &targetMeshBVTree, const pgo::Mesh::TriMeshPseudoNormal &targetMeshNormals,
  double maxConfidentDist, double maxAllowedRadiusDist, pgo::Mesh::TriMeshGeo &meshOut, const char *filename)
{
  constexpr int dumpMesh = 0;

  std::vector<int> targetMeshChoosenVertex(targetMesh.numVertices(), 0);
  tbb::enumerable_thread_specific<std::vector<std::tuple<double, double, int>>> closestDistDueryStackTLS;
  tbb::enumerable_thread_specific<std::stack<int>> lineSegmentQueryStackTLS;

  pgo::Mesh::TriMeshBVTree sphereMeshBVTree;
  sphereMeshBVTree.buildByInertiaPartition(sphereMesh);

  pgo::Mesh::TriMeshPseudoNormal sphereMeshNormal;
  sphereMeshNormal.buildPseudoNormals(sphereMesh);

  ES::V3d sphereCenter(center);
  // tbb::parallel_for(0, targetMesh.numVertices(), [&](int vi) {
  for (int vi = 0; vi < targetMesh.numVertices(); vi++) {
    ES::V3d p = targetMesh.pos(vi);
    auto &closestDistDueryStack = closestDistDueryStackTLS.local();

    // close to the sphere
    closestDistDueryStack.clear();
    auto ret = sphereMeshBVTree.closestTriangleQuery(sphereMesh, p, closestDistDueryStack);
    if (ret.dist2 < maxConfidentDist * maxConfidentDist) {
      targetMeshChoosenVertex[vi] = 1;
    }

    // in contact with the sphere
    ES::V3d dir = (ret.closestPosition - p).normalized();
    ES::V3d n = sphereMeshNormal.getPseudoNormal(sphereMesh.triangles().data(), ret.triID, ret.feature);
    if (dir.dot(n) >= 0) {
      targetMeshChoosenVertex[vi] = 1;
    }

    // dist to sphere
    ES::V3d segStart = sphereCenter;
    ES::V3d segEnd = p;

    double segLength = (segEnd - segStart).norm();
    ES::V3d segEndSelf = (segEnd - segStart) / segLength * (segLength - 1e-4) + segStart;

    auto &lineSegmentQueryStack = lineSegmentQueryStackTLS.local();
    while (!lineSegmentQueryStack.empty())
      lineSegmentQueryStack.pop();

    double segW[2];
    int retID = targetMeshBVTree.lineSegmentFirstIntersectionPointSemiExact(targetMesh, segStart, segEndSelf, segW, nullptr, &lineSegmentQueryStack);
    // does not pass through itself
    if (retID < 0) {
      while (!lineSegmentQueryStack.empty())
        lineSegmentQueryStack.pop();

      // intersect with sphere mesh
      int retID = sphereMeshBVTree.lineSegmentFirstIntersectionPointSemiExact(sphereMesh, segStart, segEnd, segW, nullptr, &lineSegmentQueryStack);
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
      retID = sphereMeshBVTree.lineSegmentFirstIntersectionPointSemiExact(sphereMesh, segStart, segEnd, segW, nullptr, &lineSegmentQueryStack);

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
    }
    else {
      targetMeshChoosenVertex[vi] = 0;
    }
  }  //);

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

  // for each triangle edge
  // find necessary edge
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
  struct Array3Less
  {
    bool operator()(const std::array<int, 3> &i0, const std::array<int, 3> &i1) const
    {
      return std::memcmp(i0.data(), i1.data(), sizeof(i0)) < 0;
    }
  };

  std::map<std::array<int, 3>, int, Array3Less> selectedTrianglesQueryByEdgeIDs;
  for (const auto &pr : selectedTrianglesWithEdges) {
    int triID = pr.first;
    std::array<int, 3> edgeIDs = { pr.second[0]->second, pr.second[1]->second, pr.second[2]->second };
    std::sort(edgeIDs.begin(), edgeIDs.end());

    auto ret = selectedTrianglesQueryByEdgeIDs.emplace(edgeIDs, triID);
    PGO_ALOG(ret.second == true);
  }

  std::map<int, int> selectedVertices;
  for (const auto &pr : selectedEdges) {
    selectedVertices.emplace(pr.first, 0);
    selectedVertices.emplace(pr.second, 1);
  }

  std::chrono::high_resolution_clock::time_point t1 = std::chrono::high_resolution_clock::now();

  // create target mesh
  // construct target mesh BV tree
  K::Point_3 centerER(center[0], center[1], center[2]);
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

  // build cutting mesh
  int inc = 0;
  std::vector<K::Point_3> cuttingMeshPoints;
  std::vector<K::Triangle_3> cuttingMeshTriangles;
  K::Point_3 sphereCenterER(sphereCenter[0], sphereCenter[1], sphereCenter[2]);
  for (auto &pr : selectedVertices) {
    pr.second = inc++;

    K::Vector_3 dir = targetMeshVerticesER[pr.first] - sphereCenterER;
    K::Point_3 p1 = sphereCenterER + dir * 2;
    cuttingMeshPoints.emplace_back(p1);
  }
  cuttingMeshPoints.emplace_back(sphereCenterER);

  pgo::Mesh::TriMeshGeo ctri;
  for (int ei = 0; ei < (int)selectedEdges.size(); ei++) {
    ES::V3i triVtxID(selectedVertices[selectedEdges[ei].first], selectedVertices[selectedEdges[ei].second], (int)cuttingMeshPoints.size() - 1);
    cuttingMeshTriangles.emplace_back(cuttingMeshPoints[triVtxID[0]], cuttingMeshPoints[triVtxID[1]], cuttingMeshPoints[triVtxID[2]]);

    ctri.addMesh(pgo::Mesh::createSingleTriangleMesh(toVec3(cuttingMeshPoints[triVtxID[0]]), toVec3(cuttingMeshPoints[triVtxID[1]]), toVec3(cuttingMeshPoints[triVtxID[2]])));
  }

  if (dumpMesh) {
    ctri.save("ctri.obj");
  }

  std::vector<K::Point_3> sphereMeshVerticesER(sphereMesh.numVertices());
  tbb::parallel_for(
    0, sphereMesh.numVertices(), [&](int vi) {
      const ES::V3d &v_d = sphereMesh.pos(vi);
      sphereMeshVerticesER[vi] = K::Point_3(v_d[0], v_d[1], v_d[2]);
    },
    tbb::static_partitioner());

  std::vector<Triangle> sphereMeshTriangles;
  for (int i = 0; i < sphereMesh.numTriangles(); i++) {
    const K::Point_3 &a = sphereMeshVerticesER[sphereMesh.triVtxID(i, 0)];
    const K::Point_3 &b = sphereMeshVerticesER[sphereMesh.triVtxID(i, 1)];
    const K::Point_3 &c = sphereMeshVerticesER[sphereMesh.triVtxID(i, 2)];
    sphereMeshTriangles.emplace_back(a, b, c);
  }

  TriangleTree sphereMeshBVTreeExact(sphereMeshTriangles.begin(), sphereMeshTriangles.end());

  // for each edge find all intersected triangles
  std::vector<std::vector<TriangleTree::Primitive_id>> edgeIntersectedTriangles(selectedEdges.size());
  tbb::parallel_for(0, (int)selectedEdges.size(), [&](int ei) {
    sphereMeshBVTreeExact.all_intersected_primitives(cuttingMeshTriangles[ei], std::back_inserter(edgeIntersectedTriangles[ei]));
  });

  // for each intersected triangle, gather intersections
  std::vector<std::vector<int>> triangleIntersectedEdges(sphereMesh.numTriangles());
  for (int ei = 0; ei < (int)selectedEdges.size(); ei++) {
    for (int i = 0; i < (int)edgeIntersectedTriangles[ei].size(); i++) {
      int tid = (int)(edgeIntersectedTriangles[ei][i] - sphereMeshTriangles.begin());
      triangleIntersectedEdges[tid].emplace_back(ei);
    }
  }

  // =====================================
  // for each intersected triangle, compute corefined triangles
  // we use cutting mesh to cut each sphere triangle
  std::vector<Patch3D> corefinedPatches(sphereMesh.numTriangles());
  tbb::parallel_for(0, sphereMesh.numTriangles(), [&](int ti) {
    // for (int ti = 0; ti < sphereMesh.numTriangles(); ti++) {
    const K::Point_3 &v0 = sphereMeshVerticesER[sphereMesh.triVtxID(ti, 0)];
    const K::Point_3 &v1 = sphereMeshVerticesER[sphereMesh.triVtxID(ti, 1)];
    const K::Point_3 &v2 = sphereMeshVerticesER[sphereMesh.triVtxID(ti, 2)];
    corefinedPatches[ti].rawTriangles.emplace_back(v0, v1, v2);

    if (triangleIntersectedEdges[ti].size() == 0ull) {
      corefinedPatches[ti].vertices = { v0, v1, v2 };
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
        const K::Triangle_3 &tri = cuttingMeshTriangles[edgeID];

        auto iter = corefinedPatches[ti].rawTriangles.begin();
        size_t count = corefinedPatches[ti].rawTriangles.size();

        for (size_t ci = 0; ci < count; ci++) {
          const auto result = CGAL::intersection(tri, *iter);
          if (result) {
            if (const K::Segment_3 *s = boost::get<K::Segment_3>(&*result)) {
              // create 3 triangles for the first point,
              // the segment must be in one of the triangle
              K::Triangle_3 remainingTri;
              bool found = false;
              for (int i = 0; i < 3; i++) {
                K::Triangle_3 tri((*s)[0], (*iter)[i], (*iter)[(i + 1) % 3]);
                if (tri.is_degenerate()) {
                  continue;
                }

                if (tri.has_on((*s)[1]) && found == false) {
                  remainingTri = tri;
                  found = true;
                }
                else {
                  corefinedPatches[ti].rawTriangles.emplace_back(tri);
                }
              }
              PGO_ALOG(remainingTri.squared_area() > 0);

              // create 3 triangles for the second point in the remaining triangle
              // the segment must be on the edge of some triangles
              for (int i = 0; i < 3; i++) {
                K::Triangle_3 tri((*s)[1], remainingTri[i], remainingTri[(i + 1) % 3]);
                if (!tri.is_degenerate()) {
                  corefinedPatches[ti].rawTriangles.emplace_back(tri);
                }
              }
            }
            else if (const K::Point_3 *p = boost::get<K::Point_3>(&*result)) {
              // create 3 triangles
              for (int i = 0; i < 3; i++) {
                K::Triangle_3 tri(*p, (*iter)[i], (*iter)[(i + 1) % 3]);

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

      K::Segment_3 e01(v0, v1), e12(v1, v2), e20(v2, v0);
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

        if (corefinedPatches[ti].vertices[vi] == v0) {
          corefinedPatches[ti].cornerVertices[0] = vi;
        }
        else if (corefinedPatches[ti].vertices[vi] == v1) {
          corefinedPatches[ti].cornerVertices[1] = vi;
        }
        else if (corefinedPatches[ti].vertices[vi] == v2) {
          corefinedPatches[ti].cornerVertices[2] = vi;
        }
      }

      corefinedPatches[ti].triangleEdgeIDs.reserve((int)corefinedPatches[ti].triangles.size());
      for (int tii = 0; tii < (int)corefinedPatches[ti].triangles.size(); tii++) {
        for (int j = 0; j < 3; j++) {
          const K::Point_3 &p0 = corefinedPatches[ti].vertices[corefinedPatches[ti].triangles[tii][j]];
          const K::Point_3 &p1 = corefinedPatches[ti].vertices[corefinedPatches[ti].triangles[tii][(j + 1) % 3]];

          for (int ei = 0; ei < (int)triangleIntersectedEdges[ti].size(); ei++) {
            if (cuttingMeshTriangles[triangleIntersectedEdges[ti][ei]].has_on(p0) &&
              cuttingMeshTriangles[triangleIntersectedEdges[ti][ei]].has_on(p1)) {
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
          K::Segment_3 seg0(cuttingMeshTriangles[edgeID][0], cuttingMeshTriangles[edgeID][2]);
          K::Segment_3 seg1(cuttingMeshTriangles[edgeID][1], cuttingMeshTriangles[edgeID][2]);

          if (seg0.has_on(corefinedPatches[ti].vertices[vi])) {
            corefinedPatches[ti].vertexIsOnTargetMesh[vi][0] = 0;
            corefinedPatches[ti].vertexIsOnTargetMesh[vi][1] = selectedEdges[edgeID].first;
          }
          else if (seg1.has_on(corefinedPatches[ti].vertices[vi])) {
            corefinedPatches[ti].vertexIsOnTargetMesh[vi][0] = 0;
            corefinedPatches[ti].vertexIsOnTargetMesh[vi][1] = selectedEdges[edgeID].second;
          }
          else if (cuttingMeshTriangles[edgeID].has_on(corefinedPatches[ti].vertices[vi])) {
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

  std::chrono::high_resolution_clock::time_point t2 = std::chrono::high_resolution_clock::now();
  SPDLOG_LOGGER_INFO(pgo::Logging::lgr(), "S1 Done. Time cost: {:.4f}s", double(std::chrono::duration_cast<std::chrono::microseconds>(t2 - t1).count()) * 1e-6);

  // ================================================================
  // next we merge all cut triangles into a closed manifold mesh
  pgo::Mesh::TriMeshNeighbor sphereMeshNeighbor(sphereMesh);

  std::vector<int> triIDQ;
  triIDQ.reserve(sphereMesh.numTriangles() * 2);

  std::vector<int> isTriVisited(sphereMesh.numTriangles(), 0);

  std::vector<K::Point_3> finalPoints;
  finalPoints.reserve(sphereMesh.numTriangles() * 10);

  std::vector<ES::V3i> finalTriangles;
  finalTriangles.reserve(sphereMesh.numTriangles() * 10);

  std::map<std::pair<int, int>, std::vector<int>> edgeVertexIDs;
  std::vector<int> newVertexIDs(sphereMesh.numVertices(), -1);

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
      std::pair<int, int> edgeID{ sphereMesh.triVtxID(triIDCur, j), sphereMesh.triVtxID(triIDCur, (j + 1) % 3) };
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
        PGO_ALOG(eitt->second.size() == corefinedPatches[triIDCur].edgeVertices[j].size());
        // find global vertex ID
        for (int k = 0; k < (int)corefinedPatches[triIDCur].edgeVertices[j].size(); k++) {
          if (vertexMasks[corefinedPatches[triIDCur].edgeVertices[j][k]] >= 0)
            continue;

          bool found = false;
          for (int r = 0; r < (int)eitt->second.size(); r++) {
            if (corefinedPatches[triIDCur].vertices[corefinedPatches[triIDCur].edgeVertices[j][k]] == finalPoints[eitt->second[r]]) {
              vertexMasks[corefinedPatches[triIDCur].edgeVertices[j][k]] = eitt->second[r];
              found = true;
              break;
            }
          }
          PGO_ALOG(found == true);
        }
      }
    }  // end for j \in [0, 1, 2]

    // address original triangle corners
    for (int j = 0; j < 3; j++) {
      int vi = corefinedPatches[triIDCur].cornerVertices[j];
      if (newVertexIDs[sphereMesh.triVtxID(triIDCur, j)] >= 0) {
        if (vertexMasks[vi] >= 0) {
          PGO_ALOG(vertexMasks[vi] == newVertexIDs[sphereMesh.triVtxID(triIDCur, j)]);
        }
        else {
          vertexMasks[vi] = newVertexIDs[sphereMesh.triVtxID(triIDCur, j)];
        }
      }
    }

    // compute all vertex ID,
    // especially for the ones that has not been visited
    for (int j = 0; j < (int)corefinedPatches[triIDCur].vertices.size(); j++) {
      if (vertexMasks[j] < 0) {
        finalPoints.emplace_back(corefinedPatches[triIDCur].vertices[j]);
        vertexMasks[j] = (int)finalPoints.size() - 1;
      }
    }

    // address original triangle corners again
    for (int j = 0; j < 3; j++) {
      int vi = corefinedPatches[triIDCur].cornerVertices[j];
      if (newVertexIDs[sphereMesh.triVtxID(triIDCur, j)] < 0) {
        newVertexIDs[sphereMesh.triVtxID(triIDCur, j)] = vertexMasks[vi];
      }
    }

    // adding delete boundary info and edge info
    for (int j = 0; j < 3; j++) {
      std::pair<int, int> edgeID{ sphereMesh.triVtxID(triIDCur, j), sphereMesh.triVtxID(triIDCur, (j + 1) % 3) };
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

      // if (corefinedPatches[triIDCur].vertexIsOnTargetMesh[j] >= 0) {
      //   vertexIsOnTargetMesh[vertexMasks[j]] = corefinedPatches[triIDCur].vertexIsOnTargetMesh[j];
      // }
    }

    // expand bfs search
    ES::V3i triN = sphereMeshNeighbor.getTriangleNeighbors(triIDCur);
    for (int j = 0; j < 3; j++) {
      if (isTriVisited[triN[j]])
        continue;

      triIDQ.emplace_back(triN[j]);
      isTriVisited[triN[j]] = 1;
    }
  }

  t2 = std::chrono::high_resolution_clock::now();
  SPDLOG_LOGGER_INFO(pgo::Logging::lgr(), "S2 Done. Time cost: {:.4f}s", double(std::chrono::duration_cast<std::chrono::microseconds>(t2 - t1).count()) * 1e-6);

  if (dumpMesh) {
    Poly finalMesh;
    pgo::CGALInterface::PolyhedronBuilderNonTriangle<K, std::vector<K::Point_3>::iterator, std::vector<ES::V3i>::iterator> builder(
      finalPoints.begin(), finalTriangles.begin(), (int)finalPoints.size(), (int)finalTriangles.size());
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
  dummyCorefinedMesh.positions().assign(finalPoints.size(), pgo::asVec3d(0.0));
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
        SPDLOG_LOGGER_DEBUG(Logger::lgr(), "Find a false positive edge set. ({}, {}, {})", edgeIDSet[0], edgeIDSet[1], edgeIDSet[2]);
      }
    }

    setID++;
  }

  if (dumpMesh) {
    std::vector<ES::V3d> kk;
    for (int i = 0; i < (int)finalPoints.size(); i++) {
      kk.emplace_back(CGAL::to_double(finalPoints[i].x()),
        CGAL::to_double(finalPoints[i].y()),
        CGAL::to_double(finalPoints[i].z()));
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

    asas.save(fmt::format("rzr00.obj"));
  }

  std::vector<int> setIDAll(setID, 0);
  std::iota(setIDAll.begin(), setIDAll.end(), 0);

  // compute the set ID of triangles that cannot be edited
  std::vector<int> invalidSetIDs;
  std::set_difference(setIDAll.begin(), setIDAll.end(),
    validSetIDs.begin(), validSetIDs.end(),
    std::back_inserter(invalidSetIDs));

  // compute the vertices that cannot be edited
  std::vector<int> unremovableVertex(finalPoints.size(), 0);
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
    std::array<int, 3> cornerVertexIDs;
    for (int vid : edgeLoops) {
      auto it = vertexIsOnTargetMesh.find(vid);
      if (it != vertexIsOnTargetMesh.end() &&
        it->second[0] == 0) {
        PGO_ALOG(ctr < 3);
        cornerVertexIDs[ctr] = vid;
        ctr++;
      }
    }

    // for each edge
    int ith = 0;
    edgeLoops1.clear();
    for (int j = 0; j < 3; j++) {
      while (edgeLoops[ith] != cornerVertexIDs[j]) {
        ith = (ith + 1) % (int)edgeLoops.size();
      }

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

    finalFaces.emplace_back(edgeLoops1);
    if (edgeLoops1.size() == 3ull) {
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
  finalPoints1.reserve(finalPoints.size());

  inc = 0;
  for (auto &v : usedVertices) {
    v.second = inc++;
    finalPoints1.emplace_back(finalPoints[v.first]);
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

      K::Point_3 triCenter = h0->vertex()->point() + e01 / 3.0 + e02 / 3.0;
      K::Vector_3 dir = triCenter - sphereCenterER;

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
  }

  if (triangulateMesh(finalMesh) == false) {
    std::cout << "First iteration cannot clear the non triangles." << std::endl;

    for (auto fi = finalMesh.facets_begin(); fi != finalMesh.facets_end(); ++fi) {
      if (fi->is_triangle())
        continue;

      Poly::Halfedge_handle h0 = fi->halfedge();
      Poly::Halfedge_handle h1 = h0;
      K::FT newCenter[3] = { 0, 0, 0 };
      int c = 0;
      do {
        const K::Point_3 &pt = h1->opposite()->vertex()->point();
        newCenter[0] += pt[0];
        newCenter[1] += pt[1];
        newCenter[2] += pt[2];

        h1 = h1->next()->opposite();
        c += 1;
      } while (h0 != h1);

      int vidCur = (int)finalMesh.size_of_vertices();
      Poly::Halfedge_handle h = finalMesh.create_center_vertex(fi->halfedge());
      h->vertex()->point() = K::Point_3(newCenter[0] / (double)c, newCenter[1] / (double)c, newCenter[2] / (double)c);
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

  if constexpr (dumpMesh) {
    std::ofstream("rzr2.off") << finalMesh;
  }

  // project vertex to the best location
  vertexProjection(finalMesh, targetMesh, targetMeshBVTreeExact,
    centerER, maxAllowedRadiusDist,
    vtxIDNew2Old, vertexIsOnTargetMesh, selectedEdges);

  t2 = std::chrono::high_resolution_clock::now();
  SPDLOG_LOGGER_INFO(pgo::Logging::lgr(), "Done. Time cost: {:.4f}s", double(std::chrono::duration_cast<std::chrono::microseconds>(t2 - t1).count()) * 1e-6);

  if constexpr (dumpMesh) {
    std::ofstream("rzr3.off") << finalMesh;
  }

  std::vector<VertexProperty> finalMeshVertexProperties(finalMesh.size_of_vertices());
  for (auto vit = finalMesh.vertices_begin(); vit != finalMesh.vertices_end(); ++vit) {
    int vid = (int)vit->id();

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

  cleanUpMesh(targetMeshPoly, targetMeshBVTreeExact, finalMesh, 5e-3, finalMeshVertexProperties, targetEdges);

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

      K::FT p[3] = { 0, 0, 0 };
      int n = 0;
      do {
        const K::Point_3 &pt = h1->opposite()->vertex()->point();
        for (int j = 0; j < 3; j++) {
          p[j] += pt[j];
        }

        h1 = h1->next()->opposite();
        n++;
      } while (h != h1);

      K::Point_3 newPt(p[0] / n, p[1] / n, p[2] / n);
      K::Vector_3 diff = newPt - vit->point();
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

int MedialAxisRepresentation::remeshSphereMeshWithTarget(const pgo::Mesh::TriMeshGeo &sphereMesh, const double center[3],
  const pgo::Mesh::TriMeshGeo &targetMesh, [[maybe_unused]] const pgo::Mesh::TriMeshBVTree &targetMeshBVTree, [[maybe_unused]] const pgo::Mesh::TriMeshPseudoNormal &targetMeshNormals,
  [[maybe_unused]] double maxConfidentDist, [[maybe_unused]] double maxAllowedRadiusDist, pgo::Mesh::TriMeshGeo &meshOut)
{
  using K = pgo::CGALInterface::KernelExact;
  // using K = CGAL::Simple_cartesian<mpq_class>;
  using Poly = pgo::CGALInterface::Polyhedron<K>;

  // double maxConfidentDist = 1e-3;
  // double maxAllowedRadiusDist = 0.1;
  int dumpMesh = 0;

  pgo::Mesh::TriMeshGeo sphereMeshEnlarged = sphereMesh;
  ES::V3d sphereCenter(center[0], center[1], center[2]);
  tbb::parallel_for(0, sphereMeshEnlarged.numVertices(), [&](int vi) {
    ES::V3d dir = sphereMeshEnlarged.pos(vi) - sphereCenter;
    double dirLen = dir.norm();
    sphereMeshEnlarged.pos(vi) = dir / dirLen * (dirLen + maxAllowedRadiusDist) + sphereCenter;
  });

  if (dumpMesh) {
    sphereMeshEnlarged.save("zra.obj");
  }

  Poly sphereMeshEnlargedPoly, targetMeshPoly, intersectedMeshPoly;
  pgo::CGALInterface::triangleMesh2Polyhedron<K>(sphereMeshEnlarged, sphereMeshEnlargedPoly);
  pgo::CGALInterface::triangleMesh2Polyhedron<K>(targetMesh, targetMeshPoly);
  CGAL::Polygon_mesh_processing::corefine_and_compute_intersection(sphereMeshEnlargedPoly, targetMeshPoly, intersectedMeshPoly);

  pgo::CGALInterface::polyhedron2TriangleMesh(intersectedMeshPoly, meshOut);
  std::ofstream("zraa.off") << intersectedMeshPoly;

  return 0;
}
