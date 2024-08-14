#pragma once

#include "CGALTemplateUtilities.h"

#include "EigenSupport.h"
#include "triKey.h"

#include <vector>
#include <unordered_map>

namespace MedialAxisRepresentation
{
using EK = CGAL::Simple_cartesian<CGAL::Gmpq>;

struct SkeletonExact
{
  struct Vertex
  {
    EK::Point_3 pos;
    std::vector<int> neighboringEdges;
    std::vector<int> neighboringTriangles;

    bool isBoundary() const { return neighboringTriangles.size() == 0ull && neighboringEdges.size() == 1ull; }
  };

  struct Edge
  {
    pgo::EigenSupport::V2i vidx;
    std::vector<int> neighboringTriangles;

    bool isBoundary() const { return neighboringTriangles.size() == 1ull; }
    bool isTriangleEdge() const { return neighboringTriangles.size() != 0ull; }
  };

  struct Triangle
  {
    pgo::EigenSupport::V3i vidx;
    pgo::EigenSupport::V3i eidx;
  };

public:
  // map
  std::unordered_map<int, Vertex> vertices;
  std::unordered_map<int, Edge> edges;
  std::unordered_map<int, Triangle> triangles;

  // inverse map
  std::unordered_map<std::pair<int, int>, int, pgo::EigenSupport::IntPairHash, pgo::EigenSupport::IntPairEqual> edgeVtxIDToEdgeID;
  std::unordered_map<pgo::Mesh::UTriKey, int> triVtxIDToTriID;

  // id counters
  int vid_counter = 0;
  int eid_counter = 0;
  int tid_counter = 0;

public:
  const Vertex &getV(int vid) const;
  Vertex &getV(int vid);

  const Edge &getE(int eid) const;
  Edge &getE(int eid);

  const Triangle &getT(int tid) const;
  Triangle &getT(int tid);

  void reset();

  void save(const char *filename) const;  // No VertexLocalProfile
  int load(const char *filename);        // No VertexLocalProfile

  void initFromInput(const std::vector<EK::Point_3> &inputVertices, const std::vector<std::pair<int, int>> &inputEdges, const std::vector<std::tuple<int, int, int>> &inputTris);
  void initFromInput(const std::vector<EK::Point_3> &inputVertices, const std::vector<std::vector<int>> &inputFaceSets);
  void exportToArray(
    std::vector<EK::Point_3> &inputVertices, 
    std::vector<std::pair<int, int>> &inputEdges, 
    std::vector<std::tuple<int, int, int>> &inputTris,
    std::vector<int> *vertexIDs = nullptr);
  
  void exportToArray(std::vector<EK::Point_3> &inputVertices, std::vector<std::pair<int, int>> &inputEdges, std::vector<std::tuple<int, int, int>> &inputTris, 
                    std::map<int, int> &vidMap, std::map<int, std::pair<int, int>> &edgeMap, std::map<int, std::tuple<int, int, int>> &triMap);
  void gatherVertexIDs(std::vector<int> &vertexIDs) const;
  void gatherTriangleIDs(std::vector<int> &triIDs) const;
  void gatherEdgeIDs(std::vector<int> &edgeIDs) const;

  void inducedSubgraph(const std::vector<int> &vertexIDs, const int &vIDOffset,
    std::map<int, int> &maVID2VID, std::vector<std::vector<int>> &edgeIDs, std::vector<std::vector<int>> &triangleIDs);

  bool canCollapseEdge(int) const;
  int collapseEdge(int eid, const EK::Point_3 &targetPos);
  int addVertexToEdge(int eid, const EK::Point_3 &targetPos);

  void saveDisplayMesh(const char *filename) const;
  bool isConnected() const;
  bool isConsistent() const;

protected:
  void eraseTriangle(int triID, int deleteMap = 1);
  void replaceVertex(int oldVID, const Vertex &oldVertex, int newVID);
};

inline const SkeletonExact::Vertex &SkeletonExact::getV(int vid) const
{
  auto iter = vertices.find(vid);
  if (iter != vertices.end()) {
    return iter->second;
  }
  else {
    throw std::domain_error("Cannot find vertex");
  }
}

inline SkeletonExact::Vertex &SkeletonExact::getV(int vid)
{
  auto iter = vertices.find(vid);
  if (iter != vertices.end()) {
    return iter->second;
  }
  else {
    throw std::domain_error("Cannot find vertex");
  }
}

inline const SkeletonExact::Edge &SkeletonExact::getE(int eid) const
{
  auto iter = edges.find(eid);
  if (iter != edges.end()) {
    return iter->second;
  }
  else {
    throw std::domain_error("Cannot find edge");
  }
}

inline SkeletonExact::Edge &SkeletonExact::getE(int eid)
{
  auto iter = edges.find(eid);
  if (iter != edges.end()) {
    return iter->second;
  }
  else {
    throw std::domain_error("Cannot find edge");
  }
}

inline const SkeletonExact::Triangle &SkeletonExact::getT(int tid) const
{
  auto iter = triangles.find(tid);
  if (iter != triangles.end()) {
    return iter->second;
  }
  else {
    throw std::domain_error("Cannot find triangle");
  }
}

inline SkeletonExact::Triangle &SkeletonExact::getT(int tid)
{
  auto iter = triangles.find(tid);
  if (iter != triangles.end()) {
    return iter->second;
  }
  else {
    throw std::domain_error("Cannot find triangle");
  }
}

inline void SkeletonExact::gatherVertexIDs(std::vector<int> &vertexIDs) const
{
  vertexIDs.clear();
  vertexIDs.reserve(vertices.size());
  for (const auto &pr_v : vertices)
    vertexIDs.emplace_back(pr_v.first);
}

inline void SkeletonExact::gatherTriangleIDs(std::vector<int> &triIDs) const
{
  triIDs.reserve(triangles.size());
  for (const auto &tri : triangles) {
    triIDs.push_back(tri.first);
  }
}

inline void SkeletonExact::gatherEdgeIDs(std::vector<int> &edgeIDs) const
{
  edgeIDs.reserve(edges.size());
  for (const auto &e : edges) {
    edgeIDs.push_back(e.first);
  }
}

// #define CCK fmt::print("Consistency check: {}\n", isConsistent())
#define CCK
}  // namespace MedialAxisRepresentation