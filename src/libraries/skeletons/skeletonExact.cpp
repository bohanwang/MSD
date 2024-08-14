#include "skeletonExact.h"

#include "pgoLogging.h"
#include "basicAlgorithms.h"
#include "createTriMesh.h"

#include <nlohmann/json.hpp>

#include <fmt/format.h>

namespace MedialAxisRepresentation
{
namespace ES = pgo::EigenSupport;
inline ES::V3d toVec3(const EK::Point_3 &pt)
{
  return ES::V3d(CGAL::to_double(pt.x()), CGAL::to_double(pt.y()), CGAL::to_double(pt.z()));
}

inline ES::V3d toVec3(const EK::Vector_3 &pt)
{
  return ES::V3d(CGAL::to_double(pt[0]), CGAL::to_double(pt[1]), CGAL::to_double(pt[2]));
}

void to_json(nlohmann::json &j, const SkeletonExact::Vertex &vtx)
{
  std::vector<char> val[3];
  for (int j = 0; j < 3; j++) {
    auto &op = vtx.pos[j].mpq();
    val[j].resize((mpz_sizeinbase(mpq_numref(op), 10) + mpz_sizeinbase(mpq_denref(op), 10) + 3) * 2);
    memset(val[j].data(), 0, sizeof(char) * val[j].size());
    mpq_get_str(val[j].data(), 10, op);
  }

  j = nlohmann::json{
    { "pos", { val[0], val[1], val[2] } },
    { "neighboring-edges", vtx.neighboringEdges },
    { "neighboring-triangles", vtx.neighboringTriangles }
  };
}

void to_json(nlohmann::json &j, const SkeletonExact::Edge &e)
{
  std::array<int, 2> vidx{ e.vidx[0], e.vidx[1] };
  j = nlohmann::json{
    { "vidx", vidx },
    { "neighboring-triangles", e.neighboringTriangles }
  };
}

void to_json(nlohmann::json &j, const SkeletonExact::Triangle &tri)
{
  std::array<int, 3> vidx{ tri.vidx[0], tri.vidx[1], tri.vidx[2] };
  std::array<int, 3> eidx{ tri.eidx[0], tri.eidx[1], tri.eidx[2] };
  j = nlohmann::json{
    { "vidx", vidx },
    { "eidx", eidx }
  };
}

void to_json(nlohmann::json &j, const SkeletonExact &ma)
{
  std::unordered_map<std::string, int> edgeVtxIDToEdgeID;
  edgeVtxIDToEdgeID.reserve(ma.edgeVtxIDToEdgeID.size());

  for (const auto &pr : ma.edgeVtxIDToEdgeID) {
    edgeVtxIDToEdgeID.emplace(
      fmt::format("{}-{}", pr.first.first, pr.first.second),
      pr.second);
  }

  std::unordered_map<std::string, int> triVtxIDToTriID;
  triVtxIDToTriID.reserve(ma.triVtxIDToTriID.size());
  for (const auto &pr : ma.triVtxIDToTriID) {
    triVtxIDToTriID.emplace(
      fmt::format("{}-{}-{}", pr.first[0], pr.first[1], pr.first[2]),
      pr.second);
  }

  j = nlohmann::json{
    { "vertices", ma.vertices },
    { "edges", ma.edges },
    { "triangles", ma.triangles },
    { "edge-vid-to-edge-id", edgeVtxIDToEdgeID },
    { "tri-vid-to-tri-id", triVtxIDToTriID },
    { "vtx-counter", ma.vid_counter },
    { "eid-counter", ma.eid_counter },
    { "tid-counter", ma.tid_counter },
  };
}

void from_json(const nlohmann::json &j, SkeletonExact::Vertex &vtx)
{
  j.at("neighboring-edges").get_to(vtx.neighboringEdges);
  j.at("neighboring-triangles").get_to(vtx.neighboringTriangles);

  std::vector<char> valStr[3];
  // for (auto it = j["pos"].begin(); it != j["pos"].end(); ++it) {
  for (int count = 0; count < 3; count++) {
    valStr[count] = j["pos"][count];
  }
  vtx.pos = EK::Point_3(CGAL::Gmpq(valStr[0].data(), 10),
    CGAL::Gmpq(valStr[1].data(), 10),
    CGAL::Gmpq(valStr[2].data(), 10));
}

void from_json(const nlohmann::json &j, SkeletonExact::Edge &e)
{
  std::array<int, 2> vidx;
  j.at("vidx").get_to(vidx);
  e.vidx[0] = vidx[0];
  e.vidx[1] = vidx[1];

  j.at("neighboring-triangles").get_to(e.neighboringTriangles);
}

void from_json(const nlohmann::json &j, SkeletonExact::Triangle &tri)
{
  std::array<int, 3> vidx;
  j.at("vidx").get_to(vidx);
  tri.vidx[0] = vidx[0];
  tri.vidx[1] = vidx[1];
  tri.vidx[2] = vidx[2];

  std::array<int, 3> eidx;
  j.at("eidx").get_to(eidx);
  tri.eidx[0] = eidx[0];
  tri.eidx[1] = eidx[1];
  tri.eidx[2] = eidx[2];
}

void from_json(const nlohmann::json &j, SkeletonExact &ma)
{
  j.at("vertices").get_to(ma.vertices);
  j.at("edges").get_to(ma.edges);
  j.at("triangles").get_to(ma.triangles);

  std::unordered_map<std::string, int> edgeVtxIDToEdgeID;
  j.at("edge-vid-to-edge-id").get_to(edgeVtxIDToEdgeID);
  ma.edgeVtxIDToEdgeID.reserve(edgeVtxIDToEdgeID.size());
  for (const auto &pr : edgeVtxIDToEdgeID) {
    auto pos = pr.first.find('-');
    ma.edgeVtxIDToEdgeID.emplace(
      std::make_pair(std::stoi(pr.first.substr(0, pos)), std::stoi(pr.first.substr(pos + 1))),
      pr.second);
  }

  std::unordered_map<std::string, int> triVtxIDToTriID;
  j.at("tri-vid-to-tri-id").get_to(triVtxIDToTriID);
  ma.triVtxIDToTriID.reserve(triVtxIDToTriID.size());
  for (const auto &pr : triVtxIDToTriID) {
    auto pos1 = pr.first.find('-');
    auto pos2 = pr.first.find('-', pos1 + 1);
    ma.triVtxIDToTriID.emplace(
      pgo::Mesh::UTriKey(std::stoi(pr.first.substr(0, pos1)),
        std::stoi(pr.first.substr(pos1 + 1, pos2 - pos1 - 1)),
        std::stoi(pr.first.substr(pos2 + 1))),
      pr.second);
  }

  j.at("vtx-counter").get_to(ma.vid_counter);
  j.at("eid-counter").get_to(ma.eid_counter);
  j.at("tid-counter").get_to(ma.tid_counter);
}

}  // namespace MedialAxisRepresentation

void MedialAxisRepresentation::SkeletonExact::reset()
{
  vertices.clear();
  edges.clear();
  triangles.clear();
  edgeVtxIDToEdgeID.clear();
  triVtxIDToTriID.clear();

  vid_counter = 0;
  eid_counter = 0;
  tid_counter = 0;
}

void MedialAxisRepresentation::SkeletonExact::save(const char *filename) const
{
  nlohmann::json jout = *this;
  std::ofstream(filename) << jout.dump(2);
}

int MedialAxisRepresentation::SkeletonExact::load(const char *filename)
{
  std::ifstream ifs(filename);
  if (!ifs) {
    std::cerr << "Cannot open filename " << filename << std::endl;
    return 1;
  }

  nlohmann::json jin;
  ifs >> jin;
  jin.get_to(*this);

  return 0;
}

void MedialAxisRepresentation::SkeletonExact::inducedSubgraph(const std::vector<int> &vertexIDs, const int &vIDOffset,
  std::map<int, int> &maVID2VID, std::vector<std::vector<int>> &edgeIDs, std::vector<std::vector<int>> &triangleIDs)
{
  for (int i = 0; i < vertexIDs.size(); i++) {
    maVID2VID[vertexIDs[i]] = i + vIDOffset;
  }

  for (auto &e : edges) {
    if (maVID2VID.find(e.second.vidx[0]) != maVID2VID.end() && maVID2VID.find(e.second.vidx[1]) != maVID2VID.end()) {
      edgeIDs.emplace_back();
      edgeIDs.back() = { maVID2VID[e.second.vidx[0]], maVID2VID[e.second.vidx[1]] };
    }
  }

  for (auto &t : triangles) {
    if (maVID2VID.find(t.second.vidx[0]) != maVID2VID.end() && maVID2VID.find(t.second.vidx[1]) != maVID2VID.end() &&
      maVID2VID.find(t.second.vidx[2]) != maVID2VID.end()) {
      triangleIDs.emplace_back();
      triangleIDs.back() = { maVID2VID[t.second.vidx[0]], maVID2VID[t.second.vidx[1]], maVID2VID[t.second.vidx[2]] };
    }
  }
}

bool MedialAxisRepresentation::SkeletonExact::canCollapseEdge(int eid) const
{
  struct Link
  {
    std::vector<int> vertexIDs;
    std::vector<int> edgeIDs;
  };

  Link lk_v[2], lk_e;
  const auto &edge = getE(eid);

  for (int i = 0; i < 2; i++) {
    const auto &v = getV(edge.vidx[i]);
    for (int eid1 : v.neighboringEdges) {
      const auto &edge1 = getE(eid1);
      for (int j = 0; j < 2; j++) {
        if (edge1.vidx[j] != edge.vidx[i]) {
          lk_v[i].vertexIDs.emplace_back(edge1.vidx[j]);
        }
      }
    }

    for (int triIdx : v.neighboringTriangles) {
      const auto &tri = getT(triIdx);
      for (int j = 0; j < 3; j++) {
        if (tri.vidx[j] != edge.vidx[i]) {
          lk_v[i].vertexIDs.emplace_back(tri.vidx[j]);
        }
      }

      for (int j = 0; j < 3; j++) {
        const auto &edge1 = getE(tri.eidx[j]);
        if (edge1.vidx[0] != edge.vidx[i] && edge1.vidx[1] != edge.vidx[i]) {
          lk_v[i].edgeIDs.emplace_back(tri.eidx[j]);
        }
      }
    }
  }

  for (int triIdx : edge.neighboringTriangles) {
    const auto &tri = getT(triIdx);
    for (int j = 0; j < 3; j++) {
      if (tri.vidx[j] != edge.vidx[0] && tri.vidx[j] != edge.vidx[1]) {
        lk_e.vertexIDs.emplace_back(tri.vidx[j]);
        break;
      }
    }
  }

  pgo::BasicAlgorithms::sortAndDeduplicateWithErase(lk_v[0].vertexIDs);
  pgo::BasicAlgorithms::sortAndDeduplicateWithErase(lk_v[0].edgeIDs);
  pgo::BasicAlgorithms::sortAndDeduplicateWithErase(lk_v[1].vertexIDs);
  pgo::BasicAlgorithms::sortAndDeduplicateWithErase(lk_v[1].edgeIDs);

  Link lk_v_inter;
  std::set_intersection(lk_v[0].vertexIDs.begin(), lk_v[0].vertexIDs.end(), lk_v[1].vertexIDs.begin(), lk_v[1].vertexIDs.end(), std::back_inserter(lk_v_inter.vertexIDs));
  std::set_intersection(lk_v[0].edgeIDs.begin(), lk_v[0].edgeIDs.end(), lk_v[1].edgeIDs.begin(), lk_v[1].edgeIDs.end(), std::back_inserter(lk_v_inter.edgeIDs));

  pgo::BasicAlgorithms::sortAndDeduplicateWithErase(lk_e.vertexIDs);

  if (lk_v_inter.edgeIDs.size() == lk_e.edgeIDs.size() &&
    lk_v_inter.vertexIDs.size() == lk_e.vertexIDs.size() &&
    std::memcmp(lk_v_inter.vertexIDs.data(), lk_e.vertexIDs.data(), lk_e.vertexIDs.size() * sizeof(int)) == 0) {
    return true;
  }
  else {
    return false;
  }
}

int MedialAxisRepresentation::SkeletonExact::collapseEdge(int eid, const EK::Point_3 &targetPos)
{
  const auto &edge = getE(eid);
  const auto &v0 = getV(edge.vidx[0]);
  const auto &v1 = getV(edge.vidx[1]);

  // if there are neighboring triangles,
  // delete them (triangles with edge eid)
  while (edge.neighboringTriangles.size() > 0) {
    int triID = edge.neighboringTriangles.back();
    eraseTriangle(triID);
  }

  Vertex newVtx;
  newVtx.pos = targetPos;

  // merge triangle neighbors
  newVtx.neighboringTriangles.insert(newVtx.neighboringTriangles.end(), v0.neighboringTriangles.begin(), v0.neighboringTriangles.end());
  newVtx.neighboringTriangles.insert(newVtx.neighboringTriangles.end(), v1.neighboringTriangles.begin(), v1.neighboringTriangles.end());

  if (newVtx.neighboringTriangles.size())
    pgo::BasicAlgorithms::sortAndDeduplicateWithErase(newVtx.neighboringTriangles);

  // delete this edge from its end vertices
  // this edge must be a pure edge after the removal of the triangles
  for (int j = 0; j < 2; j++) {
    auto &vtx = getV(edge.vidx[j]);
    auto it = std::find(vtx.neighboringEdges.begin(), vtx.neighboringEdges.end(), eid);
    PGO_ALOG(it != vtx.neighboringEdges.end());

    vtx.neighboringEdges.erase(it);
  }

  // merge edge neighbors
  newVtx.neighboringEdges.insert(newVtx.neighboringEdges.end(), v0.neighboringEdges.begin(), v0.neighboringEdges.end());
  newVtx.neighboringEdges.insert(newVtx.neighboringEdges.end(), v1.neighboringEdges.begin(), v1.neighboringEdges.end());
  if (newVtx.neighboringEdges.size())
    pgo::BasicAlgorithms::sortAndDeduplicateWithErase(newVtx.neighboringEdges);

  // erase the edge itself
  auto iter_e = edges.find(eid);
  ES::V2i oldVIDs = iter_e->second.vidx;

  // no need to check if found, as the first line of this function has done it
  edgeVtxIDToEdgeID.erase(std::make_pair(oldVIDs[0], oldVIDs[1]));
  edges.erase(iter_e);

  // erase the end vertices and replace them with the new ones
  //
  // insert the new vertex
  vertices.emplace(vid_counter++, newVtx);

  int newVID = vid_counter - 1;
  // replace v0
  replaceVertex(oldVIDs[0], v0, newVID);
  // do the same thing for v1
  replaceVertex(oldVIDs[1], v1, newVID);

  // delete those vertices
  auto iter_v = vertices.find(oldVIDs[0]);
  vertices.erase(iter_v);

  iter_v = vertices.find(oldVIDs[1]);
  vertices.erase(iter_v);

  CCK;

  return newVID;
}

int MedialAxisRepresentation::SkeletonExact::addVertexToEdge(int eid, const EK::Point_3 &targetPos)
{
  const auto &edge = getE(eid);
  // const auto &v0 = getV(edge.vidx[0]);
  // const auto &v1 = getV(edge.vidx[1]);

  std::vector<ES::V3i> deletedTriangles;
  deletedTriangles.reserve(edge.neighboringTriangles.size());

  // if there are neighboring triangles,
  // delete them (triangles with edge eid)
  while (edge.neighboringTriangles.size() > 0) {
    int triID = edge.neighboringTriangles.back();
    deletedTriangles.push_back(getT(triID).vidx);
    eraseTriangle(triID);
  }

  // create new vertex
  Vertex newVtx;
  newVtx.pos = targetPos;

  // delete this edge from its end vertices
  // this edge must be a pure edge after the removal of the triangles
  for (int j = 0; j < 2; j++) {
    auto &vtx = getV(edge.vidx[j]);
    auto it = std::find(vtx.neighboringEdges.begin(), vtx.neighboringEdges.end(), eid);
    PGO_ALOG(it != vtx.neighboringEdges.end());

    vtx.neighboringEdges.erase(it);
  }

  // insert the new vertex
  vertices.emplace(vid_counter++, newVtx);
  int newVID = vid_counter - 1;

  auto iter_e = edges.find(eid);
  ES::V2i oldVIDs = iter_e->second.vidx;

  // create and insert new edges
  Edge newEdge0, newEdge1;
  newEdge0.vidx[0] = oldVIDs[0];
  newEdge0.vidx[1] = newVID;
  if (newEdge0.vidx[0] > newEdge0.vidx[1])
    std::swap(newEdge0.vidx[0], newEdge0.vidx[1]);

  newEdge1.vidx[0] = oldVIDs[1];
  newEdge1.vidx[1] = newVID;
  if (newEdge1.vidx[0] > newEdge1.vidx[1])
    std::swap(newEdge1.vidx[0], newEdge1.vidx[1]);

  int newEID0 = eid_counter++;
  int newEID1 = eid_counter++;

  PGO_ALOG(edges.emplace(newEID0, newEdge0).second == true);
  PGO_ALOG(edges.emplace(newEID1, newEdge1).second == true);
  PGO_ALOG(edgeVtxIDToEdgeID.emplace(std::make_pair(newEdge0.vidx[0], newEdge0.vidx[1]), newEID0).second == true);
  PGO_ALOG(edgeVtxIDToEdgeID.emplace(std::make_pair(newEdge1.vidx[0], newEdge1.vidx[1]), newEID1).second == true);

  // erase the edge itself
  // no need to check if found, as the first line of this function has done it
  edgeVtxIDToEdgeID.erase(std::make_pair(oldVIDs[0], oldVIDs[1]));
  edges.erase(eid);

  // if previously the edge is not a triangle edge
  if (deletedTriangles.size() == 0) {
    getV(oldVIDs[0]).neighboringEdges.emplace_back(newEID0);
    getV(oldVIDs[1]).neighboringEdges.emplace_back(newEID1);

    getV(newVID).neighboringEdges.emplace_back(newEID0);
    getV(newVID).neighboringEdges.emplace_back(newEID1);
  }
  // otherwise, create deleted triangles now
  else {
    for (const auto &tri : deletedTriangles) {
      int vid = 0;
      if (tri[0] != oldVIDs[0] && tri[0] != oldVIDs[1]) {
        vid = tri[0];
      }
      else if (tri[1] != oldVIDs[0] && tri[1] != oldVIDs[1]) {
        vid = tri[1];
      }
      else {
        vid = tri[2];
      }

      ES::V3i triVtxIDs[2] = {
        { oldVIDs[0], vid, newVID },
        { oldVIDs[1], vid, newVID },
      };

      // add one edge
      Edge newEdge2;
      newEdge2.vidx[0] = vid;
      newEdge2.vidx[1] = newVID;
      if (newEdge2.vidx[0] > newEdge2.vidx[1])
        std::swap(newEdge2.vidx[0], newEdge2.vidx[1]);

      int newEID2 = eid_counter++;
      PGO_ALOG(edges.emplace(newEID2, newEdge2).second == true);
      PGO_ALOG(edgeVtxIDToEdgeID.emplace(std::make_pair(newEdge2.vidx[0], newEdge2.vidx[1]), newEID2).second == true);

      // add two triangles
      for (int tt = 0; tt < 2; tt++) {
        Triangle newTri;
        newTri.vidx = triVtxIDs[tt];

        // get edge id
        for (int j = 0; j < 3; j++) {
          std::pair<int, int> evidx{ newTri.vidx[j], newTri.vidx[(j + 1) % 3] };
          if (evidx.first > evidx.second)
            std::swap(evidx.first, evidx.second);
          auto it = edgeVtxIDToEdgeID.find(evidx);
          PGO_ALOG(it != edgeVtxIDToEdgeID.end());

          newTri.eidx[j] = it->second;
        }

        int newTID = tid_counter++;
        PGO_ALOG(triangles.emplace(newTID, newTri).second == true);
        PGO_ALOG(triVtxIDToTriID.emplace(pgo::Mesh::UTriKey(newTri.vidx.data()), newTID).second == true);

        // add to neighbor vertices and edges
        for (int j = 0; j < 3; j++) {
          getV(newTri.vidx[j]).neighboringTriangles.emplace_back(newTID);

          auto &e = getE(newTri.eidx[j]);
          e.neighboringTriangles.emplace_back(newTID);
          if (e.neighboringTriangles.size() > 0ull) {
            for (int z = 0; z < 2; z++) {
              auto &vvv = getV(e.vidx[z]);
              auto it = std::find(vvv.neighboringEdges.begin(), vvv.neighboringEdges.end(), newTri.eidx[j]);
              if (it != vvv.neighboringEdges.end()) {
                vvv.neighboringEdges.erase(it);
              }
            }
          }
        }
      }
    }
  }

  CCK;
  return newVID;
}

void MedialAxisRepresentation::SkeletonExact::initFromInput(
  const std::vector<EK::Point_3> &inputVertices,
  const std::vector<std::pair<int, int>> &inputEdges,
  const std::vector<std::tuple<int, int, int>> &inputTris)
{
  reset();

  for (const auto &vtx : inputVertices) {
    Vertex v;
    v.pos = vtx;
    vertices.emplace(vid_counter, v);
    vid_counter++;
  }

  for (const auto &e : inputEdges) {
    Edge edge;
    edge.vidx[0] = e.first;
    edge.vidx[1] = e.second;

    if (edge.vidx[0] > edge.vidx[1]) {
      std::swap(edge.vidx[0], edge.vidx[1]);
    }

    edges.emplace(eid_counter, edge);

    for (int i = 0; i < 2; i++) {
      auto &v = getV(edge.vidx[i]);
      v.neighboringEdges.emplace_back(eid_counter);
      pgo::BasicAlgorithms::sortAndDeduplicateWithErase(v.neighboringEdges);
    }

    edgeVtxIDToEdgeID.emplace(std::pair{ edge.vidx[0], edge.vidx[1] }, eid_counter);

    eid_counter++;
  }

  // first add edges of the triangles
  for (auto &t : inputTris) {
    std::array<int, 3> triVtxIDs = { std::get<0>(t), std::get<1>(t), std::get<2>(t) };

    for (int i = 0; i < 3; i++) {
      int v0 = triVtxIDs[i], v1 = triVtxIDs[(i + 1) % 3];
      if (v0 > v1)
        std::swap(v0, v1);
      auto it = edgeVtxIDToEdgeID.find(std::pair<int, int>{ v0, v1 });
      if (it == edgeVtxIDToEdgeID.end()) {
        // add new edge
        Edge edge;
        edge.vidx[0] = v0;
        edge.vidx[1] = v1;
        edges.emplace(eid_counter, edge);

        // for (int ei = 0; ei < 2; ei++) {
        //   auto &v = getV(edge.vidx[ei]);
        //   v.neighboringEdges.emplace_back(eid_counter);
        //   sortAndDeduplicateWithErase(v.neighboringEdges);
        // }

        PGO_ALOG(v0 < v1);
        edgeVtxIDToEdgeID.emplace(std::pair{ v0, v1 }, eid_counter);

        eid_counter++;
      }
    }
  }

  // then add triangles
  for (auto &t : inputTris) {
    Triangle tri;
    tri.vidx = { std::get<0>(t), std::get<1>(t), std::get<2>(t) };

    // get edge id
    for (int i = 0; i < 3; i++) {
      int v0 = tri.vidx[i], v1 = tri.vidx[(i + 1) % 3];
      if (v0 > v1)
        std::swap(v0, v1);
      auto it = edgeVtxIDToEdgeID.find(std::pair<int, int>{ v0, v1 });
      PGO_ALOG(it != edgeVtxIDToEdgeID.end());

      tri.eidx[i] = it->second;
    }

    triangles.emplace(tid_counter, tri);
    triVtxIDToTriID.emplace(pgo::Mesh::UTriKey(tri.vidx.data()), tid_counter);

    // add to neighbor vertices and edges
    for (int i = 0; i < 3; i++) {
      auto &v = getV(tri.vidx[i]);
      v.neighboringTriangles.emplace_back(tid_counter);
      pgo::BasicAlgorithms::sortAndDeduplicateWithErase(v.neighboringTriangles);

      auto &e = getE(tri.eidx[i]);
      e.neighboringTriangles.emplace_back(tid_counter);
      pgo::BasicAlgorithms::sortAndDeduplicateWithErase(e.neighboringTriangles);
    }
    tid_counter++;
  }
}

void MedialAxisRepresentation::SkeletonExact::initFromInput(const std::vector<EK::Point_3> &inputVertices, const std::vector<std::vector<int>> &inputFaceSets)
{
  std::vector<std::pair<int, int>> edges;
  std::vector<std::tuple<int, int, int>> triangles;
  for (auto f : inputFaceSets) {
    PGO_ALOG((int)f.size() == 3 || (int)f.size() == 2);
    if ((int)f.size() == 3) {
      triangles.emplace_back(f[0], f[1], f[2]);
    }
    else if ((int)f.size() == 2) {
      edges.emplace_back(f[0], f[1]);
    }
  }
  initFromInput(inputVertices, edges, triangles);
}

void MedialAxisRepresentation::SkeletonExact::exportToArray(
  std::vector<EK::Point_3> &inputVertices,
  std::vector<std::pair<int, int>> &inputEdges,
  std::vector<std::tuple<int, int, int>> &inputTris,
  std::vector<int> *vertexIDs)
{
  inputVertices.clear();
  inputEdges.clear();
  inputTris.clear();

  std::map<int, int> vidMap;
  for (const auto &pr : vertices) {
    inputVertices.emplace_back(pr.second.pos);
    vidMap.emplace(pr.first, (int)inputVertices.size() - 1);
  }

  for (const auto &pr : edges) {
    auto it0 = vidMap.find(pr.second.vidx[0]);
    PGO_ALOG(it0 != vidMap.end());

    auto it1 = vidMap.find(pr.second.vidx[1]);
    PGO_ALOG(it1 != vidMap.end());
    inputEdges.emplace_back(it0->second, it1->second);
  }

  for (const auto &pr : triangles) {
    auto it0 = vidMap.find(pr.second.vidx[0]);
    PGO_ALOG(it0 != vidMap.end());

    auto it1 = vidMap.find(pr.second.vidx[1]);
    PGO_ALOG(it1 != vidMap.end());

    auto it2 = vidMap.find(pr.second.vidx[2]);
    PGO_ALOG(it2 != vidMap.end());
    inputTris.emplace_back(it0->second, it1->second, it2->second);
  }

  if (vertexIDs) {
    vertexIDs->assign(inputVertices.size(), -1);
    for (const auto &pr : vidMap) {
      vertexIDs->at(pr.second) = pr.first;
    }
  }
}

void MedialAxisRepresentation::SkeletonExact::exportToArray(
  std::vector<EK::Point_3> &inputVertices, std::vector<std::pair<int, int>> &inputEdges, std::vector<std::tuple<int, int, int>> &inputTris,
  std::map<int, int> &vidMap,                        // json[vid] -> coref id
  std::map<int, std::pair<int, int>> &edgeMap,       // coref_id -> (json[vid], json[vid])
  std::map<int, std::tuple<int, int, int>> &triMap)  // coref_id -> (json[vid], json[vid], json[vid])
{
  inputVertices.clear();
  inputEdges.clear();
  inputTris.clear();

  for (const auto &pr : vertices) {
    inputVertices.emplace_back(pr.second.pos);
    vidMap.emplace(pr.first, (int)inputVertices.size() - 1);
  }

  for (const auto &pr : edges) {
    auto it0 = vidMap.find(pr.second.vidx[0]);
    PGO_ALOG(it0 != vidMap.end());

    auto it1 = vidMap.find(pr.second.vidx[1]);
    PGO_ALOG(it1 != vidMap.end());
    inputEdges.emplace_back(it0->second, it1->second);

    edgeMap.emplace((int)inputEdges.size() - 1, std::make_pair(it0->first, it1->first));
  }

  for (const auto &pr : triangles) {
    auto it0 = vidMap.find(pr.second.vidx[0]);
    PGO_ALOG(it0 != vidMap.end());

    auto it1 = vidMap.find(pr.second.vidx[1]);
    PGO_ALOG(it1 != vidMap.end());

    auto it2 = vidMap.find(pr.second.vidx[2]);
    PGO_ALOG(it2 != vidMap.end());
    inputTris.emplace_back(it0->second, it1->second, it2->second);

    triMap.emplace((int)inputTris.size() - 1, std::make_tuple(it0->first, it1->first, it2->first));
  }
}

void MedialAxisRepresentation::SkeletonExact::saveDisplayMesh(const char *filename) const
{
  std::unordered_map<int, int> vidMap;

  // write vertex
  int inc = 0;
  pgo::Mesh::TriMeshGeo mesh;
  for (const auto &v : vertices) {
    mesh.addPos(toVec3(v.second.pos));
    vidMap.emplace(v.first, inc++);
  }

  // write tri
  for (const auto &t : triangles) {
    ES::V3i idx = t.second.vidx;
    idx[0] = vidMap[idx[0]];
    idx[1] = vidMap[idx[1]];
    idx[2] = vidMap[idx[2]];

    mesh.addTri(idx);
  }

  for (const auto &e : edges) {
    if (e.second.isTriangleEdge() == false) {
      const auto &v0 = getV(e.second.vidx[0]);
      const auto &v1 = getV(e.second.vidx[1]);
      pgo::Mesh::TriMeshGeo m = pgo::Mesh::createSingleTriangleMesh(toVec3(v0.pos), toVec3(v1.pos), toVec3(v0.pos) + pgo::asVec3d(1e-5));
      mesh.addMesh(m);
    }
  }

  mesh.save(filename);
}

bool MedialAxisRepresentation::SkeletonExact::isConnected() const
{
  std::vector<int> Q;
  Q.reserve(vertices.size() * 2);

  std::unordered_set<int> visited;
  visited.reserve(vertices.size() * 2);

  visited.emplace(vertices.begin()->first);
  Q.emplace_back(vertices.begin()->first);

  size_t start = 0;
  while (start < Q.size()) {
    int curVID = Q[start++];

    const auto &v = getV(curVID);
    for (int tid : v.neighboringTriangles) {
      const auto &tri = getT(tid);
      for (int j = 0; j < 3; j++) {
        if (tri.vidx[j] == curVID)
          continue;

        if (visited.find(tri.vidx[j]) != visited.end())
          continue;

        Q.emplace_back(tri.vidx[j]);
        visited.emplace(tri.vidx[j]);
      }
    }

    for (int eid : v.neighboringEdges) {
      const auto &edge = getE(eid);
      for (int j = 0; j < 2; j++) {
        if (edge.vidx[j] == curVID)
          continue;

        if (visited.find(edge.vidx[j]) != visited.end())
          continue;

        Q.emplace_back(edge.vidx[j]);
        visited.emplace(edge.vidx[j]);
      }
    }
  }

  return visited.size() == vertices.size();
}

bool MedialAxisRepresentation::SkeletonExact::isConsistent() const
{
  struct V3iLess
  {
    bool operator()(const ES::V3i &v0, const ES::V3i &v1) const
    {
      return std::memcmp(v0.data(), v1.data(), sizeof(int) * 3) < 0;
    }
  };

  std::set<ES::V3i, V3iLess> triSet;
  for (const auto &pr_tri : triangles) {
    for (int j = 0; j < 3; j++) {
      const auto &v = getV(pr_tri.second.vidx[j]);
      auto iter = std::find(v.neighboringTriangles.begin(), v.neighboringTriangles.end(), pr_tri.first);
      if (iter == v.neighboringTriangles.end()) {
        return false;
      }

      const auto &e = getE(pr_tri.second.eidx[j]);
      iter = std::find(e.neighboringTriangles.begin(), e.neighboringTriangles.end(), pr_tri.first);
      if (iter == e.neighboringTriangles.end()) {
        return false;
      }
    }

    ES::V3i triKey = pr_tri.second.vidx;
    std::sort(triKey.data(), triKey.data() + 3);
    bool ret = triSet.emplace(triKey).second;
    if (ret == false) {
      return false;
    }
  }

  for (const auto &pr_e : edges) {
    if (pr_e.second.isTriangleEdge()) {
      for (int j = 0; j < 2; j++) {
        const auto &v = getV(pr_e.second.vidx[j]);
        auto iter = std::find(v.neighboringEdges.begin(), v.neighboringEdges.end(), pr_e.first);
        if (iter != v.neighboringEdges.end()) {
          return false;
        }
      }
    }
    else {
      for (int j = 0; j < 2; j++) {
        const auto &v = getV(pr_e.second.vidx[j]);
        auto iter = std::find(v.neighboringEdges.begin(), v.neighboringEdges.end(), pr_e.first);
        if (iter == v.neighboringEdges.end()) {
          return false;
        }
      }
    }
  }

  if (triangles.size() != triVtxIDToTriID.size())
    return false;

  for (const auto &pr : triVtxIDToTriID) {
    auto it = triangles.find(pr.second);
    if (it == triangles.end())
      return false;

    if (pgo::Mesh::UTriKey(it->second.vidx.data()) != pr.first)
      return false;
  }

  for (const auto &pr : triangles) {
    auto it = triVtxIDToTriID.find(pgo::Mesh::UTriKey(pr.second.vidx.data()));
    if (it == triVtxIDToTriID.end())
      return false;

    if (it->second != pr.first)
      return false;
  }

  if (edges.size() != edgeVtxIDToEdgeID.size())
    return false;

  for (const auto &pr : edgeVtxIDToEdgeID) {
    auto it = edges.find(pr.second);
    if (it == edges.end())
      return false;

    if (std::pair<int, int>(it->second.vidx[0], it->second.vidx[1]) != pr.first)
      return false;
  }

  for (const auto &pr : edges) {
    auto it = edgeVtxIDToEdgeID.find(std::pair<int, int>(pr.second.vidx[0], pr.second.vidx[1]));
    if (it == edgeVtxIDToEdgeID.end())
      return false;

    if (it->second != pr.first)
      return false;
  }

  return true;
}

void MedialAxisRepresentation::SkeletonExact::eraseTriangle(int triID, int deleteMap)
{
  auto iter1 = triangles.find(triID);
  PGO_ALOG(iter1 != triangles.end());

  // for debug
  // for (int j = 0; j < 3; j++) {
  //  if (iter1->second.eidx[j] == 1752) {
  //    std::cin.get();
  //  }
  //}
  // const auto &eeee1 = getE(iter1->second.eidx[0]);
  // const auto &eeee2 = getE(iter1->second.eidx[1]);
  // const auto &eeee3 = getE(iter1->second.eidx[2]);

  // for each vertex/edge of the deleted triangle,
  // we remove the neighboring info
  for (int j = 0; j < 3; j++) {
    // delete vertex neighbor triangle
    auto iter_v = vertices.find(iter1->second.vidx[j]);
    PGO_ALOG(iter_v != vertices.end());

    auto iter_t = std::find(iter_v->second.neighboringTriangles.begin(), iter_v->second.neighboringTriangles.end(), triID);
    PGO_ALOG(iter_t != iter_v->second.neighboringTriangles.end());

    iter_v->second.neighboringTriangles.erase(iter_t);

    // delete edge neighbor triangle
    auto iter_e = edges.find(iter1->second.eidx[j]);
    PGO_ALOG(iter_e != edges.end());

    iter_t = std::find(iter_e->second.neighboringTriangles.begin(), iter_e->second.neighboringTriangles.end(), triID);
    PGO_ALOG(iter_t != iter_e->second.neighboringTriangles.end());

    iter_e->second.neighboringTriangles.erase(iter_t);

    // if after modification no neighboring triangles are left, we convert this edge to pure edge
    if (iter_e->second.neighboringTriangles.size() == 0) {
      auto &vv0 = getV(iter_e->second.vidx[0]);
      vv0.neighboringEdges.emplace_back(iter1->second.eidx[j]);
      pgo::BasicAlgorithms::sortAndDeduplicateWithErase(vv0.neighboringEdges);

      auto &vv1 = getV(iter_e->second.vidx[1]);
      vv1.neighboringEdges.emplace_back(iter1->second.eidx[j]);
      pgo::BasicAlgorithms::sortAndDeduplicateWithErase(vv1.neighboringEdges);
    }
  }

  if (deleteMap) {
    size_t count = triVtxIDToTriID.erase(pgo::Mesh::UTriKey(iter1->second.vidx.data()));
    PGO_ALOG(count == 1ull);
  }

  // finally delete the triangle itself
  triangles.erase(iter1);
}

void MedialAxisRepresentation::SkeletonExact::replaceVertex(int oldVID, const Vertex &oldVertex, int newVID)
{
  // change triangle vertex ID
  for (int triID : oldVertex.neighboringTriangles) {
    auto &tri = getT(triID);
    ES::V3i oldTriVtxID = tri.vidx;

    // change the vertex IDs
    bool found = false;
    for (int j = 0; j < 3; j++) {
      if (tri.vidx[j] == oldVID) {
        tri.vidx[j] = newVID;
        found = true;
        break;
      }
    }
    PGO_ALOG(found == true);

    // change the vertex IDs in the triangle edge
    for (int j = 0; j < 3; j++) {
      auto &edge1 = getE(tri.eidx[j]);

      // if (tri.eidx[j] == 1907) {
      //   std::cout << "zz\n";
      // }

      // change the vertex IDs
      bool found = false;
      std::pair<int, int> oldEdgeVID;
      for (int k = 0; k < 2; k++) {
        if (edge1.vidx[k] == oldVID) {
          oldEdgeVID = std::make_pair(edge1.vidx[0], edge1.vidx[1]);
          edge1.vidx[k] = newVID;
          found = true;

          break;
        }
      }

      if (found) {
        if (edge1.vidx[0] > edge1.vidx[1]) {
          std::swap(edge1.vidx[0], edge1.vidx[1]);
        }

        // update mapping
        std::pair<int, int> newEdgeVID = std::make_pair(edge1.vidx[0], edge1.vidx[1]);
        auto iter_id = edgeVtxIDToEdgeID.find(oldEdgeVID);
        PGO_ALOG(iter_id != edgeVtxIDToEdgeID.end());
        edgeVtxIDToEdgeID.erase(iter_id);

        // check if there are any duplications
        iter_id = edgeVtxIDToEdgeID.find(newEdgeVID);
        // if yes
        if (iter_id != edgeVtxIDToEdgeID.end()) {
          int thisEdgeID = tri.eidx[j];
          int existingEdgeID = iter_id->second;
          // we remove it
          auto thisEdge = edges.find(thisEdgeID);
          auto existingEdge = edges.find(existingEdgeID);

          // for each triangle that is attached to it,
          // we relocate it to the existing edge
          for (int triID1 : thisEdge->second.neighboringTriangles) {
            // build connectivity to the existing one
            existingEdge->second.neighboringTriangles.emplace_back(triID1);

            // replace this edge from triangle neighbor
            auto &tri1 = getT(triID1);
            for (int k = 0; k < 3; k++) {
              if (tri1.eidx[k] == thisEdgeID) {
                tri1.eidx[k] = existingEdgeID;
              }
            }
          }
          pgo::BasicAlgorithms::sortAndDeduplicateWithErase(existingEdge->second.neighboringTriangles);

          // if it is no longer a pure edge,
          // we modify the neighboring vertex
          if (existingEdge->second.neighboringTriangles.size()) {
            for (int j = 0; j < 2; j++) {
              auto &v = getV(existingEdge->second.vidx[j]);
              auto it = std::find(v.neighboringEdges.begin(), v.neighboringEdges.end(), existingEdgeID);
              if (it != v.neighboringEdges.end()) {
                v.neighboringEdges.erase(it);
              }
            }
          }

          edges.erase(thisEdge);
        }
        // if no
        else {
          edgeVtxIDToEdgeID.emplace(newEdgeVID, tri.eidx[j]);
        }
      }
    }

    // address duplication of triangles
    auto it_old = triVtxIDToTriID.find(pgo::Mesh::UTriKey(oldTriVtxID.data()));
    PGO_ALOG(it_old != triVtxIDToTriID.end());

    auto it_new = triVtxIDToTriID.find(pgo::Mesh::UTriKey(tri.vidx.data()));

    // if this new triangle vtx ID has been accessed
    if (it_new != triVtxIDToTriID.end()) {
      // we erase this triangle
      eraseTriangle(triID, 0);
      triVtxIDToTriID.erase(it_old);
    }
    // otherwise
    else {
      triVtxIDToTriID.erase(it_old);
      triVtxIDToTriID.emplace(pgo::Mesh::UTriKey(tri.vidx.data()), triID);
    }
  }

  // change edge vertex ID
  for (int eid : oldVertex.neighboringEdges) {
    auto &edge1 = getE(eid);
    bool found = false;
    std::pair<int, int> oldEdgeVID;
    for (int j = 0; j < 2; j++) {
      if (edge1.vidx[j] == oldVID) {
        oldEdgeVID = std::make_pair(edge1.vidx[0], edge1.vidx[1]);
        edge1.vidx[j] = newVID;
        found = true;
        break;
      }
    }

    PGO_ALOG(found != false);

    if (edge1.vidx[0] > edge1.vidx[1]) {
      std::swap(edge1.vidx[0], edge1.vidx[1]);
    }

    // update mapping
    std::pair<int, int> newEdgeVID = std::make_pair(edge1.vidx[0], edge1.vidx[1]);
    auto iter_id = edgeVtxIDToEdgeID.find(oldEdgeVID);
    PGO_ALOG(iter_id != edgeVtxIDToEdgeID.end());
    edgeVtxIDToEdgeID.erase(iter_id);

    // check if there are any duplications
    iter_id = edgeVtxIDToEdgeID.find(newEdgeVID);
    // if yes
    // we remove it
    if (iter_id != edgeVtxIDToEdgeID.end()) {
      int thisEdgeID = eid;
      auto thisEdge = edges.find(thisEdgeID);

      // for each vertex, we remove this edge
      for (int j = 0; j < 2; j++) {
        if (thisEdge->second.vidx[j] == oldVID)
          continue;

        auto &v = getV(thisEdge->second.vidx[j]);
        auto it = std::find(v.neighboringEdges.begin(), v.neighboringEdges.end(), thisEdgeID);
        PGO_ALOG(it != v.neighboringEdges.end());
        v.neighboringEdges.erase(it);
      }

      // erase this edge
      edges.erase(thisEdge);
    }
    // if no
    else {
      edgeVtxIDToEdgeID.emplace(newEdgeVID, eid);
    }
  }
}