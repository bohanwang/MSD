#include "skeletonRVD.h"

#include "pgoLogging.h"
#include "EigenSupport.h"
#include "createTriMesh.h"
#include "basicAlgorithms.h"
#include "disjointSet.h"

#include <CGAL/Polygon_mesh_processing/orient_polygon_soup.h>
#include <CGAL/Polygon_mesh_processing/triangulate_faces.h>
#include <CGAL/Polygon_mesh_processing/corefinement.h>
#include <CGAL/Polygon_mesh_processing/remesh.h>
#include <CGAL/AABB_triangle_primitive.h>
#include <CGAL/Barycentric_coordinates_2/triangle_coordinates_2.h>
#include <CGAL/Polygon_mesh_processing/orient_polygon_soup_extension.h>
#include <CGAL/Cartesian_converter.h>
#include <CGAL/Side_of_triangle_mesh.h>
#include <CGAL/Exact_integer.h>
#include <CGAL/Homogeneous.h>
#include <CGAL/Extended_homogeneous.h>
#include <CGAL/cartesian_homogeneous_conversion.h>
#include <CGAL/IO/Nef_polyhedron_iostream_3.h>
#include <CGAL/Polygon_mesh_processing/triangulate_faces.h>

#include <geogram/basic/stopwatch.h>
#include <geogram/basic/file_system.h>
#include <geogram/mesh/mesh.h>
#include <geogram/mesh/mesh_io.h>
#include <geogram/mesh/mesh_reorder.h>
#include <geogram/delaunay/delaunay.h>
#include <geogram/delaunay/periodic_delaunay_3d.h>
#include <geogram/delaunay/delaunay_3d.h>
#include <geogram/voronoi/convex_cell.h>

#include <tbb/parallel_for.h>
#include <tbb/concurrent_vector.h>
#include <nlohmann/json.hpp>

#include <fmt/format.h>

#include <unordered_map>

namespace MedialAxisRepresentation
{
using HClock = std::chrono::high_resolution_clock;
using HClockPt = HClock::time_point;

namespace ES = pgo::EigenSupport;

inline double duraSecond(const HClockPt &t0, const HClockPt &t1)
{
  return double(std::chrono::duration_cast<std::chrono::microseconds>(t1 - t0).count()) * 1e-6;
}

inline ES::V3d toVec3(const EK::Point_3 &pt)
{
  return ES::V3d(CGAL::to_double(pt.x()), CGAL::to_double(pt.y()), CGAL::to_double(pt.z()));
}

inline ES::V3d toVec3(const EK::Vector_3 &pt)
{
  return ES::V3d(CGAL::to_double(pt[0]), CGAL::to_double(pt[1]), CGAL::to_double(pt[2]));
}

template<typename Kernel, typename KernelExt>
class PolyhedronConverter : public CGAL::Modifier_base<typename pgo::CGALInterface::Polyhedron<Kernel>::HalfedgeDS>
{
protected:
  const pgo::CGALInterface::Polyhedron<KernelExt> &polyIn;

public:
  PolyhedronConverter(const typename pgo::CGALInterface::Polyhedron<KernelExt> &polyIn_):
    polyIn(polyIn_)
  {
  }

  using HalfedgeDS = typename pgo::CGALInterface::Polyhedron<Kernel>::HalfedgeDS;

  void operator()(HalfedgeDS &hds)
  {
    typedef typename HalfedgeDS::Vertex Vertex;

    int numVertices = (int)polyIn.size_of_vertices();
    int numFaces = (int)polyIn.size_of_facets();

    // create a cgal incremental builder
    CGAL::Polyhedron_incremental_builder_3<HalfedgeDS> B(hds, true);
    B.begin_surface(numVertices, numFaces);

    std::map<typename pgo::CGALInterface::Polyhedron<KernelExt>::Vertex_const_handle, int> vtxMap;

    // add the polyhedron vertices
    int inc = 0;
    for (auto it = polyIn.vertices_begin(); it != polyIn.vertices_end(); ++it) {
      const auto &pt_ext = it->point();
      PGO_ALOG(pt_ext[0].degree() == 0);
      PGO_ALOG(pt_ext[1].degree() == 0);
      PGO_ALOG(pt_ext[2].degree() == 0);

      typename Kernel::Point_3 pt(pt_ext[0][0], pt_ext[1][0], pt_ext[2][0]);

      typename HalfedgeDS::Vertex_handle vit = B.add_vertex(pt);
      vtxMap.emplace(it, inc);

      vit->id() = inc++;
    }

    inc = 0;
    for (auto it = polyIn.facets_begin(); it != polyIn.facets_end(); ++it) {
      typename HalfedgeDS::Face_handle fit = B.begin_facet();
      typename pgo::CGALInterface::Polyhedron<KernelExt>::Halfedge_const_handle h0 = it->halfedge();
      typename pgo::CGALInterface::Polyhedron<KernelExt>::Halfedge_const_handle h1 = h0;

      // std::cout << "face " << inc << '\n';
      do {
        auto vit = vtxMap.find(h1->vertex());
        PGO_ALOG(vit != vtxMap.end());

        // std::cout << vit->second << ' ';
        B.add_vertex_to_facet(vit->second);

        h1 = h1->next();
      } while (h1 != h0);
      // std::cout << std::endl;

      B.end_facet();

      fit->id() = inc++;
    }

    // finish up the surface
    B.end_surface();
  }
};

void triangulatePolygon(const std::vector<EK::Point_3> &polylineVtxs,
  const pgo::Mesh::TriMeshGeo &targetMesh,
  const pgo::Mesh::TriMeshBVTree &targetMeshBVTree,
  std::vector<CGAL::Triple<int, int, int>> &triangles);
}  // namespace MedialAxisRepresentation

void MedialAxisRepresentation::computeRVD_internal_Dijkstras_withTri_refine_withEdge(const pgo::Mesh::TriMeshGeo &targetMesh, const pgo::Mesh::TriMeshBVTree &targetMeshBVTree,
  const pgo::Mesh::TriMeshGeo &targetMeshSmall, const pgo::Mesh::TriMeshBVTree &targetMeshBVTreeSmall,
  const std::vector<ES::V3d> &samplePoints,
  std::vector<EK::Point_3> &medialAxisVertices, std::vector<std::vector<int>> &medialAxisFacets,
  std::vector<ES::V3d> &finalSkeletonPoints, std::vector<std::pair<int, int>> &finalSkeletonEdges, std::vector<std::tuple<int, int, int>> &finalSkeletonTriangles,
  std::vector<ES::V3d> *cellCenters, const std::string &saveExactResults)
{
  constexpr bool dumpMesh = false;

  // create bounding box
  pgo::Mesh::BoundingBox bb(targetMesh.positions());
  bb.expand(1.1);

  HClockPt t1 = HClock::now();

  pgo::Mesh::TriMeshGeo pp;
  for (int si = 0; si < (int)samplePoints.size(); si++) {
    pp.addPos(samplePoints[si]);
  }
  pp.save("latest-sample.obj");

  // ==============================================
  // perform delaunay for the target mesh
  std::vector<double> points(samplePoints.size() * 3);
  for (size_t i = 0; i < samplePoints.size(); i++) {
    points[i * 3] = samplePoints[i][0];
    points[i * 3 + 1] = samplePoints[i][1];
    points[i * 3 + 2] = samplePoints[i][2];
  }

#if 0
  GEO::SmartPointer<GEO::PeriodicDelaunay3d> delaunay = new GEO::PeriodicDelaunay3d(false, 1.0);
  delaunay->set_keeps_infinite(true);
  delaunay->set_vertices((GEO::index_t)samplePoints.size(), points.data());
  delaunay->compute();
#else
  GEO::SmartPointer<GEO::Delaunay3d> delaunay = new GEO::Delaunay3d;
  delaunay->set_keeps_infinite(true);
  delaunay->set_vertices((GEO::index_t)samplePoints.size(), points.data());
#endif

  auto reorderIfNot = [](std::array<int, 2> &edgeID) {
    if (edgeID[0] > edgeID[1]) {
      std::swap(edgeID[0], edgeID[1]);
    }
  };

  // ==============================================
  // construct voronoi cells
  std::unordered_map<std::array<int, 2>, std::vector<int>, ES::IntArrayHash<2>, ES::IntArrayEqual<2>> edgeID_to_tetIDs;
  std::unordered_map<std::array<int, 3>, std::vector<int>, ES::IntArrayHash<3>, ES::IntArrayEqual<3>> faceID_to_tetIDs;

  std::unordered_map<int, std::vector<int>> vertexID_to_tetIDs;
  std::vector<EK::Point_3> vertices(delaunay->nb_vertices());
  std::unordered_map<int, std::array<int, 4>> tets;
  std::unordered_map<int, EK::Point_3> tetCenters;

  if constexpr (1) {
    tets.reserve(delaunay->nb_cells());

    // add vertices
    for (GEO::index_t i = 0; i < delaunay->nb_vertices(); ++i) {
      GEO::vec3 pos(delaunay->vertex_ptr(i)[0], delaunay->vertex_ptr(i)[1], delaunay->vertex_ptr(i)[2]);
      vertices[i] = EK::Point_3(pos[0], pos[1], pos[2]);
    }

    // build tets, edge -> tets, vertex -> tets, face -> tets
    for (GEO::index_t i = 0; i < delaunay->nb_cells(); i++) {
      std::array<int, 4> tetVertexIndices;
      for (GEO::index_t j = 0; j < 4; j++) {
        tetVertexIndices[j] = (int)delaunay->cell_vertex(i, j);
      }
      auto tetIt = tets.emplace((int)i, tetVertexIndices);

      for (int vi = 0; vi < 4; vi++) {
        for (int vj = vi + 1; vj < 4; vj++) {
          std::array<int, 2> edgeID{ tetIt.first->second[vi], tetIt.first->second[vj] };
          reorderIfNot(edgeID);

          auto it = edgeID_to_tetIDs.find(edgeID);
          if (it != edgeID_to_tetIDs.end()) {
            it->second.emplace_back((int)i);
          }
          else {
            auto ret = edgeID_to_tetIDs.emplace_hint(it, edgeID, std::vector<int>());
            ret->second.emplace_back((int)i);
          }
        }
      }

      for (int vi = 0; vi < 4; vi++) {
        auto it = vertexID_to_tetIDs.find(tetIt.first->second[vi]);
        if (it != vertexID_to_tetIDs.end()) {
          it->second.emplace_back((int)i);
        }
        else {
          auto ret = vertexID_to_tetIDs.emplace_hint(it, tetIt.first->second[vi], std::vector<int>());
          ret->second.emplace_back((int)i);
        }
      }

      for (int vi = 0; vi < 4; vi++) {
        for (int vj = vi + 1; vj < 4; vj++) {
          for (int vk = vj + 1; vk < 4; vk++) {
            std::array<int, 3> faceID{ tetIt.first->second[vi], tetIt.first->second[vj], tetIt.first->second[vk] };
            std::sort(faceID.begin(), faceID.end());

            auto it = faceID_to_tetIDs.find(faceID);
            if (it != faceID_to_tetIDs.end()) {
              it->second.emplace_back((int)i);
            }
            else {
              auto ret = faceID_to_tetIDs.emplace_hint(it, faceID, std::vector<int>());
              ret->second.emplace_back((int)i);
            }
          }
        }
      }
    }

    // tet centers are the voronoi cell vertices
    for (const auto &pr : tets) {
      if (delaunay->cell_is_infinite(pr.first)) {
        // tetCenters.emplace(pr.first, EK::Point_3(100, 100, 100));
        // approximate infinite cell
        std::vector<int> nonInfID;
        nonInfID.reserve(4);
        for (int i = 0; i < 4; i++) {
          if (pr.second[i] < 0) {
            continue;
          }
          else {
            nonInfID.emplace_back(pr.second[i]);
          }
        }
        PGO_ALOG((int)nonInfID.size() <= 3);

        if ((int)nonInfID.size() == 1) {
          tetCenters.emplace(pr.first, vertices[nonInfID[0]]);
        }
        else if ((int)nonInfID.size() == 2) {
          EK::Point_3 center = CGAL::midpoint(vertices[nonInfID[0]], vertices[nonInfID[1]]);
          tetCenters.emplace(pr.first, center);
        }
        else {
          PGO_ALOG((int)nonInfID.size() == 3);
          EK::Point_3 center = CGAL::circumcenter(vertices[nonInfID[0]], vertices[nonInfID[1]], vertices[nonInfID[2]]);
          tetCenters.emplace(pr.first, center);
        }
      }
      else {
        EK::Point_3 center = CGAL::circumcenter(vertices[pr.second[0]], vertices[pr.second[1]], vertices[pr.second[2]], vertices[pr.second[3]]);
        tetCenters.emplace(pr.first, center);
      }
    }

    // save
    if constexpr (0) {
      nlohmann::json jout;
      // edgeID_to_tetIDs
      for (auto &it : edgeID_to_tetIDs) {
        std::string key = fmt::format("{}_{}", it.first[0], it.first[1]);
        jout["edgeID_to_tetIDs"][key] = it.second;
      }
      // faceID_to_tetIDs
      for (auto &it : faceID_to_tetIDs) {
        std::string key = fmt::format("{}_{}_{}", it.first[0], it.first[1], it.first[2]);
        jout["faceID_to_tetIDs"][key] = it.second;
      }
      // vertexID_to_tetIDs
      for (auto &it : vertexID_to_tetIDs) {
        std::string key = fmt::format("{}", it.first);
        jout["vertexID_to_tetIDs"][key] = it.second;
      }
      // vertices
      for (int i = 0; i < vertices.size(); i++) {
        std::string key = fmt::format("{}", i);
        std::vector<char> val[3];
        for (int j = 0; j < 3; j++) {
          auto &op = vertices[i][j].mpq();
          val[j].resize((mpz_sizeinbase(mpq_numref(op), 10) + mpz_sizeinbase(mpq_denref(op), 10) + 3) * 2);
          memset(val[j].data(), 0, sizeof(char) * val[j].size());
          mpq_get_str(val[j].data(), 10, op);
        }
        jout["vertices"][key] = { val[0], val[1], val[2] };
      }
      // tets
      for (auto &it : tets) {
        std::string key = fmt::format("{}", it.first);
        jout["tets"][key] = it.second;
      }
      // tetCenters
      for (auto &it : tetCenters) {
        std::string key = fmt::format("{}", it.first);
        std::vector<char> val[3];
        for (int j = 0; j < 3; j++) {
          auto &op = it.second[j].mpq();
          val[j].resize((mpz_sizeinbase(mpq_numref(op), 10) + mpz_sizeinbase(mpq_denref(op), 10) + 3) * 2);
          memset(val[j].data(), 0, sizeof(char) * val[j].size());
          mpq_get_str(val[j].data(), 10, op);
        }
        jout["tetCenters"][key] = { val[0], val[1], val[2] };
      }

      // std::ofstream(jsonIn0.c_str()) << jout.dump(2);
      // fmt::print("saved to {}\n", jsonIn0);
    }
  }

#if 0
  if constexpr (0) {
    std::ifstream ifs(jsonIn0.c_str());
    nlohmann::json jin;
    ifs >> jin;
    for (auto it = jin["edgeID_to_tetIDs"].begin(); it != jin["edgeID_to_tetIDs"].end(); it++) {
      std::array<int, 2> edgeID;
      sscanf(it.key().c_str(), "%d_%d", &edgeID[0], &edgeID[1]);
      edgeID_to_tetIDs[edgeID] = it.value().get<std::vector<int>>();
    }
    for (auto it = jin["faceID_to_tetIDs"].begin(); it != jin["faceID_to_tetIDs"].end(); it++) {
      std::array<int, 3> faceID;
      sscanf(it.key().c_str(), "%d_%d_%d", &faceID[0], &faceID[1], &faceID[2]);
      faceID_to_tetIDs[faceID] = it.value().get<std::vector<int>>();
    }
    for (auto it = jin["vertexID_to_tetIDs"].begin(); it != jin["vertexID_to_tetIDs"].end(); it++) {
      int vertexID;
      sscanf(it.key().c_str(), "%d", &vertexID);
      vertexID_to_tetIDs[vertexID] = it.value().get<std::vector<int>>();
    }
    vertices.resize(jin["vertices"].size());
    for (auto it = jin["vertices"].begin(); it != jin["vertices"].end(); it++) {
      int idx = std::stoi(it.key());
      std::vector<char> valStr[3];
      for (int i = 0; i < 3; i++) {
        valStr[i] = it.value()[i];
      }
      vertices[idx] = EK::Point_3(CGAL::Gmpq(valStr[0].data(), 10),
        CGAL::Gmpq(valStr[1].data(), 10),
        CGAL::Gmpq(valStr[2].data(), 10));
    }
    for (auto it = jin["tets"].begin(); it != jin["tets"].end(); it++) {
      int idx = std::stoi(it.key());
      tets[idx] = it.value().get<std::array<int, 4>>();
    }
    for (auto it = jin["tetCenters"].begin(); it != jin["tetCenters"].end(); it++) {
      int idx = std::stoi(it.key());
      std::vector<char> valStr[3];
      for (int i = 0; i < 3; i++) {
        valStr[i] = it.value()[i];
      }
      tetCenters[idx] = EK::Point_3(CGAL::Gmpq(valStr[0].data(), 10),
        CGAL::Gmpq(valStr[1].data(), 10),
        CGAL::Gmpq(valStr[2].data(), 10));
    }
    fmt::print("load from {}\n", jsonIn0);
  }
#endif

  HClockPt t2 = HClock::now();

  SPDLOG_LOGGER_INFO(pgo::Logging::lgr(), "voronoi cell: {}s", duraSecond(t1, t2));

  // this block is to dump delaunay triangulation
  if constexpr (dumpMesh) {
    std::ofstream outfile("deb.obj");
    for (int vi = 0; vi < (int)vertices.size(); vi++) {
      ES::V3d p = toVec3(vertices[vi]);
      outfile << "v " << p[0] << ' ' << p[1] << ' ' << p[2] << '\n';
    }

    for (const auto &pr : tets) {
      if (delaunay->cell_is_infinite(pr.first))
        continue;

      outfile << "f " << (pr.second[0] + 1) << ' ' << (pr.second[1] + 1) << ' ' << (pr.second[2] + 1) << '\n';
      outfile << "f " << (pr.second[0] + 1) << ' ' << (pr.second[2] + 1) << ' ' << (pr.second[3] + 1) << '\n';
      outfile << "f " << (pr.second[0] + 1) << ' ' << (pr.second[1] + 1) << ' ' << (pr.second[3] + 1) << '\n';
      outfile << "f " << (pr.second[1] + 1) << ' ' << (pr.second[2] + 1) << ' ' << (pr.second[3] + 1) << '\n';
    }
    outfile.close();
  }

  auto tetsAreFaceNeighbor = [&](int ti, int tj) -> bool {
    int foundCounter = 0;
    const auto &tetVIDs_i = tets[ti];
    const auto &tetVIDs_j = tets[tj];

    for (int i = 0; i < 4; i++) {
      for (int j = 0; j < 4; j++) {
        if (tetVIDs_i[i] == tetVIDs_j[j]) {
          foundCounter++;
          break;
        }
      }
    }

    if (foundCounter == 3) {
      return true;
    }
    else {
      return false;
    }
  };

  // reorder vertices of each voronoi polygon faces
  for (auto &e_to_t : edgeID_to_tetIDs) {
    for (int ti = 0; ti < (int)e_to_t.second.size() - 1; ti++) {
      for (int tj = ti + 1; tj < (int)e_to_t.second.size(); tj++) {
        if (tetsAreFaceNeighbor(e_to_t.second[ti], e_to_t.second[tj])) {
          std::swap(e_to_t.second[ti + 1], e_to_t.second[tj]);
          break;
        }
      }
    }
  }

  // find all the edges of Voronoi cells
  // edgeID indicates the pos of tetCenter
  std::unordered_map<std::array<int, 3>, std::array<int, 2>, ES::IntArrayHash<3>, ES::IntArrayEqual<3>> delaunayTri2VoronoiCellEdges;
  for (const auto &f_to_t : faceID_to_tetIDs) {
    // ignore infinite faces
    if (f_to_t.first[0] < 0 || f_to_t.first[1] < 0 || f_to_t.first[2] < 0)
      continue;
    PGO_ALOG(f_to_t.second.size() <= 2);
    // skip infinite face
    if (f_to_t.second.size() == 2) {
      std::array<int, 2> edgeID{ f_to_t.second[0], f_to_t.second[1] };
      reorderIfNot(edgeID);
      auto it = delaunayTri2VoronoiCellEdges.find(f_to_t.first);
      PGO_ALOG(it == delaunayTri2VoronoiCellEdges.end());
      auto ret = delaunayTri2VoronoiCellEdges.emplace_hint(it, f_to_t.first, edgeID);
    }
  }
  fmt::print("delaunayTri2VoronoiCellEdges.size() = {}\n", delaunayTri2VoronoiCellEdges.size());

  using EK_Ext = CGAL::Extended_cartesian<CGAL::Gmpq>;
  using Nef = CGAL::Nef_polyhedron_3<EK_Ext>;
  using Poly_Ext = pgo::CGALInterface::Polyhedron<EK_Ext>;

  using Poly = pgo::CGALInterface::Polyhedron<EK>;
  using FaceGraphPrimitive = CGAL::AABB_face_graph_triangle_primitive<Poly>;
  using AABB_face_graph_traits = CGAL::AABB_traits<EK, FaceGraphPrimitive>;
  using FaceGraphTree = CGAL::AABB_tree<AABB_face_graph_traits>;

  HClockPt t3 = HClock::now();

  // build finite bounding box
  Nef boundingBox;
  if constexpr (1) {
    ES::V3d bmin = bb.bmin();
    ES::V3d bmax = bb.bmax();

    Nef N1(EK_Ext::Plane_3(1, 0, 0, -bmax[0]));
    Nef N2(EK_Ext::Plane_3(-1, 0, 0, bmin[0]));
    Nef N3(EK_Ext::Plane_3(0, 1, 0, -bmax[1]));
    Nef N4(EK_Ext::Plane_3(0, -1, 0, bmin[1]));
    Nef N5(EK_Ext::Plane_3(0, 0, 1, -bmax[2]));
    Nef N6(EK_Ext::Plane_3(0, 0, -1, bmin[2]));

    boundingBox = N1 * N2 * N3 * N4 * N5 * N6;
    std::cout << "bb:" << boundingBox.is_space() << ',' << boundingBox.is_simple() << std::endl;

    if (dumpMesh) {
      Poly_Ext bbPoly;
      boundingBox.convert_to_polyhedron(bbPoly);
      std::ofstream("bb.off") << bbPoly;
    }
  }

  std::vector<std::vector<int>> medialAxisSkeletonTrianglesAndEdges;
  medialAxisSkeletonTrianglesAndEdges.reserve(medialAxisFacets.size() * 20);
  enum FaceStatus : int
  {
    FS_UNADDRESSED = -1,
    FS_DELETED = -10,
    FS_INTERSECTED = -3,
  };

  std::vector<int> faceFlags;  // (medialAxisSkeletonTrianglesAndEdges.size(), FS_UNADDRESSED);

  HClockPt t5 = HClock::now();

  // std::string jsonIn1 = fmt::format("{}/jsonIn1.json", savePath);
  // std::ifstream f1(jsonIn1.c_str());
  if constexpr (1) {
    // build finite voronoi cells,
    // if there is a infinite cell, we intersect it with the bounding volume
    std::vector<Nef> voronoiCells(delaunay->nb_vertices());
    std::vector<Poly> voronoiCellPolys(delaunay->nb_vertices());

    // for (GEO::index_t vi = 0; vi < delaunay->nb_vertices(); vi++) {
    tbb::parallel_for(0, (int)delaunay->nb_vertices(), [&](int vi) {
      const EK::Point_3 &pt = vertices[(GEO::index_t)vi];
      EK_Ext::Point_3 pt_e(pt.x(), pt.y(), pt.z());

      Nef convexCell;
      bool firstTime = true;

      for (const auto &e_to_t : edgeID_to_tetIDs) {
        // ignore infinite edges
        if (e_to_t.first[0] < 0 || e_to_t.first[1] < 0)
          continue;

        // find the edge that has vi
        // if (e_to_t.first[0] == (int)vi || e_to_t.first[1] == (int)vi) {
        if (e_to_t.first[0] == vi || e_to_t.first[1] == vi) {
          // voronoi cell face is perpendicular to n all the time
          EK::Vector_3 n = vertices[e_to_t.first[0]] - vertices[e_to_t.first[1]];
          EK_Ext::Vector_3 n_e(n.x(), n.y(), n.z());

          int selTetID = -1;
          for (int ti : e_to_t.second) {
            if (delaunay->cell_is_finite(ti)) {
              selTetID = ti;
              break;
            }
          }
          PGO_ALOG(selTetID >= 0);

          // in addition, the plane goes through any tet centers // ??? why
          // EK::Point_3 p0 = tetCenters[selTetID];
          // the plane goes through the edge center
          EK::Point_3 p0 = CGAL::midpoint(vertices[e_to_t.first[0]], vertices[e_to_t.first[1]]);
          EK_Ext::Point_3 p0_e(p0.x(), p0.y(), p0.z());
          EK_Ext::Vector_3 diff_e = pt_e - p0_e;
          // nT(x - x0) = 0
          // nT x - nT x0
          if (CGAL::scalar_product(diff_e, n_e) > 0) {
            n_e *= -1;
          }

          EK_Ext::Plane_3 plane(n_e[0], n_e[1], n_e[2], -CGAL::scalar_product(n_e, EK_Ext::Vector_3(p0_e[0], p0_e[1], p0_e[2])));
          Nef halfSpace(plane);
          if (firstTime) {
            convexCell = halfSpace;
            firstTime = false;
          }
          else {
            convexCell = convexCell.intersection(halfSpace);
          }
          // finally, we intersect our cell to the bounding box to make sure it is finite and bounded.
          convexCell = convexCell * boundingBox;
        }
      }

      convexCell = convexCell * boundingBox;
      // voronoiCells.emplace_back(convexCell);
      voronoiCells[vi] = convexCell;

      Poly_Ext tempPoly;
      convexCell.convert_to_polyhedron(tempPoly);

      // construct the ek version of the cell,
      // this will be used for later cutting
      PolyhedronConverter<EK, EK_Ext> cvter(tempPoly);
      voronoiCellPolys[vi].delegate(cvter);

      if (dumpMesh) {
        std::ofstream(fmt::format("cc{}.off", vi)) << voronoiCellPolys[vi];
        pgo::Mesh::TriMeshGeo m;
        Poly triPoly = voronoiCellPolys[vi];
        CGAL::Polygon_mesh_processing::triangulate_faces(CGAL::faces(triPoly), triPoly);

        pgo::CGALInterface::polyhedron2TriangleMesh(triPoly, m);
        m.save(fmt::format("cc{}.obj", vi));
      }
      // }
    });

    HClockPt t3Half = HClock::now();
    SPDLOG_LOGGER_INFO(pgo::Logging::lgr(), "voronoi cell poly: {}s", duraSecond(t3, t3Half));

    std::vector<std::shared_ptr<FaceGraphTree>> voronoiCellBVTree(voronoiCellPolys.size());
    std::vector<std::shared_ptr<CGAL::Side_of_triangle_mesh<Poly, EK>>> voronoiSideChecker(voronoiCellPolys.size());

    tbb::parallel_for(0, (int)voronoiCellPolys.size(), [&](int vi) {
      voronoiCellBVTree[vi] = std::make_shared<FaceGraphTree>(CGAL::faces(voronoiCellPolys[vi]).first, CGAL::faces(voronoiCellPolys[vi]).second, voronoiCellPolys[vi]);
      voronoiSideChecker[vi] = std::make_shared<CGAL::Side_of_triangle_mesh<Poly, EK>>(voronoiCellPolys[vi]);
    });
    // for (const auto &p : voronoiCellPolys) {
    //   voronoiCellBVTree.emplace_back(CGAL::faces(p).first, CGAL::faces(p).second, p);
    //   voronoiSideChecker.emplace_back(p);
    // }

    HClockPt t4 = HClock::now();
    SPDLOG_LOGGER_INFO(pgo::Logging::lgr(), "voronoi cell mesh: {}s", duraSecond(t3Half, t4));

    std::vector<EK::Point_3> polylines;
    std::vector<CGAL::Triple<int, int, int>> triangles;
    for (const auto &polygon : medialAxisFacets) {
      polylines.clear();

      if (polygon.size() <= 3ull && polygon.size() > 1ull) {
        medialAxisSkeletonTrianglesAndEdges.emplace_back(polygon);
        continue;
      }
      else if (polygon.size() <= 1ull) {
        std::cout << "Encounter single vertex skeleton" << std::endl;
        continue;
      }

      for (int i : polygon)
        polylines.emplace_back(medialAxisVertices[i]);

      triangles.clear();
      CGAL::Polygon_mesh_processing::triangulate_hole_polyline(polylines, std::back_inserter(triangles));
      PGO_ALOG(triangles.size() > 0ull);

      for (const auto &tri : triangles) {
        std::vector<int> newFaces{ polygon[tri.first], polygon[tri.second], polygon[tri.third] };
        medialAxisSkeletonTrianglesAndEdges.emplace_back(std::move(newFaces));
      }
    }

    std::vector<int> faceParents(medialAxisSkeletonTrianglesAndEdges.size(), -1);
    std::unordered_map<int, std::vector<std::tuple<EK::Point_3, int>>> addedPointGlobalIndices;
    int originalEnd = (int)medialAxisVertices.size();

    auto breakSeg = [&](int ci, int fi, bool keep) {
      std::vector<FaceGraphTree::Intersection_and_primitive_id<EK::Triangle_3>::Type> intersections;
      // const EK::Point_3 &p0 = medialAxisVertices[medialAxisSkeletonTrianglesAndEdges[fi][0]];
      // const EK::Point_3 &p1 = medialAxisVertices[medialAxisSkeletonTrianglesAndEdges[fi][1]];
      const EK::Point_3 p0 = medialAxisVertices[medialAxisSkeletonTrianglesAndEdges[fi][0]];
      const EK::Point_3 p1 = medialAxisVertices[medialAxisSkeletonTrianglesAndEdges[fi][1]];

      // check if there are any intersections
      EK::Segment_3 seg(p0, p1);
      voronoiCellBVTree[ci]->all_intersections(seg, std::back_inserter(intersections));

      std::vector<std::tuple<EK::FT, EK::Point_3>> allPoints;
      allPoints.reserve(intersections.size() * 2 + 2);
      allPoints.emplace_back(0, p0);
      allPoints.emplace_back(CGAL::squared_distance(p0, p1), p1);

      // check intersection case
      for (const auto &site : intersections) {
        if (const EK::Segment_3 *s = boost::get<EK::Segment_3>(&site.first)) {
          for (int ei = 0; ei < 2; ei++) {
            bool found = false;
            for (const auto &pr : allPoints) {
              if (std::get<1>(pr) == s->vertex(ei)) {
                found = true;
                break;
              }
            }
            if (!found) {
              allPoints.emplace_back(CGAL::squared_distance(p0, s->vertex(ei)), s->vertex(ei));
            }
          }
        }
        // if the intersection is one point, it means the point is on the boundary of the cell but the others in the cell
        else if (const EK::Point_3 *p = boost::get<EK::Point_3>(&site.first)) {
          bool found = false;
          for (const auto &pr : allPoints) {
            if (std::get<1>(pr) == *p) {
              found = true;
              break;
            }
          }

          if (!found) {
            allPoints.emplace_back(CGAL::squared_distance(p0, *p), *p);
          }
        }
        // if the intersection is more than one points, then the segment must be intersected with the cell in a way
        // that some parts of the segment are outside of the cell
        else if (const std::vector<EK::Point_3> *pts = boost::get<std::vector<EK::Point_3>>(&site.first)) {
          for (const auto &p : *pts) {
            bool found = false;
            for (const auto &pr : allPoints) {
              if (std::get<1>(pr) == p) {
                found = true;
                break;
              }
            }

            if (!found) {
              allPoints.emplace_back(CGAL::squared_distance(p0, p), p);
            }
          }
        }
      }

      if (allPoints.size() == 2ull) {
        if (keep)
          faceFlags[fi] = ci;

        return false;
      }
      else {
        std::sort(allPoints.begin(), allPoints.end(),
          [](const std::tuple<EK::FT, EK::Point_3> &v1, const std::tuple<EK::FT, EK::Point_3> &v2) {
            return std::get<0>(v1) < std::get<0>(v2);
          });

        faceFlags[fi] = FS_DELETED;

        std::vector<int> vtxIDs(allPoints.size());
        vtxIDs[0] = medialAxisSkeletonTrianglesAndEdges[fi][0];
        vtxIDs.back() = medialAxisSkeletonTrianglesAndEdges[fi][1];

        for (size_t i = 1; i < allPoints.size() - 1; i++) {
          bool found = false;
          for (int vi = originalEnd; vi < (int)medialAxisVertices.size(); vi++) {
            if (std::get<1>(allPoints[i]) == medialAxisVertices[vi]) {
              vtxIDs[i] = vi;
              found = true;
              break;
            }
          }

          if (found == false) {
            medialAxisVertices.emplace_back(std::get<1>(allPoints[i]));
            vtxIDs[i] = (int)medialAxisVertices.size() - 1;
          }
        }

        for (size_t i = 0; i < allPoints.size() - 1; i++) {
          std::vector<int> fvid = { vtxIDs[i], vtxIDs[i + 1] };
          medialAxisSkeletonTrianglesAndEdges.emplace_back(fvid);
          faceFlags.emplace_back(FS_UNADDRESSED);
        }

        return true;
      }
    };

    auto breakFace = [&](int ci, int fi, bool keep) -> bool {
      Poly &cell = voronoiCellPolys[ci];
      // const EK::Point_3 &p0 = medialAxisVertices[medialAxisSkeletonTrianglesAndEdges[fi][0]];
      // const EK::Point_3 &p1 = medialAxisVertices[medialAxisSkeletonTrianglesAndEdges[fi][1]];
      // const EK::Point_3 &p2 = medialAxisVertices[medialAxisSkeletonTrianglesAndEdges[fi][2]];
      // I dont know why but the above code sometimes causes a crash
      const EK::Point_3 p0 = medialAxisVertices[medialAxisSkeletonTrianglesAndEdges[fi][0]];
      const EK::Point_3 p1 = medialAxisVertices[medialAxisSkeletonTrianglesAndEdges[fi][1]];
      const EK::Point_3 p2 = medialAxisVertices[medialAxisSkeletonTrianglesAndEdges[fi][2]];

      Poly mesh;
      mesh.make_triangle(p0, p1, p2);

      CGAL::Polygon_mesh_processing::corefine(mesh, cell,
        CGAL::parameters::do_not_modify(false),
        CGAL::parameters::do_not_modify(true));

      if (mesh.size_of_vertices() == 3ull) {
        if (keep)
          faceFlags[fi] = ci;

        return false;
      }
      else {
        std::map<int, int> vertexIDMap;
        int inc = 0;
        for (auto vit = mesh.vertices_begin(); vit != mesh.vertices_end(); ++vit) {
          vit->id() = inc++;

          if (vit->point() == p0) {
            vertexIDMap[vit->id()] = medialAxisSkeletonTrianglesAndEdges[fi][0];
          }
          else if (vit->point() == p1) {
            vertexIDMap[vit->id()] = medialAxisSkeletonTrianglesAndEdges[fi][1];
          }
          else if (vit->point() == p2) {
            vertexIDMap[vit->id()] = medialAxisSkeletonTrianglesAndEdges[fi][2];
          }
          else {
            if constexpr (1) {
              bool found = false;
              for (int vi = originalEnd; vi < (int)medialAxisVertices.size(); vi++) {
                if (vit->point() == medialAxisVertices[vi]) {
                  vertexIDMap[vit->id()] = vi;
                  found = true;
                  break;
                }
              }

              if (found == false) {
                medialAxisVertices.emplace_back(vit->point());
                vertexIDMap[vit->id()] = (int)medialAxisVertices.size() - 1;
              }
            }
            else {
              int parentID = fi;
              while (faceParents[parentID] >= 0) {
                parentID = faceParents[parentID];
              }

              auto cache_it = addedPointGlobalIndices.find(parentID);
              if (cache_it != addedPointGlobalIndices.end()) {
                bool found = false;
                for (const auto &pr : cache_it->second) {
                  if (vit->point() == std::get<0>(pr)) {
                    found = true;
                    vertexIDMap[vit->id()] = std::get<1>(pr);
                    break;
                  }
                }

                if (found == false) {
                  medialAxisVertices.emplace_back(vit->point());
                  vertexIDMap[vit->id()] = (int)medialAxisVertices.size() - 1;
                  cache_it->second.emplace_back(vit->point(), (int)medialAxisVertices.size() - 1);
                }
              }
              else {
                medialAxisVertices.emplace_back(vit->point());
                vertexIDMap[vit->id()] = (int)medialAxisVertices.size() - 1;
                auto ret = addedPointGlobalIndices.emplace(parentID, std::vector<std::tuple<EK::Point_3, int>>());
                ret.first->second.emplace_back(vit->point(), (int)medialAxisVertices.size() - 1);
              }
            }
          }
        }

        faceFlags[fi] = FS_DELETED;

        for (auto fit = mesh.facets_begin(); fit != mesh.facets_end(); ++fit) {
          Poly::Halfedge_handle h0 = fit->halfedge();
          Poly::Halfedge_handle h1 = h0->next();
          Poly::Halfedge_handle h2 = h1->next();

          std::vector<int> fvid(3);
          auto it = vertexIDMap.find(h0->vertex()->id());
          PGO_ALOG(it != vertexIDMap.end());
          fvid[0] = it->second;

          it = vertexIDMap.find(h1->vertex()->id());
          PGO_ALOG(it != vertexIDMap.end());
          fvid[1] = it->second;

          it = vertexIDMap.find(h2->vertex()->id());
          PGO_ALOG(it != vertexIDMap.end());
          fvid[2] = it->second;

          medialAxisSkeletonTrianglesAndEdges.emplace_back(fvid);
          faceFlags.emplace_back(FS_UNADDRESSED);
          faceParents.emplace_back(fi);
        }

        return true;
      }
    };

    faceFlags.resize(medialAxisSkeletonTrianglesAndEdges.size(), FS_UNADDRESSED);
    faceFlags.reserve(medialAxisSkeletonTrianglesAndEdges.size() * 10);

    for (int ci = 0; ci < (int)voronoiCellPolys.size(); ci++) {
      int offset = 0;
      while (1) {
        tbb::parallel_for(offset, (int)medialAxisSkeletonTrianglesAndEdges.size(), [&](int fi) {
          thread_local static std::vector<FaceGraphTree::Intersection_and_primitive_id<EK::Triangle_3>::Type> intersections;
          // for(int fi = offset; fi < (int)medialAxisSkeletonTrianglesAndEdges.size(); fi++) {
          //   static std::vector<FaceGraphTree::Intersection_and_primitive_id<EK::Triangle_3>::Type> intersections;

          if (faceFlags[fi] == FS_UNADDRESSED) {
            if (medialAxisSkeletonTrianglesAndEdges[fi].size() == 2ull) {
              // const EK::Point_3 &p0 = medialAxisVertices[medialAxisSkeletonTrianglesAndEdges[fi][0]];
              // const EK::Point_3 &p1 = medialAxisVertices[medialAxisSkeletonTrianglesAndEdges[fi][1]];
              const EK::Point_3 p0 = medialAxisVertices[medialAxisSkeletonTrianglesAndEdges[fi][0]];
              const EK::Point_3 p1 = medialAxisVertices[medialAxisSkeletonTrianglesAndEdges[fi][1]];

              // if all two point is not outside
              if (voronoiSideChecker[ci]->operator()(p0) != CGAL::ON_UNBOUNDED_SIDE &&
                voronoiSideChecker[ci]->operator()(p1) != CGAL::ON_UNBOUNDED_SIDE) {
                // check if there are any intersections
                EK::Segment_3 seg(p0, p1);
                intersections.clear();
                voronoiCellBVTree[ci]->all_intersections(seg, std::back_inserter(intersections));

                // no intersection
                // it means the entire triangle is inside
                if (intersections.size() == 0ull) {
                  faceFlags[fi] = ci;
                }
                else {
                  bool isInside = true;
                  // check intersection case
                  for (const auto &site : intersections) {
                    // if the intersection is a segment, it means the entire segment is on the boundary of the cell
                    if (const EK::Segment_3 *s = boost::get<EK::Segment_3>(&site.first)) {
                    }
                    // if the intersection is one point, it means the point is on the boundary of the cell but the others in the cell
                    else if (const EK::Point_3 *p = boost::get<EK::Point_3>(&site.first)) {
                    }
                    // if the intersection is more than one points, then the segment must be intersected with the cell in a way
                    // that some parts of the segment are outside of the cell
                    else if (const std::vector<EK::Point_3> *pts = boost::get<std::vector<EK::Point_3>>(&site.first)) {
                      if (pts->size() > 1ull) {
                        isInside = false;
                        break;
                      }
                    }
                  }

                  if (isInside) {
                    faceFlags[fi] = ci;
                  }
                  else {
                    faceFlags[fi] = FS_INTERSECTED;
                  }
                }
              }
              else {
                EK::Segment_3 seg(p0, p1);
                if (voronoiCellBVTree[ci]->do_intersect(seg)) {
                  faceFlags[fi] = FS_INTERSECTED;
                }
              }
            }
            else if (medialAxisSkeletonTrianglesAndEdges[fi].size() == 3ull) {
              // const EK::Point_3 &p0 = medialAxisVertices[medialAxisSkeletonTrianglesAndEdges[fi][0]];
              // const EK::Point_3 &p1 = medialAxisVertices[medialAxisSkeletonTrianglesAndEdges[fi][1]];
              // const EK::Point_3 &p2 = medialAxisVertices[medialAxisSkeletonTrianglesAndEdges[fi][2]];
              const EK::Point_3 p0 = medialAxisVertices[medialAxisSkeletonTrianglesAndEdges[fi][0]];
              const EK::Point_3 p1 = medialAxisVertices[medialAxisSkeletonTrianglesAndEdges[fi][1]];
              const EK::Point_3 p2 = medialAxisVertices[medialAxisSkeletonTrianglesAndEdges[fi][2]];

              // if all three point is not outside
              if (voronoiSideChecker[ci]->operator()(p0) != CGAL::ON_UNBOUNDED_SIDE &&
                voronoiSideChecker[ci]->operator()(p1) != CGAL::ON_UNBOUNDED_SIDE &&
                voronoiSideChecker[ci]->operator()(p2) != CGAL::ON_UNBOUNDED_SIDE) {
                // check if there are any intersections
                EK::Triangle_3 tri(p0, p1, p2);
                intersections.clear();
                voronoiCellBVTree[ci]->all_intersections(tri, std::back_inserter(intersections));

                // no intersection
                // it means the entire triangle is inside
                if (intersections.size() == 0ull) {
                  faceFlags[fi] = ci;
                }
                else {
                  bool isInside = true;
                  // check intersection case
                  for (const auto &site : intersections) {
                    Poly::Facet_const_handle f = site.second;
                    const auto &fp0 = f->halfedge()->vertex()->point();
                    const auto &fp1 = f->halfedge()->next()->vertex()->point();
                    const auto &fp2 = f->halfedge()->next()->next()->vertex()->point();

                    // if the intersection is a triangle, it means the entire triangle is on the boundary of the cell
                    if (const EK::Triangle_3 *t = boost::get<EK::Triangle_3>(&site.first)) {
                    }
                    // if the intersection is a segment, it means the entire segment is on the boundary of the cell
                    else if (const EK::Segment_3 *s = boost::get<EK::Segment_3>(&site.first)) {
                    }
                    // if the intersection is one point, it means the point is on the boundary of the cell but others are in the cell
                    else if (const EK::Point_3 *p = boost::get<EK::Point_3>(&site.first)) {
                    }
                    // if the intersection is more than one points, then the triangle must be intersected with the cell in a way
                    // that some parts of the triangle are outside of the cell
                    else if (const std::vector<EK::Point_3> *pts = boost::get<std::vector<EK::Point_3>>(&site.first)) {
                      if (pts->size() > 1ull) {
                        isInside = false;
                        break;
                      }
                    }
                  }

                  if (isInside) {
                    faceFlags[fi] = ci;
                  }
                  else {
                    faceFlags[fi] = FS_INTERSECTED;
                  }
                }
              }
              else {
                EK::Triangle_3 tri(p0, p1, p2);
                if (voronoiCellBVTree[ci]->do_intersect(tri)) {
                  faceFlags[fi] = FS_INTERSECTED;
                }
              }
            }
            else {
              throw std::runtime_error("impossible");
            }
          }
        });

        int oldSize = (int)medialAxisSkeletonTrianglesAndEdges.size();

        int counter = 0;
        for (int fi = offset; fi < (int)medialAxisSkeletonTrianglesAndEdges.size(); fi++) {
          if (faceFlags[fi] == FS_INTERSECTED) {
            bool ret = false;
            std::vector<FaceGraphTree::Intersection_and_primitive_id<EK::Triangle_3>::Type> intersections;
            if (medialAxisSkeletonTrianglesAndEdges[fi].size() == 2ull) {
              // std::cout << "Encounter edge skeleton." << std::endl;
              ret = breakSeg(ci, fi, false);
            }
            else if (medialAxisSkeletonTrianglesAndEdges[fi].size() == 3ull) {
              ret = breakFace(ci, fi, false);
            }
            else {
              throw std::runtime_error("impossible");
            }

            if (ret)
              counter++;
            else {
              faceFlags[fi] = FS_UNADDRESSED;
            }
          }
        }

        if (counter == 0)
          break;

        offset = oldSize;
      }
    }

    // save
    // std::vector<EK::Point_3> medialAxisVertices;
    // std::vector<std::vector<int>> medialAxisSkeletonTrianglesAndEdges;
    // std::vector<int> faceFlags(medialAxisSkeletonTrianglesAndEdges.size(), FS_UNADDRESSED);
    if constexpr (0) {
      nlohmann::json jout;
      for (int i = 0; i < (int)medialAxisVertices.size(); i++) {
        std::string key = fmt::format("{}", i);
        std::vector<char> val[3];
        for (int j = 0; j < 3; j++) {
          auto &op = medialAxisVertices[i][j].mpq();
          val[j].resize((mpz_sizeinbase(mpq_numref(op), 10) + mpz_sizeinbase(mpq_denref(op), 10) + 3) * 2);
          memset(val[j].data(), 0, sizeof(char) * val[j].size());
          mpq_get_str(val[j].data(), 10, op);
        }
        jout["medialAxisVertices"][key] = { val[0], val[1], val[2] };
      }
      for (int i = 0; i < (int)medialAxisSkeletonTrianglesAndEdges.size(); i++) {
        std::string key = fmt::format("{}", i);
        jout["medialAxisSkeletonTrianglesAndEdges"][key] = medialAxisSkeletonTrianglesAndEdges[i];
      }
      for (int i = 0; i < (int)faceFlags.size(); i++) {
        std::string key = fmt::format("{}", i);
        jout["faceFlags"][key] = faceFlags[i];
      }
      // std::ofstream(jsonIn1.c_str()) << jout.dump(2);
      // fmt::print("saved to {}\n", jsonIn1);
    }
  }

#if 0
  if constexpr (0) {
    std::ifstream ifs(jsonIn1.c_str());
    nlohmann::json jin;
    ifs >> jin;
    medialAxisVertices.resize(jin["medialAxisVertices"].size());
    for (auto it = jin["medialAxisVertices"].begin(); it != jin["medialAxisVertices"].end(); it++) {
      int idx = std::stoi(it.key());
      std::vector<char> valStr[3];
      for (int i = 0; i < 3; i++) {
        valStr[i] = it.value()[i];
      }
      medialAxisVertices[idx] = EK::Point_3(CGAL::Gmpq(valStr[0].data(), 10),
        CGAL::Gmpq(valStr[1].data(), 10),
        CGAL::Gmpq(valStr[2].data(), 10));
    }
    medialAxisSkeletonTrianglesAndEdges.resize(jin["medialAxisSkeletonTrianglesAndEdges"].size());
    for (auto it = jin["medialAxisSkeletonTrianglesAndEdges"].begin(); it != jin["medialAxisSkeletonTrianglesAndEdges"].end(); it++) {
      int idx = std::stoi(it.key());
      medialAxisSkeletonTrianglesAndEdges[idx].resize(it.value().size());
      for (int i = 0; i < (int)it.value().size(); i++) {
        medialAxisSkeletonTrianglesAndEdges[idx][i] = it.value()[i];
      }
    }
    faceFlags.resize(jin["faceFlags"].size());
    for (auto it = jin["faceFlags"].begin(); it != jin["faceFlags"].end(); it++) {
      int idx = std::stoi(it.key());
      faceFlags[idx] = it.value();
    }
    fmt::print("load from {}\n", jsonIn1);
  }
#endif

  if constexpr (dumpMesh) {
    tbb::parallel_for(0, (int)samplePoints.size(), [&](int ci) {
      // for (int ci = 0; ci < (int)samplePoints.size(); ci++) {
      pgo::Mesh::TriMeshGeo mm;
      for (size_t i = 0; i < medialAxisVertices.size(); i++) {
        ES::V3d p = toVec3(medialAxisVertices[i]);
        mm.addPos(p);
      }

      for (int fi = 0; fi < (int)medialAxisSkeletonTrianglesAndEdges.size(); fi++) {
        if (faceFlags[fi] == ci) {
          if (medialAxisSkeletonTrianglesAndEdges[fi].size() == 3ull) {
            mm.addTri(ES::V3i(medialAxisSkeletonTrianglesAndEdges[fi][0], medialAxisSkeletonTrianglesAndEdges[fi][1], medialAxisSkeletonTrianglesAndEdges[fi][2]));
          }
          else if (medialAxisSkeletonTrianglesAndEdges[fi].size() == 2ull) {
            mm.addPos(mm.pos(medialAxisSkeletonTrianglesAndEdges[fi][0]) + pgo::asVec3d(1e-6));
            mm.addTri(ES::V3i(medialAxisSkeletonTrianglesAndEdges[fi][0], medialAxisSkeletonTrianglesAndEdges[fi][1], mm.numVertices() - 1));
          }
          else {
            abort();
          }
        }
      }
      mm.save(fmt::format("a{}.obj", ci));
    });
  }

  HClockPt t6 = HClock::now();
  SPDLOG_LOGGER_INFO(pgo::Logging::lgr(), "cutting: {}s", duraSecond(t5, t6));
  //////

  // post process medialAxisSkeletonTrianglesAndEdges, remove medialAxisSkeletonTrianglesAndEdges that has (negative flag and COULD be WRONG) duplicated triangles/edges
  std::unordered_map<std::array<int, 3>, int, ES::IntArrayHash<3>, ES::IntArrayEqual<3>> newMedialAxisSkeletonTriangles;
  std::unordered_map<std::array<int, 2>, int, ES::IntArrayHash<2>, ES::IntArrayEqual<2>> newMedialAxisSkeletonEdges;
  std::vector<int> newFaceFlags;
  int idx = 0;
  for (int i = 0; i < (int)medialAxisSkeletonTrianglesAndEdges.size(); i++) {
    if (faceFlags[i] >= 0) {
      std::vector<int> tri = medialAxisSkeletonTrianglesAndEdges[i];
      std::sort(tri.begin(), tri.end());
      if ((int)tri.size() == 2) {
        std::array<int, 2> edge = { tri[0], tri[1] };
        if (newMedialAxisSkeletonEdges.find(edge) == newMedialAxisSkeletonEdges.end()) {
          newMedialAxisSkeletonEdges[edge] = idx;
          newFaceFlags.emplace_back(faceFlags[i]);
          idx++;
        }
      }
      else if ((int)tri.size() == 3) {
        std::array<int, 3> triangle = { tri[0], tri[1], tri[2] };
        if (newMedialAxisSkeletonTriangles.find(triangle) == newMedialAxisSkeletonTriangles.end()) {
          newMedialAxisSkeletonTriangles[triangle] = idx;
          newFaceFlags.emplace_back(faceFlags[i]);
          idx++;
        }
      }
      else {
        PGO_ALOG(false);
      }
    }
  }
  PGO_ALOG(newMedialAxisSkeletonTriangles.size() + newMedialAxisSkeletonEdges.size() == newFaceFlags.size());
  PGO_ALOG((int)newFaceFlags.size() == idx);
  std::vector<std::vector<int>> newMedialAxisSkeletonTrianglesAndEdges(newMedialAxisSkeletonTriangles.size() + newMedialAxisSkeletonEdges.size());

  for (auto &it : newMedialAxisSkeletonTriangles) {
    std::vector<int> &vec = newMedialAxisSkeletonTrianglesAndEdges[it.second];
    vec = { it.first[0], it.first[1], it.first[2] };
  }

  for (auto &it : newMedialAxisSkeletonEdges) {
    std::vector<int> &vec = newMedialAxisSkeletonTrianglesAndEdges[it.second];
    vec = { it.first[0], it.first[1] };
  }

  std::swap(medialAxisSkeletonTrianglesAndEdges, newMedialAxisSkeletonTrianglesAndEdges);
  std::swap(faceFlags, newFaceFlags);

  std::unordered_map<int, std::vector<int>> vtxTriangleNeighbors;
  for (int ti = 0; ti < (int)medialAxisSkeletonTrianglesAndEdges.size(); ti++) {
    PGO_ALOG(medialAxisSkeletonTrianglesAndEdges[ti].size() == 3ull || medialAxisSkeletonTrianglesAndEdges[ti].size() == 2ull);

    for (int j = 0; j < (int)medialAxisSkeletonTrianglesAndEdges[ti].size(); j++) {
      auto it = vtxTriangleNeighbors.find(medialAxisSkeletonTrianglesAndEdges[ti][j]);
      if (it != vtxTriangleNeighbors.end()) {
        it->second.emplace_back(ti);
      }
      else {
        vtxTriangleNeighbors.emplace(medialAxisSkeletonTrianglesAndEdges[ti][j], std::vector<int>{ ti });
      }
    }
  }

  std::unordered_map<int, std::vector<int>> triangleVtxNeighbors;
  for (const auto &pr : vtxTriangleNeighbors) {
    for (int i = 0; i < (int)pr.second.size(); i++) {
      for (int j = i + 1; j < (int)pr.second.size(); j++) {
        triangleVtxNeighbors[pr.second[i]].emplace_back(pr.second[j]);
        triangleVtxNeighbors[pr.second[j]].emplace_back(pr.second[i]);
      }
    }
  }
  for (auto &pr : triangleVtxNeighbors) {
    pgo::BasicAlgorithms::sortAndDeduplicateWithErase(pr.second);
  }

  // remove triangles that has no neighbors
  std::vector<int> removeTriIDs;
  if (triangleVtxNeighbors.size() != medialAxisSkeletonTrianglesAndEdges.size()) {
    for (int i = 0; i < (int)medialAxisSkeletonTrianglesAndEdges.size(); i++) {
      if (triangleVtxNeighbors.find(i) == triangleVtxNeighbors.end()) {
        std::cout << "triangle " << i << " has no neighbors." << std::endl;
        removeTriIDs.emplace_back(i);

        if constexpr (dumpMesh) {
          pgo::Mesh::TriMeshGeo debugMesh;
          if ((int)medialAxisSkeletonTrianglesAndEdges[i].size() == 2) {
            debugMesh.addPos(toVec3(medialAxisVertices[medialAxisSkeletonTrianglesAndEdges[i][0]]));
            debugMesh.addPos(toVec3(medialAxisVertices[medialAxisSkeletonTrianglesAndEdges[i][1]]));
            debugMesh.addPos(toVec3(medialAxisVertices[medialAxisSkeletonTrianglesAndEdges[i][1]]) + pgo::asVec3d(1e-6));
          }
          else {
            debugMesh.addPos(toVec3(medialAxisVertices[medialAxisSkeletonTrianglesAndEdges[i][0]]));
            debugMesh.addPos(toVec3(medialAxisVertices[medialAxisSkeletonTrianglesAndEdges[i][1]]));
            debugMesh.addPos(toVec3(medialAxisVertices[medialAxisSkeletonTrianglesAndEdges[i][2]]));
          }
          debugMesh.addTri(ES::V3i(0, 1, 2));
          debugMesh.save("debug.obj");
        }
      }
    }
  }

  std::unordered_map<int, std::vector<int>> triangleNeighbors;
  if (removeTriIDs.size() > 0) {
    // erase triangles that has no neighborrs in medialAxisSkeletonTrianglesAndEdges and newFaceFlags
    pgo::BasicAlgorithms::sortAndDeduplicateWithErase(removeTriIDs);
    std::cout << "remove " << removeTriIDs.size() << " triangles." << std::endl;
    std::cout << "before remove: " << medialAxisSkeletonTrianglesAndEdges.size() << std::endl;
    for (int i = (int)removeTriIDs.size() - 1; i >= 0; i--) {
      int triID = removeTriIDs[i];
      medialAxisSkeletonTrianglesAndEdges.erase(medialAxisSkeletonTrianglesAndEdges.begin() + triID);
      faceFlags.erase(faceFlags.begin() + triID);
    }
    std::cout << "after remove: " << medialAxisSkeletonTrianglesAndEdges.size() << std::endl;

    // recalculate neigbors
    std::unordered_map<int, std::vector<int>> vtxTriangleFinalNeighbors;
    for (int ti = 0; ti < (int)medialAxisSkeletonTrianglesAndEdges.size(); ti++) {
      PGO_ALOG(medialAxisSkeletonTrianglesAndEdges[ti].size() == 3ull || medialAxisSkeletonTrianglesAndEdges[ti].size() == 2ull);
      for (int j = 0; j < (int)medialAxisSkeletonTrianglesAndEdges[ti].size(); j++) {
        auto it = vtxTriangleFinalNeighbors.find(medialAxisSkeletonTrianglesAndEdges[ti][j]);
        if (it != vtxTriangleFinalNeighbors.end()) {
          it->second.emplace_back(ti);
        }
        else {
          vtxTriangleFinalNeighbors.emplace(medialAxisSkeletonTrianglesAndEdges[ti][j], std::vector<int>{ ti });
        }
      }
    }

    for (const auto &pr : vtxTriangleFinalNeighbors) {
      for (int i = 0; i < (int)pr.second.size(); i++) {
        for (int j = i + 1; j < (int)pr.second.size(); j++) {
          triangleNeighbors[pr.second[i]].emplace_back(pr.second[j]);
          triangleNeighbors[pr.second[j]].emplace_back(pr.second[i]);
        }
      }
    }
    for (auto &pr : triangleNeighbors) {
      pgo::BasicAlgorithms::sortAndDeduplicateWithErase(pr.second);
    }

    if constexpr (dumpMesh) {
      tbb::parallel_for(0, (int)samplePoints.size(), [&](int ci) {
        std::ofstream tempfile(fmt::format("aFinal{}.obj", ci));
        pgo::Mesh::TriMeshGeo dumpMesh_a;
        for (int fi = 0; fi < (int)medialAxisSkeletonTrianglesAndEdges.size(); fi++) {
          if (faceFlags[fi] == ci) {
            if ((int)medialAxisSkeletonTrianglesAndEdges[fi].size() == 2) {
              dumpMesh_a.addPos(toVec3(medialAxisVertices[medialAxisSkeletonTrianglesAndEdges[fi][0]]));
              dumpMesh_a.addPos(toVec3(medialAxisVertices[medialAxisSkeletonTrianglesAndEdges[fi][1]]));
              dumpMesh_a.addPos(toVec3(medialAxisVertices[medialAxisSkeletonTrianglesAndEdges[fi][1]]) + pgo::asVec3d(1e-6));
            }
            else {
              dumpMesh_a.addPos(toVec3(medialAxisVertices[medialAxisSkeletonTrianglesAndEdges[fi][0]]));
              dumpMesh_a.addPos(toVec3(medialAxisVertices[medialAxisSkeletonTrianglesAndEdges[fi][1]]));
              dumpMesh_a.addPos(toVec3(medialAxisVertices[medialAxisSkeletonTrianglesAndEdges[fi][2]]));
            }
            dumpMesh_a.addTri(ES::V3i(dumpMesh_a.numVertices() - 3, dumpMesh_a.numVertices() - 2, dumpMesh_a.numVertices() - 1));
          }
        }
        dumpMesh_a.save(fmt::format("aFinal{}.obj", ci));
      });
    }
  }
  else {
    triangleNeighbors = triangleVtxNeighbors;
  }

  // std::unordered_map<int, std::vector<int>> triangleNeighbors = triangleVtxNeighbors;
  std::cout << "#tri: " << triangleNeighbors.size() << std::endl;

  //////////
  //////////
  using Tri = EK::Triangle_3;
  using Seg = EK::Segment_3;
  using TrianagleIterator = std::vector<Tri>::iterator;
  using TrianaglePrimitive = CGAL::AABB_triangle_primitive<EK, TrianagleIterator>;
  using AABB_triangle_traits = CGAL::AABB_traits<EK, TrianaglePrimitive>;
  using TriangleTree = CGAL::AABB_tree<AABB_triangle_traits>;

  // triangles
  std::vector<int> seedTriangles(samplePoints.size());
  tbb::parallel_for(0, (int)samplePoints.size(), [&](int si) {
    EK::Point_3 samplePt(samplePoints[si][0], samplePoints[si][1], samplePoints[si][2]);

    EK::FT minDist = 100;
    int selTriID = -1;
    for (int ti = 0; ti < (int)medialAxisSkeletonTrianglesAndEdges.size(); ti++) {
      if (faceFlags[ti] != si) {
        continue;
      }

      EK::FT dist;
      if ((int)medialAxisSkeletonTrianglesAndEdges[ti].size() == 2) {
        Seg seg(medialAxisVertices[medialAxisSkeletonTrianglesAndEdges[ti][0]],
          medialAxisVertices[medialAxisSkeletonTrianglesAndEdges[ti][1]]);
        dist = CGAL::squared_distance(samplePt, seg);
      }
      else {
        PGO_ALOG((int)medialAxisSkeletonTrianglesAndEdges[ti].size() == 3);
        Tri tri(medialAxisVertices[medialAxisSkeletonTrianglesAndEdges[ti][0]],
          medialAxisVertices[medialAxisSkeletonTrianglesAndEdges[ti][1]],
          medialAxisVertices[medialAxisSkeletonTrianglesAndEdges[ti][2]]);
        dist = CGAL::squared_distance(samplePt, tri);
      }

      if (minDist > dist) {
        minDist = dist;
        selTriID = ti;
      }
    }
    // if (selTriID == -1) {
    //   std::cout << "no triangle selected for seed " << si << std::endl;
    //   pgo::Mesh::TriMeshGeo triMeshDebug0;
    //   triMeshDebug0.addPos(toVec3(samplePt));
    //   triMeshDebug0.save("noTriangleSelected.obj");
    // }
    seedTriangles[si] = selTriID;
    // std::cout << "s" << si << ": " << std::sqrt(CGAL::to_double(minDist)) << std::endl;
  });

  HClockPt t7 = HClock::now();
  SPDLOG_LOGGER_INFO(pgo::Logging::lgr(), "dist: {}s", duraSecond(t6, t7));

  if constexpr (dumpMesh) {
    // Debug faceFlags
    pgo::Mesh::TriMeshGeo triMeshDebug0;
    for (int i = 0; i < (int)medialAxisSkeletonTrianglesAndEdges.size(); i++) {
      if (faceFlags[i] < 0) {
        std::cout << "face flag" << i << " is not selected." << std::endl;
      }
      if ((int)medialAxisSkeletonTrianglesAndEdges[i].size() == 2) {
        triMeshDebug0.addMesh(pgo::Mesh::createSingleTriangleMesh(
          toVec3(medialAxisVertices[medialAxisSkeletonTrianglesAndEdges[i][0]]),
          toVec3(medialAxisVertices[medialAxisSkeletonTrianglesAndEdges[i][1]]),
          toVec3(medialAxisVertices[medialAxisSkeletonTrianglesAndEdges[i][1]]) + pgo::asVec3d(1e-6)));
      }
      else {
        triMeshDebug0.addMesh(pgo::Mesh::createSingleTriangleMesh(
          toVec3(medialAxisVertices[medialAxisSkeletonTrianglesAndEdges[i][0]]),
          toVec3(medialAxisVertices[medialAxisSkeletonTrianglesAndEdges[i][1]]),
          toVec3(medialAxisVertices[medialAxisSkeletonTrianglesAndEdges[i][2]])));
      }
    }
    triMeshDebug0.save("debugFaceFlags.obj");
  }

  // find one connected component for each seed triangle
  std::vector<int> faceSelFlags(medialAxisSkeletonTrianglesAndEdges.size(), -1);
  std::vector<int> Q;
  Q.reserve(medialAxisSkeletonTrianglesAndEdges.size());

  size_t start = 0;
  for (int si = 0; si < (int)samplePoints.size(); si++) {
    start = 0;
    Q.clear();

    Q.emplace_back(seedTriangles[si]);
    faceSelFlags[seedTriangles[si]] = si;

    while (Q.size() > start) {
      int curID = Q[start++];

      auto it = triangleNeighbors.find(curID);
      PGO_ALOG(it != triangleNeighbors.end());

      for (const auto &ntri : it->second) {
        if (faceSelFlags[ntri] >= 0)
          continue;

        if (faceFlags[ntri] != faceFlags[seedTriangles[si]])
          continue;

        Q.emplace_back(ntri);
        faceSelFlags[ntri] = si;
      }
    }
  }

  std::vector<std::vector<double>> distances(samplePoints.size());
  std::vector<std::vector<int>> prevTriangleIDs(samplePoints.size());

  std::vector<ES::V3d> triangleCenters(medialAxisSkeletonTrianglesAndEdges.size());
  std::vector<EK::Point_3> triangleCentersEK(medialAxisSkeletonTrianglesAndEdges.size());

  tbb::parallel_for(0, (int)medialAxisSkeletonTrianglesAndEdges.size(), [&](int ti) {
    if ((int)medialAxisSkeletonTrianglesAndEdges[ti].size() == 2) {
      EK::Point_3 pEK0 = medialAxisVertices[medialAxisSkeletonTrianglesAndEdges[ti][0]];
      EK::Point_3 pEK1 = medialAxisVertices[medialAxisSkeletonTrianglesAndEdges[ti][1]];
      triangleCentersEK[ti] = CGAL::midpoint(pEK0, pEK1);
      triangleCenters[ti] = toVec3(triangleCentersEK[ti]);
    }
    else {
      PGO_ALOG((int)medialAxisSkeletonTrianglesAndEdges[ti].size() == 3);
      EK::Point_3 pEK0 = medialAxisVertices[medialAxisSkeletonTrianglesAndEdges[ti][0]];
      EK::Point_3 pEK1 = medialAxisVertices[medialAxisSkeletonTrianglesAndEdges[ti][1]];
      EK::Point_3 pEK2 = medialAxisVertices[medialAxisSkeletonTrianglesAndEdges[ti][2]];
      triangleCentersEK[ti] = CGAL::centroid(pEK0, pEK1, pEK2);
      triangleCenters[ti] = toVec3(triangleCentersEK[ti]);
    }
  });

  // for (int si = 0; si < (int)samplePoints.size(); si++) {
  tbb::parallel_for(0, (int)samplePoints.size(), [&](int si) {
    int srcTriID = seedTriangles[si];
    std::vector<int> &prev = prevTriangleIDs[si];
    prev.assign(medialAxisSkeletonTrianglesAndEdges.size(), -1);

    std::vector<double> &dist = distances[si];
    dist.assign(medialAxisSkeletonTrianglesAndEdges.size(), std::numeric_limits<double>::max());
    dist[srcTriID] = 0;

    std::priority_queue<std::pair<double, int>, std::vector<std::pair<double, int>>, std::greater<std::pair<double, int>>> Q;
    Q.emplace(0, srcTriID);

    while (!Q.empty()) {
      int curID = Q.top().second;
      ES::V3d centerCurID = triangleCenters[curID];

      Q.pop();

      auto it = triangleNeighbors.find(curID);
      PGO_ALOG(it != triangleNeighbors.end());

      for (const auto &ntri : it->second) {
        if (ntri == curID)
          continue;

        ES::V3d centerNextID = triangleCenters[ntri];
        double weight = (centerCurID - centerNextID).norm();

        // ES::V3d midPt = (centerCurID + centerNextID) * 0.5;
        // auto queryRet = targetMeshBVTree.closestTriangleQuery(targetMesh, midPt);
        double w_scaled = 1.0;  // / (std::sqrt(queryRet.dist2) + 1e-4);

        weight *= w_scaled;

        if (dist[ntri] > dist[curID] + weight) {
          dist[ntri] = dist[curID] + weight;
          prev[ntri] = curID;
          Q.emplace(dist[ntri], ntri);
        }
      }
    }
  });

  HClockPt t8 = HClock::now();
  SPDLOG_LOGGER_INFO(pgo::Logging::lgr(), "graph shortest distance {}s", duraSecond(t7, t8));

  // post processing faceSelFlags
  // Wrong, resulting in disconnected components
  // for (int i = 0; i < (int)medialAxisSkeletonTrianglesAndEdges.size(); i++) {
  //   if (faceSelFlags[i] < 0) {
  //     double minDist = std::numeric_limits<double>::max();
  //     int minDistID = -1;
  //     for (int si = 0; si < (int)samplePoints.size(); si++) {
  //       if (distances[si][i] < minDist) {
  //         minDist = distances[si][i];
  //         minDistID = si;
  //       }
  //     }
  //     faceSelFlags[i] = minDistID;
  //   }
  // }

  std::vector<int> negFaceSelFlagTris;
  negFaceSelFlagTris.reserve(medialAxisSkeletonTrianglesAndEdges.size());
  for (int i = 0; i < (int)medialAxisSkeletonTrianglesAndEdges.size(); i++) {
    if (faceSelFlags[i] < 0) {
      negFaceSelFlagTris.emplace_back(i);
    }
  }

  // for every connected component that is not flagged, assign it to the nearest seed triangle that is first encountered by BFS
  while ((int)negFaceSelFlagTris.size() > 0) {
    start = 0;
    Q.clear();
    int selTri = negFaceSelFlagTris.back();
    Q.emplace_back(selTri);

    int encounteredFlag = -1;
    std::vector<int> allNeighborTypes;
    std::vector<int> allNeighborFlags;
    while (Q.size() > start) {
      int curID = Q[start++];
      auto it = triangleNeighbors.find(curID);
      PGO_ALOG(it != triangleNeighbors.end());

      for (const auto &ntri : it->second) {
        std::set<int> vtxSet0;
        for (int vi : medialAxisSkeletonTrianglesAndEdges[curID]) {
          vtxSet0.insert(vi);
        }
        std::set<int> vtxSet1;
        for (int vi : medialAxisSkeletonTrianglesAndEdges[ntri]) {
          vtxSet1.insert(vi);
        }

        std::vector<int> sharedVtx;
        std::set_intersection(vtxSet0.begin(), vtxSet0.end(), vtxSet1.begin(), vtxSet1.end(), std::back_inserter(sharedVtx));

        PGO_ALOG((int)sharedVtx.size() == 1 || (int)sharedVtx.size() == 2);

        if (faceSelFlags[ntri] >= 0) {
          allNeighborTypes.emplace_back((int)sharedVtx.size());
          allNeighborFlags.emplace_back(faceSelFlags[ntri]);
        }

        if (faceSelFlags[ntri] >= 0 && (int)sharedVtx.size() > 1) {
          // only when it is edge neighbor
          encounteredFlag = faceSelFlags[ntri];
          break;
        }
        else {
          if (faceSelFlags[ntri] < 0) {
            // only when it is not flagged
            if (std::find(Q.begin(), Q.end(), ntri) == Q.end()) {
              Q.emplace_back(ntri);
            }
          }
        }
      }

      if (encounteredFlag >= 0) {
        break;
      }
    }

    if (encounteredFlag < 0) {
      // check if all neighbors are vertex neighbors
      bool isAllVtxNeighbor = true;
      for (int i = 0; i < (int)allNeighborTypes.size(); i++) {
        if (allNeighborTypes[i] == 2) {
          isAllVtxNeighbor = false;
          break;
        }
      }
      if (isAllVtxNeighbor) {
        // assign it to the first encountered flag
        encounteredFlag = allNeighborFlags[0];
      }
    }

    if (encounteredFlag < 0) {
      pgo::Mesh::TriMeshGeo debugMesh;
      for (int i = 0; i < (int)Q.size(); i++) {
        if ((int)medialAxisSkeletonTrianglesAndEdges[Q[i]].size() == 2) {
          debugMesh.addMesh(pgo::Mesh::createSingleTriangleMesh(toVec3(medialAxisVertices[medialAxisSkeletonTrianglesAndEdges[Q[i]][0]]),
            toVec3(medialAxisVertices[medialAxisSkeletonTrianglesAndEdges[Q[i]][1]]),
            toVec3(medialAxisVertices[medialAxisSkeletonTrianglesAndEdges[Q[i]][0]]) + pgo::asVec3d(1e-6)));
        }
        else {
          debugMesh.addMesh(pgo::Mesh::createSingleTriangleMesh(toVec3(medialAxisVertices[medialAxisSkeletonTrianglesAndEdges[Q[i]][0]]),
            toVec3(medialAxisVertices[medialAxisSkeletonTrianglesAndEdges[Q[i]][1]]),
            toVec3(medialAxisVertices[medialAxisSkeletonTrianglesAndEdges[Q[i]][2]])));
        }
      }
      debugMesh.save("debug.obj");
    }

    PGO_ALOG(encounteredFlag >= 0);
    for (int i = 0; i < (int)Q.size(); i++) {
      faceSelFlags[Q[i]] = encounteredFlag;
    }
    // remove all the triangles that are already processed
    for (int i = 0; i < (int)negFaceSelFlagTris.size(); i++) {
      if (faceSelFlags[negFaceSelFlagTris[i]] >= 0) {
        negFaceSelFlagTris.erase(negFaceSelFlagTris.begin() + i);
        i--;
      }
    }
  }

  if (cellCenters) {
    cellCenters->resize(samplePoints.size());
    for (int si = 0; si < (int)samplePoints.size(); si++) {
      ES::V3d center(0, 0, 0);
      int count = 0;
      for (int ti = 0; ti < (int)medialAxisSkeletonTrianglesAndEdges.size(); ti++) {
        if (faceSelFlags[ti] == si) {
          center += triangleCenters[ti];
          count += 1;
        }
      }
      PGO_ALOG(count > 0);

      center /= count;
      cellCenters->at(si) = ES::V3d(center[0], center[1], center[2]);
    }
  }

  if constexpr (dumpMesh) {
    // Debug faceSelFlags
    pgo::Mesh::TriMeshGeo triMeshDebug;
    for (int i = 0; i < (int)medialAxisSkeletonTrianglesAndEdges.size(); i++) {
      if (faceSelFlags[i] < 0) {
        std::cout << "face selflag" << i << " is not selected." << std::endl;
      }
      if ((int)medialAxisSkeletonTrianglesAndEdges[i].size() == 2) {
        triMeshDebug.addMesh(pgo::Mesh::createSingleTriangleMesh(
          toVec3(medialAxisVertices[medialAxisSkeletonTrianglesAndEdges[i][0]]),
          toVec3(medialAxisVertices[medialAxisSkeletonTrianglesAndEdges[i][1]]),
          toVec3(medialAxisVertices[medialAxisSkeletonTrianglesAndEdges[i][0]]) + pgo::asVec3d(1e-6)));
      }
      else {
        triMeshDebug.addMesh(pgo::Mesh::createSingleTriangleMesh(
          toVec3(medialAxisVertices[medialAxisSkeletonTrianglesAndEdges[i][0]]),
          toVec3(medialAxisVertices[medialAxisSkeletonTrianglesAndEdges[i][1]]),
          toVec3(medialAxisVertices[medialAxisSkeletonTrianglesAndEdges[i][2]])));
      }
    }
    triMeshDebug.save("debugFaceSelFlags.obj");

    tbb::parallel_for(0, (int)samplePoints.size(), [&](int si) {
      pgo::Mesh::TriMeshGeo debugMesh_s;
      for (int fi = 0; fi < (int)medialAxisSkeletonTrianglesAndEdges.size(); fi++) {
        if (faceSelFlags[fi] == si) {
          if ((int)medialAxisSkeletonTrianglesAndEdges[fi].size() == 2) {
            debugMesh_s.addMesh(pgo::Mesh::createSingleTriangleMesh(toVec3(medialAxisVertices[medialAxisSkeletonTrianglesAndEdges[fi][0]]),
              toVec3(medialAxisVertices[medialAxisSkeletonTrianglesAndEdges[fi][1]]),
              toVec3(medialAxisVertices[medialAxisSkeletonTrianglesAndEdges[fi][0]]) + pgo::asVec3d(1e-6)));
          }
          else {
            debugMesh_s.addMesh(pgo::Mesh::createSingleTriangleMesh(toVec3(medialAxisVertices[medialAxisSkeletonTrianglesAndEdges[fi][0]]),
              toVec3(medialAxisVertices[medialAxisSkeletonTrianglesAndEdges[fi][1]]),
              toVec3(medialAxisVertices[medialAxisSkeletonTrianglesAndEdges[fi][2]])));
          }
        }
      }
      debugMesh_s.save(fmt::format("s{}.obj", si));
    });
  }

  // Now recalculate the voronoi diagram for all connected components, convert to triangulation for adding edges and triangles
  // neighbors of each connected component, which could be different from delauny neighbors
  // find all the edges of Voronoi cells
  std::vector<std::array<int, 2>> allDelaunayEdges;
  allDelaunayEdges.reserve((int)samplePoints.size() * (int)samplePoints.size());
  for (int si = 0; si < (int)samplePoints.size(); si++) {
    for (int ti = 0; ti < (int)medialAxisSkeletonTrianglesAndEdges.size(); ti++) {
      if (faceSelFlags[ti] != si)
        continue;

      auto it = triangleNeighbors.find(ti);
      PGO_ALOG(it != triangleNeighbors.end());
      for (const auto &ntri : it->second) {
        PGO_ALOG(faceSelFlags[ntri] >= 0);

        if (faceSelFlags[ntri] == si)
          continue;

        int vtx0 = si, vtx1 = faceSelFlags[ntri];
        if (vtx0 > vtx1)
          std::swap(vtx0, vtx1);

        allDelaunayEdges.emplace_back(std::array<int, 2>{ vtx0, vtx1 });
      }
    }
  }

  pgo::BasicAlgorithms::sortAndDeduplicateWithErase(allDelaunayEdges);

  std::map<int, std::vector<int>> maVtxNeigborSamplePts;
  for (int ti = 0; ti < (int)medialAxisSkeletonTrianglesAndEdges.size(); ti++) {
    for (int vi : medialAxisSkeletonTrianglesAndEdges[ti]) {
      maVtxNeigborSamplePts[vi].emplace_back(faceSelFlags[ti]);
    }
  }

  for (auto &v : maVtxNeigborSamplePts) {
    pgo::BasicAlgorithms::sortAndDeduplicateWithErase(v.second);
  }

  // it is possible a vertex is connected to more than 3 sample points, in this case, we need to split the triangle
  std::map<int, std::vector<int>> maVtx2SkeletonSurface;
  for (auto &v : maVtxNeigborSamplePts) {
    if ((int)v.second.size() < 3) {
      if ((int)v.second.size() == 2) {
        int vid0 = v.second[0], vid1 = v.second[1];
        if (vid0 > vid1)
          std::swap(vid0, vid1);
        PGO_ALOG(std::find(allDelaunayEdges.begin(), allDelaunayEdges.end(), std::array<int, 2>{ vid0, vid1 }) != allDelaunayEdges.end());
      }
      continue;
    }
    if ((int)v.second.size() == 3) {
      for (int j = 0; j < 3; j++) {
        int vid0 = v.second[j], vid1 = v.second[(j + 1) % 3];
        if (vid0 > vid1)
          std::swap(vid0, vid1);
        PGO_ALOG(std::find(allDelaunayEdges.begin(), allDelaunayEdges.end(), std::array<int, 2>{ vid0, vid1 }) != allDelaunayEdges.end());
      }
      maVtx2SkeletonSurface[v.first] = v.second;
    }
    else {
      // // pick one vertex order make it a polygon
      // std::vector<std::array<int, 2>> edgePairs;
      // edgePairs.reserve(v.second.size() * v.second.size());
      // for (int j = 0; j < (int)v.second.size(); j++) {
      //   for (int k = j + 1; k < (int)v.second.size(); k++) {
      //     int vid0 = v.second[j], vid1 = v.second[k];
      //     if (vid0 > vid1)
      //       std::swap(vid0, vid1);

      //     if (std::find(allDelaunayEdges.begin(), allDelaunayEdges.end(), std::array<int, 2>{ vid0, vid1 }) != allDelaunayEdges.end()) {
      //       edgePairs.emplace_back(std::array<int, 2>{vid0, vid1});
      //     }
      //   }
      // }
      // PGO_ALOG(edgePairs.size() >= v.second.size());

      // TODO make it a surface mesh by computing medial axis of the volume
      for (int j = 0; j < (int)v.second.size(); j++) {
        for (int k = j + 1; k < (int)v.second.size(); k++) {
          int vid0 = v.second[j], vid1 = v.second[k];
          if (vid0 > vid1)
            std::swap(vid0, vid1);
          PGO_ALOG(std::find(allDelaunayEdges.begin(), allDelaunayEdges.end(), std::array<int, 2>{ vid0, vid1 }) != allDelaunayEdges.end());
        }
      }
      maVtx2SkeletonSurface[v.first] = v.second;
    }
  }

  if constexpr (dumpMesh) {
    // Debug and check
    for (int i = 0; i < (int)allDelaunayEdges.size(); i++) {
      for (int j = i + 1; j < (int)allDelaunayEdges.size(); j++) {
        PGO_ALOG(allDelaunayEdges[i] != allDelaunayEdges[j]);
      }
    }

    pgo::Mesh::TriMeshGeo zz;
    for (int e = 0; e < (int)allDelaunayEdges.size(); e++) {
      auto &edge = allDelaunayEdges[e];
      int si = edge[0], sj = edge[1];
      zz.addMesh(pgo::Mesh::createSingleTriangleMesh(samplePoints[si], samplePoints[sj], samplePoints[si] + pgo::asVec3d(1e-6)));
    }
    zz.save("sampleNeighbors.obj");
  }

  // find triangles between every pair of connected components
  std::map<std::pair<int, int>, std::vector<int>> interfaceTriangles;

  for (const auto &edge : allDelaunayEdges) {
    int si = edge[0], sj = edge[1];

    for (int ti = 0; ti < (int)medialAxisSkeletonTrianglesAndEdges.size(); ti++) {
      if (faceSelFlags[ti] != si && faceSelFlags[ti] != sj)
        continue;

      auto it = triangleNeighbors.find(ti);
      PGO_ALOG(it != triangleNeighbors.end());
      for (const auto &ntri : it->second) {
        bool isInterfaceTriangle = (faceSelFlags[ti] == si && faceSelFlags[ntri] == sj) || (faceSelFlags[ti] == sj && faceSelFlags[ntri] == si);

        if (isInterfaceTriangle) {
          interfaceTriangles[std::pair(si, sj)].emplace_back(ti);
        }
      }
    }

    if constexpr (dumpMesh) {
      if (interfaceTriangles.find(std::pair(si, sj)) == interfaceTriangles.end())
        continue;
      pgo::Mesh::TriMeshGeo zz;
      for (int tri : interfaceTriangles[std::pair(si, sj)]) {
        if ((int)medialAxisSkeletonTrianglesAndEdges[tri].size() == 2) {
          zz.addMesh(pgo::Mesh::createSingleTriangleMesh(toVec3(medialAxisVertices[medialAxisSkeletonTrianglesAndEdges[tri][0]]),
            toVec3(medialAxisVertices[medialAxisSkeletonTrianglesAndEdges[tri][1]]),
            toVec3(medialAxisVertices[medialAxisSkeletonTrianglesAndEdges[tri][0]]) + pgo::asVec3d(1e-6)));
        }
        else {
          zz.addMesh(pgo::Mesh::createSingleTriangleMesh(toVec3(medialAxisVertices[medialAxisSkeletonTrianglesAndEdges[tri][0]]),
            toVec3(medialAxisVertices[medialAxisSkeletonTrianglesAndEdges[tri][1]]),
            toVec3(medialAxisVertices[medialAxisSkeletonTrianglesAndEdges[tri][2]])));
        }
      }
      zz.save(fmt::format("interfaceTriangles_{}_{}.obj", si, sj));
    }
  }

  // For each delauny edge, build "connected components" for all the interface triangles
  // connecting centers of triangles, if one edge is outside of the mesh, then it should be discarded
  std::map<std::pair<int, int>, std::vector<std::vector<int>>> interfaceTrianglesByCC;

  for (auto &pr : interfaceTriangles) {
    pgo::BasicAlgorithms::sortAndDeduplicateWithErase(pr.second);

    // std::vector<int> Q;
    std::set<int> visitedTri;
    std::vector<int> Q;

    while (1) {
      if (visitedTri.size() == pr.second.size())
        break;

      Q.clear();
      int selTriID = -1;
      for (int vi = 0; vi < (int)pr.second.size(); vi++) {
        if (visitedTri.find(pr.second[vi]) != visitedTri.end()) {
          continue;
        }

        visitedTri.emplace(pr.second[vi]);
        selTriID = pr.second[vi];
        Q.emplace_back(selTriID);
        break;
      }

      for (int vi = 0; vi < (int)pr.second.size(); vi++) {
        if (visitedTri.find(pr.second[vi]) != visitedTri.end()) {
          continue;
        }

        if (selTriID == pr.second[vi]) {
          continue;
        }

        // create line segment
        ES::V3d start, end, p0, p1, p2;
        if ((int)medialAxisSkeletonTrianglesAndEdges[selTriID].size() == 2) {
          p0 = toVec3(medialAxisVertices[medialAxisSkeletonTrianglesAndEdges[selTriID][0]]);
          p1 = toVec3(medialAxisVertices[medialAxisSkeletonTrianglesAndEdges[selTriID][1]]);
          start = (p0 + p1) / 2.0;
        }
        else {
          p0 = toVec3(medialAxisVertices[medialAxisSkeletonTrianglesAndEdges[selTriID][0]]);
          p1 = toVec3(medialAxisVertices[medialAxisSkeletonTrianglesAndEdges[selTriID][1]]);
          p2 = toVec3(medialAxisVertices[medialAxisSkeletonTrianglesAndEdges[selTriID][2]]);
          start = (p0 + p1 + p2) / 3.0;
        }

        if ((int)medialAxisSkeletonTrianglesAndEdges[pr.second[vi]].size() == 2) {
          p0 = toVec3(medialAxisVertices[medialAxisSkeletonTrianglesAndEdges[pr.second[vi]][0]]);
          p1 = toVec3(medialAxisVertices[medialAxisSkeletonTrianglesAndEdges[pr.second[vi]][1]]);
          end = (p0 + p1) / 2.0;
        }
        else {
          p0 = toVec3(medialAxisVertices[medialAxisSkeletonTrianglesAndEdges[pr.second[vi]][0]]);
          p1 = toVec3(medialAxisVertices[medialAxisSkeletonTrianglesAndEdges[pr.second[vi]][1]]);
          p2 = toVec3(medialAxisVertices[medialAxisSkeletonTrianglesAndEdges[pr.second[vi]][2]]);
          end = (p0 + p1 + p2) / 3.0;
        }

        // bool hasIntersection = targetMeshBVTree.hasLineSegmentIntersectionExact(targetMesh, start, end); // sometimes crushed w/ exact
        int hasIntersection = targetMeshBVTree.lineSegmentFirstIntersectionPoint(targetMesh, start, end);
        if (hasIntersection == -1) {
          Q.emplace_back(pr.second[vi]);
          visitedTri.emplace(pr.second[vi]);
        }
      }

      // If the newly added Q is connected with the existing Qs, then merge them
      std::vector<int> neighborQIDs;
      neighborQIDs.reserve(interfaceTrianglesByCC[pr.first].size());

      for (int qi : Q) {
        auto it = triangleNeighbors.find(qi);
        PGO_ALOG(it != triangleNeighbors.end());
        for (int prevQi = 0; std::vector<int> & prevQ : interfaceTrianglesByCC[pr.first]) {
          for (const int &prevQii : prevQ) {
            auto itQi = std::find(it->second.begin(), it->second.end(), prevQii);
            if (itQi != it->second.end()) {
              neighborQIDs.emplace_back(prevQi);
              break;
            }
          }
          prevQi++;
        }
      }

      if ((int)neighborQIDs.size() > 0) {
        pgo::BasicAlgorithms::sortAndDeduplicateWithErase(neighborQIDs);
        std::vector<int> newQ;
        newQ.insert(newQ.end(), Q.begin(), Q.end());
        for (int qi : neighborQIDs) {
          newQ.insert(newQ.end(), interfaceTrianglesByCC[pr.first][qi].begin(), interfaceTrianglesByCC[pr.first][qi].end());
        }
        Q = newQ;
        // remove neighborQIDs
        while (neighborQIDs.size() > 0) {
          interfaceTrianglesByCC[pr.first].erase(interfaceTrianglesByCC[pr.first].begin() + neighborQIDs.back());
          neighborQIDs.pop_back();
        }
      }

      interfaceTrianglesByCC[pr.first].emplace_back(Q);
    }
  }

  if constexpr (dumpMesh) {
    pgo::Mesh::TriMeshGeo zz;

    for (auto &pr : interfaceTrianglesByCC) {
      for (int pi = 0; pi < (int)pr.second.size(); pi++) {
        for (int j = 0; j < (int)pr.second[pi].size(); j++) {
          ES::V3d center;
          if ((int)medialAxisSkeletonTrianglesAndEdges[pr.second[pi][j]].size() == 2) {
            ES::V3d p0 = toVec3(medialAxisVertices[medialAxisSkeletonTrianglesAndEdges[pr.second[pi][j]][0]]);
            ES::V3d p1 = toVec3(medialAxisVertices[medialAxisSkeletonTrianglesAndEdges[pr.second[pi][j]][1]]);
            center = (p0 + p1) / 2;
          }
          else {
            ES::V3d p0 = toVec3(medialAxisVertices[medialAxisSkeletonTrianglesAndEdges[pr.second[pi][j]][0]]);
            ES::V3d p1 = toVec3(medialAxisVertices[medialAxisSkeletonTrianglesAndEdges[pr.second[pi][j]][1]]);
            ES::V3d p2 = toVec3(medialAxisVertices[medialAxisSkeletonTrianglesAndEdges[pr.second[pi][j]][2]]);
            center = (p0 + p1 + p2) / 3;
          }
          zz.addPos(center);
        }
      }
    }
    zz.save("cc.obj");

    // Debug store all the keys of interfaceTrianglesByCC
    std::vector<std::pair<int, int>> interfaceTrianglesByCCKeys;
    interfaceTrianglesByCCKeys.reserve(interfaceTrianglesByCC.size());
    for (const auto &pr : interfaceTrianglesByCC) {
      interfaceTrianglesByCCKeys.emplace_back(pr.first);
    }
  }

  HClockPt t9 = HClock::now();
  SPDLOG_LOGGER_INFO(pgo::Logging::lgr(), "cc analysis: {}s", duraSecond(t8, t9));

  finalSkeletonPoints.clear();
  finalSkeletonEdges.clear();
  finalSkeletonTriangles.clear();

  std::map<std::pair<int, int>, std::vector<std::vector<int>>> delaunayEdge2Polyline;

  pgo::Mesh::TriMeshGeo zzz, zzz1;

  for (auto &pr : interfaceTrianglesByCC) {
    int sampledPtID0 = pr.first.first, sampledPtID1 = pr.first.second;
    if (sampledPtID0 > sampledPtID1) {
      std::swap(sampledPtID0, sampledPtID1);
    }

    for (int j = 0; j < (int)pr.second.size(); j++) {
      double minPathDist = 1e100;
      int selTri = -1;

      for (int ti : pr.second[j]) {
        if (distances[sampledPtID0][ti] + distances[sampledPtID1][ti] < minPathDist) {
          minPathDist = distances[sampledPtID0][ti] + distances[sampledPtID1][ti];
          selTri = ti;
        }
      }

      PGO_ALOG(selTri >= 0);

      std::vector<int> path;
      path.reserve(medialAxisSkeletonTrianglesAndEdges.size());

      int curID = selTri;
      do {
        path.emplace_back(curID);
        curID = prevTriangleIDs[sampledPtID0][curID];
      } while (curID >= 0);

      std::reverse(path.begin(), path.end());

      curID = prevTriangleIDs[sampledPtID1][selTri];
      if (curID < 0) {
        // means the selTri is the same as the sampledPtID1, samplePtID1 is on interfaceTrianglesByCC
      }
      else {
        do {
          path.emplace_back(curID);
          curID = prevTriangleIDs[sampledPtID1][curID];
        } while (curID >= 0);
      }

      // check if any segment in the path is outside the target mesh, keep subdivide and project to MA // TODO project along the middle plane
      // for (auto it = path.begin(); it != path.end() - 1;) {
      //   while (1) {
      //     ES::V3d p0 = triangleCenters[*it], p1 = triangleCenters[*(it + 1)];
      //     bool hasIntersection = targetMeshBVTreeSmall.hasLineSegmentIntersectionExact(targetMeshSmall, p0, p1);
      //     if (!hasIntersection) {
      //       break;
      //     }
      //     else {
      //       ES::V3d midPt = (p0 + p1) * 0.5;
      //       double minDist = 1e20;
      //       int projTriID = -1;
      //       for (auto itTri = triangleCenters.begin(); itTri != triangleCenters.end(); itTri++) {
      //         if (len(*itTri - midPt) < minDist) {
      //           minDist = len(*itTri - midPt);
      //           projTriID = itTri - triangleCenters.begin();
      //         }
      //       }
      //       if (projTriID == *it || projTriID == *(it + 1)) {
      //         // try to project along the middle plane
      //         projTriID = -1;
      //         minDist = 1e20;
      //         ES::V3d normal = norm(p1 - p0);
      //         for (auto itTri = triangleCenters.begin(); itTri != triangleCenters.end(); itTri++) {
      //           double dist = dot(normal, (*itTri - midPt));
      //           if (dist < 1e-4) {
      //             if (len(*itTri - midPt) < minDist) {
      //               minDist = len(*itTri - midPt);
      //               projTriID = itTri - triangleCenters.begin();
      //             }
      //           }
      //         }
      //         PGO_ALOG(projTriID != *it && projTriID != *(it + 1));
      //       }
      //       path.insert(it + 1, projTriID);
      //     }
      //   }
      //   it++;
      // }

      std::list<std::pair<int, int>> edgeQ;
      edgeQ.emplace_back(0, (int)path.size() - 1);

      std::vector<double> pathAccDistances;
      pathAccDistances.emplace_back(0);

      for (int pi = 1; pi < (int)path.size(); pi++) {
        pathAccDistances.emplace_back((triangleCenters[path[pi]] - triangleCenters[path[pi - 1]]).norm() + pathAccDistances[pi - 1]);
      }

      pgo::Mesh::TriMeshGeo zzzz1;
      for (int pi = 0; pi < (int)path.size() - 1; pi++) {
        zzzz1.addMesh(pgo::Mesh::createSingleTriangleMesh(triangleCenters[path[pi]], triangleCenters[path[pi + 1]], triangleCenters[path[pi]] + pgo::asVec3d(1e-5)));
      }
      zzzz1.save("zzzz1.obj");

      for (auto it = edgeQ.begin(); it != edgeQ.end();) {
        ES::V3d p = triangleCenters[path[it->first]];
        ES::V3d q = triangleCenters[path[it->second]];

        bool hasIntersection = targetMeshBVTreeSmall.hasLineSegmentIntersectionExact(targetMeshSmall, p, q);

        if (hasIntersection) {
          int startID = it->first;
          int endID = it->second;

          double dist = pathAccDistances[endID] - pathAccDistances[startID];
          double halfDist = dist * 0.5;
          auto it1 = std::lower_bound(pathAccDistances.begin() + startID, pathAccDistances.begin() + endID + 1, halfDist + pathAccDistances[startID]);
          int midID;
          if (it1 != pathAccDistances.begin() + endID + 1) {
            midID = std::min((int)(it1 - pathAccDistances.begin()), endID - 1);
          }
          else {
            midID = endID - 1;
          }

          // fmt::print("startID {} midID {} endID {}\n", startID, midID, endID);
          // pgo::Mesh::TriMeshGeo zzzz;
          // for (auto seg : edgeQ) {
          //   zzzz.addMesh(pgo::Mesh::createSingleTriangleMesh(triangleCenters[path[seg.first]], triangleCenters[path[seg.second]], triangleCenters[path[seg.first]] + pgo::asVec3d(1e-5)));
          // }
          // zzzz.save("zzzz.obj");
          if (midID == startID || midID == endID) {
            // it means the primitives of the two centers are neighbors and are on the shortest path, and their connection is intersecting with the target mesh
            auto itStartID = triangleNeighbors.find(path[startID]);
            // pgo::Mesh::TriMeshGeo pts;
            // pts.addMesh(pgo::Mesh::createSingleTriangleMesh(triangleCenters[path[startID]], triangleCenters[path[endID]], triangleCenters[path[startID]] + pgo::asVec3d(1e-5)));
            // pts.save("debugpts.obj");
            PGO_ALOG(std::find(itStartID->second.begin(), itStartID->second.end(), path[endID]) != itStartID->second.end());

            // find the midpoint of the shared primitives of these two primitives
            int primitiveTypeStart = (int)medialAxisSkeletonTrianglesAndEdges[path[startID]].size();
            int primitiveTypeEnd = (int)medialAxisSkeletonTrianglesAndEdges[path[endID]].size();
            PGO_ALOG(primitiveTypeStart == 2 || primitiveTypeStart == 3);
            PGO_ALOG(primitiveTypeEnd == 2 || primitiveTypeEnd == 3);

            std::vector<int> startPrimitiveIDs = medialAxisSkeletonTrianglesAndEdges[path[startID]];
            std::vector<int> endPrimitiveIDs = medialAxisSkeletonTrianglesAndEdges[path[endID]];
            std::vector<int> sharedPrimitiveIDs;
            std::set_intersection(startPrimitiveIDs.begin(), startPrimitiveIDs.end(), endPrimitiveIDs.begin(), endPrimitiveIDs.end(), std::back_inserter(sharedPrimitiveIDs));
            PGO_ALOG((int)sharedPrimitiveIDs.size() == 1 || (int)sharedPrimitiveIDs.size() == 2);

            ES::V3d newAddingPt;
            EK::Point_3 newAddingPtEK;
            if ((int)sharedPrimitiveIDs.size() == 1) {
              newAddingPtEK = medialAxisVertices[sharedPrimitiveIDs[0]];
              newAddingPt = toVec3(newAddingPtEK);
            }
            else {
              EK::Point_3 newAddingPtEK0 = medialAxisVertices[sharedPrimitiveIDs[0]];
              EK::Point_3 newAddingPtE1 = medialAxisVertices[sharedPrimitiveIDs[1]];
              newAddingPtEK = CGAL::midpoint(newAddingPtEK0, newAddingPtE1);
              newAddingPt = toVec3(newAddingPtEK);
            }

            // Add a triangle center to triangleCenters and triangleCentersEK
            triangleCenters.emplace_back(newAddingPt);
            triangleCentersEK.emplace_back(newAddingPtEK);

            auto itPathBegin = path.begin();
            path.insert(itPathBegin + endID, (int)triangleCenters.size() - 1);

            // update edgeQ
            for (auto itEdgeQ = edgeQ.begin(); itEdgeQ != edgeQ.end(); ++itEdgeQ) {
              if (itEdgeQ->first >= endID) {
                itEdgeQ->first++;
              }
              if (itEdgeQ->second >= endID) {
                itEdgeQ->second++;
              }
            }

            edgeQ.emplace_back(startID, startID + 1);
            edgeQ.emplace_back(startID + 1, endID + 1);
            auto it2 = std::next(it);
            auto it3 = std::next(it2);
            edgeQ.erase(it);
            it = it3;
            continue;
          }

          PGO_ALOG(midID > startID);
          PGO_ALOG(endID > midID);

          edgeQ.emplace_back(startID, midID);
          edgeQ.emplace_back(midID, endID);

          auto it2 = std::next(it);
          edgeQ.erase(it);
          it = it2;
        }
        else {
          ++it;
        }
      }

      if (dumpMesh) {
        for (int pi = 0; pi < (int)path.size() - 1; pi++) {
          zzz.addMesh(pgo::Mesh::createSingleTriangleMesh(triangleCenters[path[pi]], triangleCenters[path[pi + 1]], triangleCenters[path[pi]] + pgo::asVec3d(1e-5)));
        }

        for (auto seg : edgeQ) {
          zzz1.addMesh(pgo::Mesh::createSingleTriangleMesh(triangleCenters[path[seg.first]], triangleCenters[path[seg.second]], triangleCenters[path[seg.first]] + pgo::asVec3d(1e-5)));
        }
      }

      std::vector<int> edgeQIDs;
      for (int id = 0; auto seg : edgeQ) {
        finalSkeletonEdges.emplace_back(path[seg.first], path[seg.second]);
        if (path[seg.first] == -1 || path[seg.second] == -1) {
          std::cout << "debug" << std::endl;
        }
        edgeQIDs.emplace_back(seg.first);
        edgeQIDs.emplace_back(seg.second);
        id++;
      }

      pgo::BasicAlgorithms::sortAndDeduplicateWithErase(edgeQIDs);
      // PGO_ALOG(delaunayEdge2Polyline.find(pr.first) == delaunayEdge2Polyline.end()); // the delaunay edge should not be in the map yet
      std::vector<int> pathIDs(edgeQIDs.size());
      for (int i = 0; i < (int)edgeQIDs.size(); i++) {
        pathIDs[i] = path[edgeQIDs[i]];
      }

      auto it = delaunayEdge2Polyline.find(pr.first);
      if (it != delaunayEdge2Polyline.end()) {
        // PGO_ALOG((int)it->second.size() == j);

        // add a safe (filtering out some bug cases if any)
        bool isExisting = false;
        for (auto &ply : it->second) {
          std::vector<int> plyNew = pathIDs;
          pgo::BasicAlgorithms::sortAndDeduplicateWithErase(plyNew);
          std::vector<int> plyOld = ply;
          pgo::BasicAlgorithms::sortAndDeduplicateWithErase(plyOld);
          if (plyOld == plyNew) {
            isExisting = true;
          }
        }

        if (isExisting == false) {
          it->second.emplace_back(pathIDs);
        }
      }
      else {
        auto ret = delaunayEdge2Polyline.emplace_hint(it, pr.first, std::vector<std::vector<int>>());
        ret->second.emplace_back(pathIDs);
      }
    }
  }

  if (dumpMesh) {
    zzz.save("rr.obj");
    zzz1.save("rr1.obj");
  }

  std::map<int, int> usedVertices;
  for (const auto &e : finalSkeletonEdges) {
    usedVertices.emplace(e.first, 0);
    usedVertices.emplace(e.second, 0);
  }

  std::vector<EK::Point_3> finalSkeletonPtsEK;
  finalSkeletonPtsEK.reserve(usedVertices.size());
  for (auto &pr : usedVertices) {
    finalSkeletonPoints.emplace_back(triangleCenters[pr.first]);
    pr.second = (int)finalSkeletonPoints.size() - 1;
    finalSkeletonPtsEK.emplace_back(triangleCentersEK[pr.first]);
  }

  for (auto &e : finalSkeletonEdges) {
    auto it = usedVertices.find(e.first);
    PGO_ALOG(it != usedVertices.end());

    e.first = it->second;

    it = usedVertices.find(e.second);
    PGO_ALOG(it != usedVertices.end());

    e.second = it->second;
  }

  // check the keys of interfaceTrianglesByCC are always in increasing order
  for (const auto &pr : interfaceTrianglesByCC) {
    PGO_ALOG(pr.first.first < pr.first.second);
  }

  if constexpr (dumpMesh) {
    // Debug
    // show all tetCenters
    pgo::Mesh::TriMeshGeo tetCentersMesh;
    for (const auto &tetCenter : tetCenters) {
      tetCentersMesh.addPos(toVec3(tetCenter.second));
    }
    tetCentersMesh.save("tetCenters.obj");
  }

  // Add triangles
  if constexpr (1) {
    for (const auto &maVtx2SklSurf : maVtx2SkeletonSurface) {
      // add triangle
      // TODO check if is outside and postprocessing
      std::vector<int> polylineFinalMedialAxisVtxIDs;

      std::vector<std::array<int, 3>> allAddingTriangles;
      // create triangle // for any triplet of end points
      if ((int)maVtx2SklSurf.second.size() == 3) {
        allAddingTriangles.emplace_back(std::array<int, 3>{ maVtx2SklSurf.second[0], maVtx2SklSurf.second[1], maVtx2SklSurf.second[2] });
      }
      else {
        PGO_ALOG((int)maVtx2SklSurf.second.size() > 3);
        // select one point, the rest of them form a polygon, pick every edge of the polygon with the selected point and form a triangle
        int selTopID = 0;
        for (int ti = 1; ti < (int)maVtx2SklSurf.second.size(); ti++) {
          int tj = ti + 1;
          if (tj > (int)maVtx2SklSurf.second.size() - 1) {
            tj = 1;
          }
          allAddingTriangles.emplace_back(std::array<int, 3>{ maVtx2SklSurf.second[selTopID], maVtx2SklSurf.second[ti], maVtx2SklSurf.second[tj] });
        }
      }

      // for (int ti = 0; ti < (int)maVtx2SklSurf.second.size(); ti++) {
      //   for (int tj = ti + 1; tj < (int)maVtx2SklSurf.second.size(); tj++) {
      //     for (int tk = tj + 1; tk < (int)maVtx2SklSurf.second.size(); tk++) {
      //       allAddingTriangles.emplace_back(std::array<int, 3>{ maVtx2SklSurf.second[ti], maVtx2SklSurf.second[tj], maVtx2SklSurf.second[tk] });
      //     }
      //   }
      // }

      for (const auto &addTri : allAddingTriangles) {
        // int numEnds = maVtx2SklSurf.second.size();
        int numEnds = 3;
        std::vector<std::vector<int>> polylines(numEnds);
        for (int i = 0; i < numEnds; i++) {
          int vtxID0 = addTri[i], vtxID1 = addTri[(i + 1) % numEnds];
          if (vtxID0 > vtxID1) {
            std::swap(vtxID0, vtxID1);
          }

          const std::vector<std::vector<int>> &polyline_list = delaunayEdge2Polyline[std::make_pair(vtxID0, vtxID1)];
          PGO_ALOG((int)polyline_list.size() > 0);

          if ((int)polyline_list.size() == 1) {
            polylines[i] = polyline_list[0];
          }
          else {
            // determine which one of polyline_list to use, check the line segs connecting points on the path with triangle center are inside the mesh
            // TODO change to find the connected interface triangle
            int selPlyID = -1;
            const std::vector<std::vector<int>> &interfaceTris = interfaceTrianglesByCC[std::make_pair(vtxID0, vtxID1)];
            PGO_ALOG(interfaceTris.size() == polyline_list.size());
            for (int plyID = 0; plyID < (int)interfaceTris.size(); plyID++) {
              const std::vector<int> &inTri = interfaceTris[plyID];
              auto neighborTris = vtxTriangleNeighbors[maVtx2SklSurf.first];
              bool isNeighbor = false;
              for (int nTri : neighborTris) {
                if (std::find(inTri.begin(), inTri.end(), nTri) != inTri.end()) {
                  isNeighbor = true;
                  break;
                }
              }
              if (isNeighbor) {
                selPlyID = plyID;
                break;
              }
            }
            // ES::V3d start = toVec3(medialAxisVertices[maVtx2SklSurf.first]);
            // for (int plyID = 0; auto &ply : polyline_list) {
            //   bool noIntersection = true;
            //   for (int id = 0; id < (int)ply.size(); id++) {
            //     ES::V3d end = triangleCenters[ply[id]];

            //     // Debug
            //     if constexpr (dumpMesh) {
            //       pgo::Mesh::TriMeshGeo polylineDebugMesh;
            //       polylineDebugMesh.addPos(start);
            //       polylineDebugMesh.addPos(end);
            //       polylineDebugMesh.addPos(start + pgo::asVec3d(1e-5));
            //       polylineDebugMesh.addTri(ES::V3i(0, 1, 2));
            //       polylineDebugMesh.save(fmt::format("polylineDebugMesh{}.obj", plyID));
            //     }

            //     bool hasIntersection = targetMeshBVTree.hasLineSegmentIntersectionExact(targetMesh, start, end);

            //     if (hasIntersection) {
            //       noIntersection = false;
            //       break;
            //     }
            //   }

            //   if (noIntersection) {
            //     selPlyID = plyID;
            //     break;  // maybe comment and check if selPlyID is unique
            //   }

            //   plyID++;
            // }

            PGO_ALOG(selPlyID >= 0);
            polylines[i] = polyline_list[selPlyID];
          }

          if constexpr (dumpMesh) {
            pgo::Mesh::TriMeshGeo polylineDebugMesh;
            for (int id = 0; id < (int)polylines[i].size() - 1; id++) {
              polylineDebugMesh.addPos(triangleCenters[polylines[i][id]]);
              polylineDebugMesh.addPos(triangleCenters[polylines[i][(id + 1)]]);
              polylineDebugMesh.addPos(triangleCenters[polylines[i][id]] + pgo::asVec3d(1e-5));
              polylineDebugMesh.addTri(ES::V3i(id * 3, id * 3 + 1, id * 3 + 2));
            }

            polylineDebugMesh.save(fmt::format("polylineDebugMesh{}.obj", i));
          }
        }

        polylineFinalMedialAxisVtxIDs.insert(polylineFinalMedialAxisVtxIDs.end(), polylines[0].begin(), polylines[0].end());
        int endID;
        for (int i = 0; i < numEnds - 1; i++) {
          int lastEndID = polylineFinalMedialAxisVtxIDs.back();
          int nextPolylineID = (i + 1) % numEnds;
          PGO_ALOG(lastEndID == polylines[nextPolylineID].front() || lastEndID == polylines[nextPolylineID].back());
          int excludeEnd = 0;
          if (i == numEnds - 2) {
            excludeEnd = 1;
          }
          if (lastEndID == polylines[nextPolylineID].back()) {
            polylineFinalMedialAxisVtxIDs.insert(polylineFinalMedialAxisVtxIDs.end(), polylines[nextPolylineID].rbegin() + 1, polylines[nextPolylineID].rend() - excludeEnd);
            endID = polylines[nextPolylineID].front();
          }
          else {
            polylineFinalMedialAxisVtxIDs.insert(polylineFinalMedialAxisVtxIDs.end(), polylines[nextPolylineID].begin() + 1, polylines[nextPolylineID].end() - excludeEnd);
            endID = polylines[nextPolylineID].back();
          }
        }
        PGO_ALOG(endID == polylines[0].front());

        // triangulate polyline
        std::vector<EK::Point_3> polylineFinalMedialAxisVtxs;
        polylineFinalMedialAxisVtxs.reserve(polylineFinalMedialAxisVtxIDs.size());
        for (int i = 0; i < (int)polylineFinalMedialAxisVtxIDs.size(); i++) {
          EK::Point_3 p = triangleCentersEK[polylineFinalMedialAxisVtxIDs[i]];
          polylineFinalMedialAxisVtxs.emplace_back(std::move(p));
        }

        // for (int i = 0; i < (int)polylineFinalMedialAxisVtxs.size(); i++) {
        //   pgo::Mesh::TriMeshGeo debugMesh;
        //   debugMesh.addPos(toVec3(polylineFinalMedialAxisVtxs[i]));
        //   debugMesh.save(fmt::format("debugMeshPts{}.obj", i));
        // }

        // pgo::Mesh::TriMeshGeo debugMesh;
        // for (int i = 0; i < (int)polylineFinalMedialAxisVtxs.size(); i++) {
        //   if (i == 0) {
        //     debugMesh.addPos(toVec3(polylineFinalMedialAxisVtxs[i]));
        //   } else {
        //     debugMesh.addPos(toVec3(polylineFinalMedialAxisVtxs[i]));
        //     debugMesh.addPos(toVec3(polylineFinalMedialAxisVtxs[i]) + pgo::asVec3d(1e-5));
        //     debugMesh.addTri(ES::V3i((i - 1) * 2, (i - 1) * 2 + 1, (i - 1) * 2 + 2));
        //   }
        //   debugMesh.save(fmt::format("debugMesh{}.obj", i));
        // }
        std::vector<CGAL::Triple<int, int, int>> triangles;
        if ((int)polylineFinalMedialAxisVtxs.size() == 3) {
          triangles.emplace_back(0, 1, 2);
        }
        else if ((int)polylineFinalMedialAxisVtxs.size() > 3) {
          // CGAL::Polygon_mesh_processing::triangulate_hole_polyline(polylineFinalMedialAxisVtxs, std::back_inserter(triangles));
          triangulatePolygon(polylineFinalMedialAxisVtxs, targetMesh, targetMeshBVTree, triangles);
        }

        // pgo::Mesh::TriMeshGeo debugTriMesh;
        // for (int t = 0; t < (int)triangles.size(); t++) {
        //   debugTriMesh.addPos(toVec3(polylineFinalMedialAxisVtxs[triangles[t].get<0>()]));
        //   debugTriMesh.addPos(toVec3(polylineFinalMedialAxisVtxs[triangles[t].get<1>()]));
        //   debugTriMesh.addPos(toVec3(polylineFinalMedialAxisVtxs[triangles[t].get<2>()]));
        //   debugTriMesh.addTri(ES::V3i(t * 3, t * 3 + 1, t * 3 + 2));
        // }
        // debugTriMesh.save("debugTriMesh.obj");

        // debug triangulation
        for (const auto &tri : triangles) {
          if (polylineFinalMedialAxisVtxIDs[tri.get<0>()] == polylineFinalMedialAxisVtxIDs[tri.get<1>()]) {
            continue;
          }
          if (polylineFinalMedialAxisVtxIDs[tri.get<0>()] == polylineFinalMedialAxisVtxIDs[tri.get<2>()]) {
            continue;
          }
          if (polylineFinalMedialAxisVtxIDs[tri.get<1>()] == polylineFinalMedialAxisVtxIDs[tri.get<2>()]) {
            continue;
          }
          finalSkeletonTriangles.emplace_back(usedVertices[polylineFinalMedialAxisVtxIDs[tri.get<0>()]],
            usedVertices[polylineFinalMedialAxisVtxIDs[tri.get<1>()]],
            usedVertices[polylineFinalMedialAxisVtxIDs[tri.get<2>()]]);
        }
      }
    }
  }

  // post process finalSkeletonEdges, if it exists in finalSkeletonTriangles, remove it
  fmt::print("before remove finalSkeletonEdges.size() = {}\n", finalSkeletonEdges.size());
  for (const auto &tri : finalSkeletonTriangles) {
    std::array<int, 3> vIdx = { std::get<0>(tri), std::get<1>(tri), std::get<2>(tri) };
    for (int i = 0; i < 3; i++) {
      auto itErase = std::find(finalSkeletonEdges.begin(), finalSkeletonEdges.end(), std::pair<int, int>{ vIdx[i], vIdx[(i + 1) % 3] });
      if (itErase != finalSkeletonEdges.end()) {
        finalSkeletonEdges.erase(itErase);
      }
      itErase = std::find(finalSkeletonEdges.begin(), finalSkeletonEdges.end(), std::pair<int, int>{ vIdx[(i + 1) % 3], vIdx[i] });
      if (itErase != finalSkeletonEdges.end()) {
        finalSkeletonEdges.erase(itErase);
      }
    }
  }
  fmt::print("after remove finalSkeletonEdges.size() = {}\n", finalSkeletonEdges.size());

  pgo::Mesh::TriMeshGeo finalSkeleton;
  for (int i = 0; const auto &e : finalSkeletonEdges) {
    finalSkeleton.addPos(finalSkeletonPoints[e.first]);
    finalSkeleton.addPos(finalSkeletonPoints[e.second]);
    finalSkeleton.addPos(finalSkeletonPoints[e.second] + pgo::asVec3d(1e-6));
    finalSkeleton.addTri(ES::V3i(3 * i, 3 * i + 1, 3 * i + 2));
    i++;
  }
  finalSkeleton.save("rrr.obj");

  pgo::Mesh::TriMeshGeo finalSkeletonTri;
  for (int t = 0; const auto &tri : finalSkeletonTriangles) {
    finalSkeletonTri.addPos(finalSkeletonPoints[std::get<0>(tri)]);
    finalSkeletonTri.addPos(finalSkeletonPoints[std::get<1>(tri)]);
    finalSkeletonTri.addPos(finalSkeletonPoints[std::get<2>(tri)]);
    finalSkeletonTri.addTri(ES::V3i(3 * t, 3 * t + 1, 3 * t + 2));
    t++;
  }
  finalSkeletonTri.save("rrrTri.obj");

  // save simplified medial axis

  if (saveExactResults.length() > 0) {
    if (std::filesystem::exists(saveExactResults) == false) {
      std::filesystem::create_directories(saveExactResults);
    }

    SkeletonExact maSimplified;
    maSimplified.initFromInput(finalSkeletonPtsEK, finalSkeletonEdges, finalSkeletonTriangles);

    maSimplified.saveDisplayMesh(fmt::format("{}/maSimplified.obj", saveExactResults).c_str());
    maSimplified.save(fmt::format("{}/maInitSimplified.json", saveExactResults).c_str());
  }

  // // Debug load and save
  // MedialAxisSkeletonExact maSimplifiedLoad;
  // maSimplifiedLoad.load(fmt::format("{}/maInitSimplified.json", savePath).c_str());
  // maSimplifiedLoad.saveDisplayMesh(fmt::format("{}/maSimplifiedLoad.obj", savePath).c_str());
}

void MedialAxisRepresentation::triangulatePolygon(const std::vector<EK::Point_3> &polylineVtxs,
  const pgo::Mesh::TriMeshGeo &targetMesh,
  const pgo::Mesh::TriMeshBVTree &targetMeshBVTree,
  std::vector<CGAL::Triple<int, int, int>> &triangles)
{
  PGO_ALOG((int)polylineVtxs.size() > 3);

  std::vector<ES::V3d> polylineVtxsV3d(polylineVtxs.size());
  for (int i = 0; i < (int)polylineVtxs.size(); i++) {
    polylineVtxsV3d[i] = toVec3(polylineVtxs[i]);
  }

  std::vector<int> polylineVtxIDs(polylineVtxs.size());
  std::iota(polylineVtxIDs.begin(), polylineVtxIDs.end(), 0);

  int numTriLastIter = -1;
  while (true) {
    // find one ear vertex
    int earID = -1, prevID = -1, nextID = -1;
    for (int i = 0; i < (int)polylineVtxIDs.size(); i++) {
      prevID = (i - 1 + polylineVtxIDs.size()) % polylineVtxIDs.size();
      nextID = (i + 1) % polylineVtxIDs.size();

      // check it is a convex vertex
      ES::V3d v0 = polylineVtxsV3d[polylineVtxIDs[prevID]] - polylineVtxsV3d[polylineVtxIDs[i]];
      ES::V3d v1 = polylineVtxsV3d[polylineVtxIDs[nextID]] - polylineVtxsV3d[polylineVtxIDs[i]];
      // angle of two vec
      double angle = std::acos(v0.dot(v1) / ((v0).norm() * (v1).norm()));
      if (angle > M_PI) {
        continue;
      }

      // check if it is an ear
      if (polylineVtxs[polylineVtxIDs[prevID]] == polylineVtxs[polylineVtxIDs[nextID]]) {
        earID = i;
        break;
      }
      ES::V3d start = polylineVtxsV3d[polylineVtxIDs[prevID]];
      ES::V3d end = polylineVtxsV3d[polylineVtxIDs[nextID]];
      bool hasIntersection = targetMeshBVTree.hasLineSegmentIntersectionExact(targetMesh, start, end);
      if (!hasIntersection) {
        earID = i;
        break;
      }
    }

    if (earID == -1) {
      break;
    }

    if ((int)polylineVtxIDs.size() == 2) {
      break;
    }

    PGO_ALOG(prevID != earID);
    PGO_ALOG(nextID != earID);
    PGO_ALOG(prevID != nextID);
    PGO_ALOG(polylineVtxIDs[prevID] != polylineVtxIDs[nextID] && polylineVtxIDs[prevID] != polylineVtxIDs[earID] && polylineVtxIDs[nextID] != polylineVtxIDs[earID]);
    // add triangle
    triangles.emplace_back(polylineVtxIDs[prevID], polylineVtxIDs[earID], polylineVtxIDs[nextID]);

    // update polylineVtxIDs
    polylineVtxIDs.erase(polylineVtxIDs.begin() + earID);
  }
}

void MedialAxisRepresentation::preprocessingMA(
  const std::vector<EK::Point_3> &medialAxisVertices, const std::vector<std::vector<int>> &medialAxisFacets,
  const pgo::Mesh::TriMeshGeo &mesh, const pgo::Mesh::TriMeshBVTree &meshBVTree,
  std::vector<EK::Point_3> &finalMedialAxisVertices, std::vector<std::vector<int>> &finalMedialAxisFacets)
{
  finalMedialAxisFacets.reserve(medialAxisFacets.size());

  for (const auto &f : medialAxisFacets) {
    if ((int)f.size() == 2) {
      ES::V3d p = toVec3(medialAxisVertices[f[0]]);
      ES::V3d q = toVec3(medialAxisVertices[f[1]]);
      bool hasIntersection = meshBVTree.hasLineSegmentIntersectionExact(mesh, p, q);
      if (!hasIntersection)
        finalMedialAxisFacets.emplace_back(std::vector<int>{ f[0], f[1] });
    }
    else if ((int)f.size() == 3) {
      bool hasIntersection = false;
      for (int j = 0; j < 3; j++) {
        ES::V3d p = toVec3(medialAxisVertices[f[j]]);
        ES::V3d q = toVec3(medialAxisVertices[f[(j + 1) % 3]]);
        hasIntersection = meshBVTree.hasLineSegmentIntersectionExact(mesh, p, q);
        if (hasIntersection)
          break;
      }
      if (!hasIntersection)
        finalMedialAxisFacets.emplace_back(f);
    }
    else {
      PGO_ALOG(false);
    }
  }

  // raw MA clean up, remove all necessary disconnected components
  std::unordered_map<int, std::vector<int>> vtxTriNeighbor;
  for (int ti = 0; ti < (int)finalMedialAxisFacets.size(); ti++) {
    for (int j = 0; j < (int)finalMedialAxisFacets[ti].size(); j++) {
      auto it = vtxTriNeighbor.find(finalMedialAxisFacets[ti][j]);
      if (it != vtxTriNeighbor.end()) {
        it->second.emplace_back(ti);
      }
      else {
        vtxTriNeighbor.emplace(finalMedialAxisFacets[ti][j], std::vector<int>{ ti });
      }
    }
  }

  std::unordered_map<int, std::vector<int>> faceNeighbors;
  for (const auto &pr : vtxTriNeighbor) {
    for (int i = 0; i < (int)pr.second.size(); i++) {
      for (int j = i + 1; j < (int)pr.second.size(); j++) {
        faceNeighbors[pr.second[i]].emplace_back(pr.second[j]);
        faceNeighbors[pr.second[j]].emplace_back(pr.second[i]);
      }
    }
  }
  for (auto &pr : faceNeighbors) {
    pgo::BasicAlgorithms::sortAndDeduplicateWithErase(pr.second);
  }

  std::cout << "finalMedialAxisFacets.size(): " << finalMedialAxisFacets.size() << std::endl;

  std::vector<int> faceVisited(finalMedialAxisFacets.size(), 0);
  std::vector<int> toRemoveFaces;
  std::vector<int> faceQueue;
  int start = 0;
  while (1) {
    start = 0;
    faceQueue.clear();
    int unvisitedFace = -1;
    for (int i = 0; i < (int)faceVisited.size(); i++) {
      if (faceVisited[i] == 0) {
        unvisitedFace = i;
        break;
      }
    }
    PGO_ALOG(unvisitedFace >= 0);

    std::memset(faceVisited.data(), 0, sizeof(int) * faceVisited.size());
    faceQueue.emplace_back(unvisitedFace);
    faceVisited[unvisitedFace] = 1;

    // BFS
    while ((int)faceQueue.size() > start) {
      int curFace = faceQueue[start++];

      auto it = faceNeighbors.find(curFace);
      // ALOG(it != faceNeighbors.end());
      if (it == faceNeighbors.end())
        continue;

      for (const auto &fi : it->second) {
        if (faceVisited[fi] == 0) {
          faceQueue.emplace_back(fi);
          faceVisited[fi] = 1;
        }
      }
    }

    if ((int)faceQueue.size() > (int)finalMedialAxisFacets.size() / 2) {
      for (int i = 0; i < (int)faceVisited.size(); i++) {
        if (faceVisited[i] == 0) {
          toRemoveFaces.emplace_back(i);
        }
      }
      break;
    }
  }

  // Debug
  for (int i = 0; i < (int)finalMedialAxisFacets.size(); i++) {
    if (faceNeighbors.find(i) == faceNeighbors.end()) {
      toRemoveFaces.emplace_back(i);
    }
  }

  pgo::BasicAlgorithms::sortAndDeduplicateWithErase(toRemoveFaces);
  std::cout << "toRemoveFaces.size(): " << toRemoveFaces.size() << std::endl;

  for (int fi = (int)toRemoveFaces.size() - 1; fi >= 0; fi--) {
    finalMedialAxisFacets.erase(finalMedialAxisFacets.begin() + toRemoveFaces[fi]);
  }
  std::cout << "finalMedialAxisFacets.size(): " << finalMedialAxisFacets.size() << std::endl;

  finalMedialAxisVertices = medialAxisVertices;
  std::vector<int> toRemoveVtx;
  for (int vi = 0; vi < (int)finalMedialAxisVertices.size(); vi++) {
    bool found = false;
    for (int fi = 0; fi < (int)finalMedialAxisFacets.size(); fi++) {
      if (std::find(finalMedialAxisFacets[fi].begin(), finalMedialAxisFacets[fi].end(), vi) != finalMedialAxisFacets[fi].end()) {
        found = true;
        break;
      }
    }
    if (!found) {
      toRemoveVtx.emplace_back(vi);
    }
  }

  // remove duplicated vertices
  std::vector<std::pair<int, int>> duplicatedVtxPairs;
  for (int vi = 0; vi < (int)finalMedialAxisVertices.size(); vi++) {
    for (int vj = vi + 1; vj < (int)finalMedialAxisVertices.size(); vj++) {
      if (finalMedialAxisVertices[vi] == finalMedialAxisVertices[vj]) {
        duplicatedVtxPairs.emplace_back(vi, vj);
      }
    }
  }

  for (int fi = 0; fi < (int)finalMedialAxisFacets.size(); fi++) {
    for (int vi = 0; vi < (int)finalMedialAxisFacets[fi].size(); vi++) {
      for (const auto &pr : duplicatedVtxPairs) {
        if (finalMedialAxisFacets[fi][vi] == pr.second) {
          finalMedialAxisFacets[fi][vi] = pr.first;
          std::set<int> tmpSet(finalMedialAxisFacets[fi].begin(), finalMedialAxisFacets[fi].end());
          PGO_ALOG((int)tmpSet.size() == (int)finalMedialAxisFacets[fi].size());  // no duplicated ids
        }
      }
    }
  }

  for (const auto &pr : duplicatedVtxPairs) {
    toRemoveVtx.emplace_back(pr.second);
  }

  pgo::BasicAlgorithms::sortAndDeduplicateWithErase(toRemoveVtx);
  std::cout << "toRemoveVtx.size(): " << toRemoveVtx.size() << std::endl;

  for (int vi = (int)toRemoveVtx.size() - 1; vi >= 0; vi--) {
    finalMedialAxisVertices.erase(finalMedialAxisVertices.begin() + toRemoveVtx[vi]);
  }

  // change finalMedialAxisFacets
  for (int fi = 0; fi < (int)finalMedialAxisFacets.size(); fi++) {
    for (int vi = 0; vi < (int)finalMedialAxisFacets[fi].size(); vi++) {
      int oldVtx = finalMedialAxisFacets[fi][vi];
      int newVtx = std::lower_bound(toRemoveVtx.begin(), toRemoveVtx.end(), oldVtx) - toRemoveVtx.begin();
      finalMedialAxisFacets[fi][vi] = oldVtx - newVtx;
    }
    pgo::BasicAlgorithms::sortAndDeduplicateWithErase(finalMedialAxisFacets[fi]);
  }
}