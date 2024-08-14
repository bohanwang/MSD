
#include "pgoLogging.h"
#include "initPredicates.h"
#include "geogramInterface.h"

#include <iostream>

int main(int argc, char *argv[])
{
  using namespace MedialAxisRepresentation;

  if (argc != 5) {
    std::cerr << argv[0] << " <mesh> <#pt> <#iter> <#+pt>" << std::endl;
    return 1;
  }

  std::ofstream("cmd.txt") << argv[0] << ' ' << argv[1] << ' ' << argv[2] << ' ' << argv[3] << ' ' << argv[4];

  pgo::Mesh::initPredicates();
  pgo::

  GeogramUtilities::initGEO();
  Logger::init();
  VegaFEM::Logging::init();

  std::cout << "Filename: " << argv[1] << '\n'
            << "#tgt points: " << argv[2] << '\n'
            << "#iter: " << argv[3] << '\n'
            << "#added pt: " << argv[4] << std::endl;

  int numTargetPoints = std::stoi(argv[2]);
  int numIter = std::stoi(argv[3]);
  int numAddPtsPerIter = std::stoi(argv[4]);

  std::set_terminate([]() {
    std::cout << "Unhandled exception\n"
              << std::flush;
    std::abort();
  });

  TriMeshGeo inputMesh;
  std::string filename = argv[1];
  int ret = 1;
  if (!filename.compare(filename.size() - 4, 4, ".obj")) {
    ret = LibIGLInterface::loadObjMesh(argv[1], inputMesh);
  }
  else if (!filename.compare(filename.size() - 4, 4, ".off")) {
    ret = LibIGLInterface::loadOffMesh(argv[1], inputMesh);
  }

  if (ret != 0)
    return 1;

  // BoundingBox bb(inputMesh.positions());
  // double l = bb.longestSide().second;

  // for (int vi = 0; vi < inputMesh.numVertices(); vi++) {
  //   Vec3d diff = inputMesh.pos(vi) - bb.bmin();
  //   diff /= l;

  //   inputMesh.pos(vi) = diff;
  // }
  inputMesh.save("m.obj");

  TriMeshBVTree inputMeshBVTree;
  inputMeshBVTree.buildByInertiaPartition(inputMesh);

  TriMeshPseudoNormal inputMeshNormal;
  inputMeshNormal.buildPseudoNormals(inputMesh);

#if 0
    std::vector<EK::Point_3> vtx;
    std::vector<std::vector<int>> facets;
    MedialAxisRepresentation::computeMedialAxis(inputMesh, vtx, facets);

    std::vector<std::pair<int, int>> edges;
    std::vector<std::tuple<int, int, int>> triangles;

    std::vector<EK::Point_3> polylines;
    std::vector<CGAL::Triple<int, int, int>> trianglesTemp;
    for (const auto &polygon : facets) {
      polylines.clear();

      if (polygon.size() == 2ull) {
        edges.emplace_back(polygon[0], polygon[1]);
        if (edges.back().first > edges.back().second) {
          std::swap(edges.back().first, edges.back().second);
        }
      }
      else if (polygon.size() == 3ull) {
        triangles.emplace_back(polygon[0], polygon[1], polygon[2]);
      }
      else if (polygon.size() > 3ull) {
        for (int vi : polygon)
          polylines.emplace_back(vtx[vi]);

        trianglesTemp.clear();
        CGAL::Polygon_mesh_processing::triangulate_hole_polyline(polylines, std::back_inserter(trianglesTemp));
        ALOG(trianglesTemp.size() > 0ull);

        for (const auto &tri : trianglesTemp) {
          std::array<int, 3> newFaces{ polygon[tri.first], polygon[tri.second], polygon[tri.third] };
          triangles.emplace_back(newFaces[0], newFaces[1], newFaces[2]);
        }
      }
    }

    MedialAxisSkeletonExact ma;
    ma.initFromInput(vtx, edges, triangles);
    std::cout << "is one cc: " << ma.isConnected() << std::endl;
    std::cout << "is consistent: " << ma.isConsistent() << std::endl;

    ma.saveDisplayMesh("ma1.obj");
#endif

  std::vector<ES::V3d> skeletonVtx;
  std::vector<std::pair<int, int>> skeletonEdges;
  std::vector<std::tuple<int, int, int>> skeletonTriangles;

  MedialAxisSkeletonExact maSk;
  if (maSk.load("opt.ma.s.json") == 0) {
    std::vector<EK::Point_3> skeletonEk;
    maSk.exportToArray(skeletonEk, skeletonEdges, skeletonTriangles);

    for (const auto &p : skeletonEk) {
      skeletonVtx.emplace_back(CGAL::to_double(p[0]), CGAL::to_double(p[1]), CGAL::to_double(p[2]));
    }
  }

  if (std::filesystem::exists("tt") == false) {
    std::filesystem::create_directories("tt");
  }

  if (std::filesystem::exists("opt") == false) {
    std::filesystem::create_directories("opt");
  }

  if (skeletonVtx.size() && skeletonEdges.size()) {
  }
  else {
    solveSkeleton(inputMesh, inputMeshBVTree, inputMeshNormal,
      inputMesh, inputMeshBVTree, inputMeshNormal,
      numTargetPoints, numIter, numAddPtsPerIter, skeletonVtx, skeletonEdges, skeletonTriangles);

    std::vector<std::array<double, 3>> p(skeletonVtx.size());
    for (int vi = 0; vi < (int)skeletonVtx.size(); vi++) {
      p[vi][0] = skeletonVtx[vi][0];
      p[vi][1] = skeletonVtx[vi][1];
      p[vi][2] = skeletonVtx[vi][2];
    }

    std::vector<std::array<int, 2>> edge(skeletonEdges.size());
    for (int i = 0; i < (int)skeletonEdges.size(); i++) {
      edge[i] = std::array<int, 2>{ skeletonEdges[i].first, skeletonEdges[i].second };
    }

    std::vector<std::array<int, 3>> tri(skeletonTriangles.size());

    nlohmann::json jj;
    jj["vtx"] = p;
    jj["edge"] = edge;
    jj["tri"] = tri;

    std::ofstream("opt.json") << jj.dump(2);
  }

  postProcessing(inputMesh, inputMeshBVTree, inputMeshNormal,
    skeletonVtx, skeletonEdges, skeletonTriangles);

  exit(1);

#if 0
  std::vector<std::vector<int>> facets1;
  for (const auto &f : facets) {
    if (f.size() > 2ull) {
      facets1.emplace_back(f);
    }
    else {
      vtx.emplace_back(vtx[f[0]] + EK::Vector_3(1e-6, 1e-6, 1e-6));
      int idnew = (int)vtx.size() - 1;
      facets1.emplace_back(f);
      facets1.back().emplace_back(idnew);
    }
  }

  std::ofstream outfile("ma1.obj");
  for (int i = 0; i < (int)vtx.size(); i++) {
    outfile << fmt::format("v {} {} {}", CGAL::to_double(vtx[i][0]), CGAL::to_double(vtx[i][1]), CGAL::to_double(vtx[i][2])) << std::endl;
  }

  for (int i = 0; i < (int)facets1.size(); i++) {
    ALOG(facets1[i].size() > 2ull);
    outfile << "f";
    for (int vi : facets1[i]) {
      outfile << " " << (vi + 1);
    }
    outfile << std::endl;
  }

  outfile.close();
#endif
  return 0;
}
