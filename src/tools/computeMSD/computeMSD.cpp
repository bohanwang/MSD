#include "skeletons/skeletonExact.h"
#include "skeletons/medialSkeletalDiagramOpt.h"
#include "primitives/primitiveICP.h"
#include "primitives/sphereCorefinement.h"
#include "primitives/cylinderCorefinement.h"
#include "primitives/prismCorefinement.h"

#include "pgoLogging.h"
#include "initPredicates.h"
#include "geogramInterface.h"
#include "triMeshGeo.h"
#include "triMeshPseudoNormal.h"
#include "boundingVolumeTree.h"
#include "EigenSupport.h"
#include "geometryQuery.h"
#include "createTriMesh.h"

#include <nlohmann/json.hpp>

#include <tbb/parallel_for.h>

#include <argparse/argparse.hpp>

#include <iostream>

static void skeletonPostProcessing(const pgo::Mesh::TriMeshGeo &inputMesh, const pgo::Mesh::TriMeshBVTree &inputMeshBVTree, const pgo::Mesh::TriMeshPseudoNormal &inputMeshNormal, bool removePointOnLine,
  std::vector<pgo::EigenSupport::V3d> &finalSkeletonPoints, std::vector<std::pair<int, int>> &finalSkeletonEdges, std::vector<std::tuple<int, int, int>> &finalSkeletonTriangles,
  std::map<int, int> &vidMap, std::map<int, std::pair<int, int>> &edgeMap, std::map<int, std::tuple<int, int, int>> &triMap);

int main(int argc, char *argv[])
{
  namespace ES = pgo::EigenSupport;
  using namespace MedialAxisRepresentation;

  argparse::ArgumentParser program("Compute Medial Skeletal Diagram");

  program.add_argument("-i", "--input-mesh")
    .help("Input surface mesh filename")
    .required()
    .metavar("PATH");
  ;

  program.add_argument("-a", "--input-ma")
    .help("Input surface mesh medial axis filename")
    .required()
    .metavar("PATH");

  program.add_argument("-n", "--num-vertices")
    .help("The number of skeleton vertices requested")
    .required()
    .scan<'i', int>()
    .metavar("INT");

  program.add_argument("-e", "--num-additional-solve")
    .help("The number of additional solving after the first phase. By default (=1), it does not do additional solving.")
    .default_value(1)
    .scan<'i', int>()
    .metavar("INT");

  program.add_argument("-N", "--num-additional-vertices")
    .help("The number of additional vertices to be added per iteration during additional solving after the first phase.")
    .default_value(1)
    .scan<'i', int>()
    .metavar("INT");

  try {
    program.parse_args(argc, argv);  // Example: ./main --color orange
  }
  catch (const std::exception &err) {
    std::cerr << err.what() << std::endl;
    std::cerr << program;
    return 1;
  }

  std::string maFilename = program.get<std::string>("--input-ma");
  std::string inputMeshFilename = program.get<std::string>("--input-mesh");
  int numTargetPoints = program.get<int>("--num-vertices");
  int numIter = program.get<int>("--num-additional-solve");
  int numAddPtsPerIter = program.get<int>("--num-additional-vertices");

  if constexpr (1) {
    std::ofstream cfile("cmd.txt");
    for (int i = 0; i < argc; i++) {
      cfile << argv[i] << ' ';
    }
  }

  pgo::Logging::init();
  pgo::Mesh::initPredicates();
  pgo::GeogramInterface::initGEO();

  std::cout << "Surface Mesh Filename: " << inputMeshFilename << '\n'
            << "MA Filename: " << maFilename << '\n'
            << "#tgt points: " << numTargetPoints << '\n'
            << "#iter: " << numIter << '\n'
            << "#added pt: " << numAddPtsPerIter << std::endl;

  std::set_terminate([]() {
    std::cout << "Unhandled exception\n"
              << std::flush;
    std::abort();
  });

  pgo::Mesh::TriMeshGeo inputMesh;
  if (inputMesh.load(inputMeshFilename) != true) {
    std::cerr << "Cannot load " << inputMeshFilename << std::endl;
    return 1;
  }

  pgo::Mesh::TriMeshBVTree inputMeshBVTree;
  inputMeshBVTree.buildByInertiaPartition(inputMesh);

  pgo::Mesh::TriMeshPseudoNormal inputMeshNormal;
  inputMeshNormal.buildPseudoNormals(inputMesh);

  std::vector<ES::V3d> skeletonVtx;
  std::vector<std::pair<int, int>> skeletonEdges;
  std::vector<std::tuple<int, int, int>> skeletonTriangles;

  SkeletonExact maSk;
  if (maSk.load("opt.ma.json") == 0) {
    std::vector<EK::Point_3> skeletonEk;
    maSk.exportToArray(skeletonEk, skeletonEdges, skeletonTriangles);

    for (const auto &p : skeletonEk) {
      skeletonVtx.emplace_back(CGAL::to_double(p[0]), CGAL::to_double(p[1]), CGAL::to_double(p[2]));
    }
  }

  if (skeletonVtx.size() && (skeletonEdges.size() || skeletonTriangles.size())) {
    std::cout << "Skeleton loaded." << std::endl;
  }
  else {
    solveSkeleton(inputMesh, inputMeshBVTree, inputMeshNormal,
      inputMesh, inputMeshBVTree, inputMeshNormal,
      maFilename, numTargetPoints, numIter, numAddPtsPerIter,
      skeletonVtx, skeletonEdges, skeletonTriangles);
  }

  std::map<int, int> vidMap;
  std::map<int, std::pair<int, int>> edgeMap;
  std::map<int, std::tuple<int, int, int>> triMap;

  skeletonPostProcessing(inputMesh, inputMeshBVTree, inputMeshNormal, true,
    skeletonVtx, skeletonEdges, skeletonTriangles,
    vidMap, edgeMap, triMap);

  return 0;
}

void skeletonPostProcessing(const pgo::Mesh::TriMeshGeo &inputMesh, const pgo::Mesh::TriMeshBVTree &inputMeshBVTree, const pgo::Mesh::TriMeshPseudoNormal &inputMeshNormal, bool removePointOnLine,
  std::vector<pgo::EigenSupport::V3d> &finalSkeletonPoints, std::vector<std::pair<int, int>> &finalSkeletonEdges, std::vector<std::tuple<int, int, int>> &finalSkeletonTriangles,
  std::map<int, int> &vidMap, std::map<int, std::pair<int, int>> &edgeMap, std::map<int, std::tuple<int, int, int>> &triMap)
{
  namespace ES = pgo::EigenSupport;
  using namespace MedialAxisRepresentation;

  SkeletonExact maPost;
  std::vector<pgo::EigenSupport::V3d> skeletonPoints = finalSkeletonPoints;

  if (maPost.load("opt.ma.post.json") == 0) {
  }
  else {
    for (int i = 0; i < (int)finalSkeletonEdges.size(); i++) {
      if (finalSkeletonEdges[i].first > finalSkeletonEdges[i].second) {
        std::swap(finalSkeletonEdges[i].first, finalSkeletonEdges[i].second);
      }
    }

    std::set<std::pair<int, int>> edgeSet;
    for (const auto &e : finalSkeletonEdges)
      edgeSet.emplace(e);

    std::set<std::array<int, 3>, ES::IntArrayLess<3>> triSet;
    for (int ti = 0; ti < (int)finalSkeletonTriangles.size(); ti++) {
      std::array<int, 3> vids{
        std::get<0>(finalSkeletonTriangles[ti]),
        std::get<1>(finalSkeletonTriangles[ti]),
        std::get<2>(finalSkeletonTriangles[ti]),
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

    finalSkeletonEdges.assign(edgeSet.begin(), edgeSet.end());
    finalSkeletonTriangles.clear();
    for (const auto &tri : triSet) {
      finalSkeletonTriangles.emplace_back(tri[0], tri[1], tri[2]);
    }
    std::cout << "# duplicated edge removed:" << edgeDel << std::endl;

    std::vector<EK::Point_3> vtxEK;
    for (int i = 0; i < finalSkeletonPoints.size(); i++) {
      vtxEK.emplace_back(finalSkeletonPoints[i][0], finalSkeletonPoints[i][1], finalSkeletonPoints[i][2]);
    }

    maPost.initFromInput(vtxEK, finalSkeletonEdges, finalSkeletonTriangles);
    std::cout << "is connected:" << maPost.isConnected() << std::endl;
    std::cout << "is consistent:" << maPost.isConsistent() << std::endl;

    if (removePointOnLine) {
      while (1) {
        bool modified = false;
        for (const auto &pr : maPost.vertices) {
          if (pr.second.neighboringEdges.size() == 2ull &&
            pr.second.neighboringTriangles.size() == 0ull) {
            EK::Point_3 p0 = pr.second.pos;
            ES::V3d pt(CGAL::to_double(p0[0]), CGAL::to_double(p0[1]), CGAL::to_double(p0[2]));

            const auto &e0 = maPost.getE(pr.second.neighboringEdges[0]);
            const auto &e1 = maPost.getE(pr.second.neighboringEdges[1]);

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

            EK::Point_3 p1 = maPost.getV(vids[0]).pos;
            EK::Point_3 p2 = maPost.getV(vids[1]).pos;

            EK::Vector_3 d0 = p1 - p0;
            EK::Vector_3 d1 = p2 - p0;

            double d0_norm = std::sqrt(CGAL::to_double(d0.squared_length()));
            double d1_norm = std::sqrt(CGAL::to_double(d1.squared_length()));
            double dot_d0_d1 = CGAL::to_double(CGAL::scalar_product(d0, d1));
            double cosAngle = dot_d0_d1 / (d0_norm * d1_norm);

            if (cosAngle < std::cos(179.0 / 180.0 * M_PI)) {
              continue;
            }

            ES::V3d v0(CGAL::to_double(p1[0]), CGAL::to_double(p1[1]), CGAL::to_double(p1[2]));
            ES::V3d v1(CGAL::to_double(p2[0]), CGAL::to_double(p2[1]), CGAL::to_double(p2[2]));

            if (inputMeshBVTree.hasLineSegmentIntersectionExact(inputMesh, v0, v1)) {
              continue;
            }

            if (maPost.canCollapseEdge(pr.second.neighboringEdges[0]) == false) {
              continue;
            }

            maPost.collapseEdge(pr.second.neighboringEdges[0], p1);
            modified = true;
            break;
          }
        }

        if (modified == false)
          break;
      }
    }  // end if

    maPost.save("opt.ma.post.json");
    maPost.saveDisplayMesh("post-sk.obj");
  }

  if constexpr (1) {
    std::vector<EK::Point_3> vtxEK;
    std::vector<std::pair<int, int>> se;
    std::vector<std::tuple<int, int, int>> st;

    maPost.exportToArray(vtxEK, se, st, vidMap, edgeMap, triMap);

    skeletonPoints.clear();
    skeletonPoints.resize(vtxEK.size());
    for (size_t i = 0; i < vtxEK.size(); i++) {
      skeletonPoints[i][0] = CGAL::to_double(vtxEK[i][0]);
      skeletonPoints[i][1] = CGAL::to_double(vtxEK[i][1]);
      skeletonPoints[i][2] = CGAL::to_double(vtxEK[i][2]);
    }

    finalSkeletonEdges = se;
    finalSkeletonTriangles = st;

    // export vidMap to json
    nlohmann::json j_map(vidMap);
    std::ofstream o(fmt::format("vidMap.json"));
    o << std::setw(4) << j_map << std::endl;

    // export edgeMap to json
    nlohmann::json j_edgeMap(edgeMap);
    std::ofstream o2(fmt::format("edgeMap.json"));
    o2 << std::setw(4) << j_edgeMap << std::endl;

    // export triMap to json
    nlohmann::json j_triMap(triMap);
    std::ofstream o3(fmt::format("triMap.json"));
    o3 << std::setw(4) << j_triMap << std::endl;
  }
}

void primitiveFitting(const pgo::Mesh::TriMeshGeo &inputMesh, const pgo::Mesh::TriMeshBVTree &inputMeshBVTree, const pgo::Mesh::TriMeshPseudoNormal &inputMeshNormal,
  const std::vector<pgo::EigenSupport::V3d> &skeletonPoints, std::vector<std::pair<int, int>> finalSkeletonEdges, std::vector<std::tuple<int, int, int>> finalSkeletonTriangles,
  const std::map<int, int> &vidMap, const std::map<int, std::pair<int, int>> &edgeMap, const std::map<int, std::tuple<int, int, int>> &triMap)
{
  namespace ES = pgo::EigenSupport;
  using namespace MedialAxisRepresentation;

  pgo::Mesh::TriMeshGeo finalMeshes;
  pgo::Mesh::TriMeshGeo corefFinalMeshes;

  double avgLength = 0;
  double count = 0;
  for (int ti = 0; ti < inputMesh.numTriangles(); ti++) {
    for (int j = 0; j < 3; j++) {
      avgLength += (inputMesh.pos(ti, j) - inputMesh.pos(ti, (j + 1) % 3)).norm();
      count += 1.0;
    }
  }
  avgLength /= count;

  for (int i = 0; i < (int)finalSkeletonEdges.size(); i++) {
    if (finalSkeletonEdges[i].first > finalSkeletonEdges[i].second) {
      std::swap(finalSkeletonEdges[i].first, finalSkeletonEdges[i].second);
    }
  }

  int edgeDel = 0;
  for (int ti = 0; ti < (int)finalSkeletonTriangles.size(); ti++) {
    int vids[3] = {
      std::get<0>(finalSkeletonTriangles[ti]),
      std::get<1>(finalSkeletonTriangles[ti]),
      std::get<2>(finalSkeletonTriangles[ti]),
    };

    for (int j = 0; j < 3; j++) {
      std::pair edgeID{ vids[j], vids[(j + 1) % 3] };
      if (edgeID.first > edgeID.second) {
        std::swap(edgeID.first, edgeID.second);
      }

      for (auto it = finalSkeletonEdges.begin(); it != finalSkeletonEdges.end(); ++it) {
        if (*it == edgeID) {
          finalSkeletonEdges.erase(it);
          edgeDel++;
          break;
        }
      }
    }
  }

  int nIter = 25;
  double sc = 1e-1;
  double dilateSize = avgLength * 2;

  int vEnd = (int)skeletonPoints.size();
  int eEnd = (int)finalSkeletonEdges.size();
  int tEnd = (int)finalSkeletonTriangles.size();
  int coref = 1;

  std::cout << vEnd << ',' << eEnd << ',' << tEnd << std::endl;

  for (int vi = 0; vi < vEnd; vi++) {
    int numCenters = 1;
    ES::MXd centers(numCenters, 3);
    ES::VXd centerRadii(numCenters);

    centers.row(0) = skeletonPoints[vi];

    double avgR = 0;
    for (int j = 0; j < numCenters; j++) {
      ES::V3d p = centers.row(j);
      auto ret = inputMeshBVTree.closestTriangleQuery(inputMesh, p);
      centerRadii[j] = std::sqrt(ret.dist2);
      avgR += centerRadii[j];
    }
    avgR /= (double)numCenters;

    MedialAxisRepresentation::PrimitiveICPSolverParameters primitiveICPParams;
    primitiveICPParams.sphereRadius = avgR;
    primitiveICPParams.useCorefine = false;
    primitiveICPParams.verbose = false;
    primitiveICPParams.maxPDNumIter = 15;
    primitiveICPParams.maxNumIter = nIter;
    primitiveICPParams.smoothnessCoeff = sc;
    primitiveICPParams.expansionCoeff = 2;
    primitiveICPParams.contactCoeff = 1e4;

    pgo::Mesh::TriMeshGeo fittingMesh;
    MedialAxisRepresentation::PrimitiveICPSolver solver(centers, centerRadii, primitiveICPParams, numCenters);
    int ret = solver.sim(centers, inputMesh, inputMeshBVTree, inputMeshNormal, fittingMesh);
    if (ret == 0) {
      finalMeshes.addMesh(fittingMesh);

      fittingMesh.save("last-icp.obj");
      std::ofstream("last-icp.txt") << skeletonPoints[vi][0] << ',' << skeletonPoints[vi][1] << ',' << skeletonPoints[vi][2];

      fittingMesh.save(fmt::format("fit.1-v{:04d}.obj", vi));

      if (coref) {
        pgo::Mesh::TriMeshGeo meshOut;
        corefineSphereMeshWithTarget(fittingMesh, skeletonPoints[vi].data(), inputMesh,
          inputMeshBVTree, inputMeshNormal, 1e-3, dilateSize, meshOut);

        meshOut.save(fmt::format("coref-v{:04d}.obj", vi));
        corefFinalMeshes.addMesh(meshOut);
      }
    }
  }

  // finalSkeletonEdges.size()
  for (int ei = 0; ei < eEnd; ei++) {
    int numCenters = 2;
    ES::MXd centers(numCenters, 3);
    ES::VXd centerRadii(numCenters);

    centers.row(0) = skeletonPoints[finalSkeletonEdges[ei].first];
    centers.row(1) = skeletonPoints[finalSkeletonEdges[ei].second];

    ES::V3d p0 = skeletonPoints[finalSkeletonEdges[ei].first];
    ES::V3d p1 = skeletonPoints[finalSkeletonEdges[ei].second];

    if (inputMeshBVTree.hasLineSegmentIntersectionExact(inputMesh, p0, p1)) {
      std::cout << "skeleton intersected with surface. Skip" << std::endl;
      pgo::Mesh::createSingleTriangleMesh(p0, p1, p0 + pgo::asVec3d(1e-6)).save("inter.obj");
      continue;
    }

    double avgR = 0;
    double minR = 1e100;
    int nSamples = 20;
    for (int j = 0; j < nSamples; j++) {
      double ratio = double(j) / double(nSamples - 1);
      ES::V3d p = centers.row(0) * (1 - ratio) + centers.row(1) * ratio;
      auto ret = inputMeshBVTree.closestTriangleQuery(inputMesh, p);
      // avgR = std::min(avgR, std::sqrt(ret.dist2));

      if (j == 0)
        centerRadii(0) = std::sqrt(ret.dist2);
      else if (j == nSamples - 1)
        centerRadii(1) = std::sqrt(ret.dist2);

      avgR += std::sqrt(ret.dist2);
      minR = std::min(minR, std::sqrt(ret.dist2));
    }

    avgR /= nSamples;

    // minR = std::max(minR, 0.01);
    // centerRadii(0) = minR;
    // centerRadii(1) = minR;

    // centerRadii(0) = avgR;
    // centerRadii(1) = avgR;

    MedialAxisRepresentation::PrimitiveICPSolverParameters primitiveICPParams;
    primitiveICPParams.sphereRadius = (centerRadii(0) + centerRadii(1)) * 0.5;
    // primitiveICPParams.sphereRadius = minR;
    primitiveICPParams.useCorefine = false;
    primitiveICPParams.verbose = false;
    primitiveICPParams.maxPDNumIter = 15;
    primitiveICPParams.maxNumIter = nIter;
    primitiveICPParams.smoothnessCoeff = sc;
    primitiveICPParams.expansionCoeff = 2;
    primitiveICPParams.contactCoeff = 1e4;

    pgo::Mesh::TriMeshGeo fittingMesh;
    MedialAxisRepresentation::PrimitiveICPSolver solver(centers, centerRadii, primitiveICPParams, numCenters);
    int ret = solver.sim(centers, inputMesh, inputMeshBVTree, inputMeshNormal, fittingMesh);

    if (ret == 0) {
      finalMeshes.addMesh(fittingMesh);

      ES::V3d centers_v[2] = {
        ES::V3d(centers(0, 0), centers(0, 1), centers(0, 2)),
        ES::V3d(centers(1, 0), centers(1, 1), centers(1, 2))
      };

      fittingMesh.save("last-icp.obj");
      std::ofstream("last-icp.txt")
        << centers_v[0][0] << ',' << centers_v[0][1] << ',' << centers_v[0][2] << '\n'
        << centers_v[1][0] << ',' << centers_v[1][1] << ',' << centers_v[1][2];

      fittingMesh.save(fmt::format("fit.1-e{:04d}.obj", ei));

      pgo::Mesh::TriMeshGeo mOut;
      fillCylinder(fittingMesh, centers_v, mOut);
      mOut.save(fmt::format("fit.1-e{:04d}.obj", ei));

      if (coref) {
        pgo::Mesh::TriMeshGeo meshOut;
        corefineCylinderMeshWithTarget(fittingMesh, centers_v, inputMesh,
          inputMeshBVTree, inputMeshNormal, 1e-3, dilateSize, meshOut);

        meshOut.save(fmt::format("coref-e{:04d}.obj", ei));

        corefFinalMeshes.addMesh(meshOut);

        // std::vector<double> centersDump = {
        //   centers_v[0][0],
        //   centers_v[0][1],
        //   centers_v[0][2],
        //   centers_v[1][0],
        //   centers_v[1][1],
        //   centers_v[1][2],
        // };

        // nlohmann::json jOut;
        // jOut["centers"] = centersDump;
        // std::ofstream ofs(fmt::format("coref-e{:04d}.json", ei));
        // ofs << jOut << std::endl;
      }
    }
  }

  // finalSkeletonTriangles.size()
  for (int ti = 0; ti < tEnd; ti++) {
    int numCenters = 3;
    ES::MXd centers(numCenters, 3);
    ES::VXd centerRadii(numCenters);

    centers.row(0) = skeletonPoints[std::get<0>(finalSkeletonTriangles[ti])];
    centers.row(1) = skeletonPoints[std::get<1>(finalSkeletonTriangles[ti])];
    centers.row(2) = skeletonPoints[std::get<2>(finalSkeletonTriangles[ti])];

    pgo::Mesh::TriMeshGeo aa = pgo::Mesh::createSingleTriangleMesh(
      skeletonPoints[std::get<0>(finalSkeletonTriangles[ti])],
      skeletonPoints[std::get<1>(finalSkeletonTriangles[ti])],
      skeletonPoints[std::get<2>(finalSkeletonTriangles[ti])]);

    std::vector<int> l;
    inputMeshBVTree.triangleIntersectionExact(inputMesh, aa.pos(0), aa.pos(1), aa.pos(2), l);
    if (l.size()) {
      std::cout << "skeleton intersected with surface. Skip" << std::endl;
      aa.save("mad.obj");

      continue;
    }

    double avgR = 1e100;
    for (int j = 0; j < 3; j++) {
      ES::V3d p = centers.row(j);
      auto ret = inputMeshBVTree.closestTriangleQuery(inputMesh, p);
      // avgR = std::min(avgR, std::sqrt(ret.dist2));
      centerRadii(j) = std::sqrt(ret.dist2);
    }
    // centerRadii(0) = avgR;
    // centerRadii(1) = avgR;

    MedialAxisRepresentation::PrimitiveICPSolverParameters primitiveICPParams;
    primitiveICPParams.sphereRadius = (centerRadii[0] + centerRadii[1] + centerRadii[2]) / 3.0;
    primitiveICPParams.useCorefine = false;
    primitiveICPParams.verbose = false;
    primitiveICPParams.maxPDNumIter = 15;
    primitiveICPParams.maxNumIter = nIter;
    primitiveICPParams.smoothnessCoeff = sc;
    primitiveICPParams.expansionCoeff = 2;
    primitiveICPParams.contactCoeff = 1e4;

    pgo::Mesh::TriMeshGeo fittingMesh;
    MedialAxisRepresentation::PrimitiveICPSolver solver(centers, centerRadii, primitiveICPParams, numCenters);
    int ret = solver.sim(centers, inputMesh, inputMeshBVTree, inputMeshNormal, fittingMesh);

    if (ret == 0) {
      finalMeshes.addMesh(fittingMesh);

      ES::V3d centers_v[3] = {
        ES::V3d(centers(0, 0), centers(0, 1), centers(0, 2)),
        ES::V3d(centers(1, 0), centers(1, 1), centers(1, 2)),
        ES::V3d(centers(2, 0), centers(2, 1), centers(2, 2)),
      };

      fittingMesh.save("last-icp.obj");
      std::ofstream("last-icp.txt")
        << centers_v[0][0] << ',' << centers_v[0][1] << ',' << centers_v[0][2] << '\n'
        << centers_v[1][0] << ',' << centers_v[1][1] << ',' << centers_v[1][2] << '\n'
        << centers_v[2][0] << ',' << centers_v[2][1] << ',' << centers_v[2][2];

      fittingMesh.save(fmt::format("fit.1-t{:04d}.obj", ti));

      pgo::Mesh::TriMeshGeo mOut;
      fillPrism(fittingMesh, centers_v, mOut);
      mOut.save(fmt::format("fit.1-t{:04d}.obj", ti));

      if (coref) {
        pgo::Mesh::TriMeshGeo meshOut;
        corefinePrismMeshWithTarget(fittingMesh, centers_v, inputMesh,
          inputMeshBVTree, inputMeshNormal, 1e-3, dilateSize, meshOut);
        meshOut.save(fmt::format("coref-t{:04d}.obj", ti));

        corefFinalMeshes.addMesh(meshOut);

        // std::vector<double> centersDump = {
        //   centers_v[0][0],
        //   centers_v[0][1],
        //   centers_v[0][2],
        //   centers_v[1][0],
        //   centers_v[1][1],
        //   centers_v[1][2],
        //   centers_v[2][0],
        //   centers_v[2][1],
        //   centers_v[2][2],
        // };

        // nlohmann::json jOut;
        // jOut["centers"] = centersDump;
        // std::ofstream ofs(fmt::format("coref-t{:04d}.json", ti));
        // ofs << jOut << std::endl;
      }
    }
  }

  finalMeshes.save("exp.obj");
  corefFinalMeshes.save("exp-coref.obj");
}
