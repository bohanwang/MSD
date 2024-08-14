#include "templateSpheres.h"

#include "cgalInterface.h"
#include "pgoLogging.h"

#include "createTriMesh.h"

#include <fmt/format.h>

#include <cstring>
#include <filesystem>

namespace MedialAxisRepresentation
{
namespace ES = pgo::EigenSupport;
}

MedialAxisRepresentation::TemplateSpheres::TemplateSpheres(int numVerticesMin, int numVerticesMax, SphereRefinementMethod srm, const char *cacheFolder)
{
  if (cacheFolder && std::strlen(cacheFolder)) {
    for (int i = 0;; i++) {
      std::string filename = fmt::format("{}/sph{:03d}.obj", cacheFolder, i);
      if (std::filesystem::exists(filename) == false)
        break;

      sphereMeshes.emplace_back();
      bool ret = sphereMeshes.back().load(filename);
      PGO_ALOG(ret == true);
    }
  }

  if (sphereMeshes.size() == 0ull) {
    if (srm == SphereRefinementMethod::ICOSAHEDRON) {
      generateIcosahedronSphereMeshes(numVerticesMin, numVerticesMax);
    }
    else if (srm == SphereRefinementMethod::ISOTROPIC) {
      generateIsotropicSphereMeshes(numVerticesMin, numVerticesMax);
    }
    else {
      throw std::invalid_argument("Unknown refinement method");
    }

    if (cacheFolder && std::strlen(cacheFolder)) {
      if (!std::filesystem::exists(cacheFolder))
        std::filesystem::create_directories(cacheFolder);

      for (int i = 0; i < (int)sphereMeshes.size(); i++) {
        std::string filename = fmt::format("{}/sph{:03d}.obj", cacheFolder, i);
        sphereMeshes[i].save(filename);
      }
    }
  }

  for (const auto &m : sphereMeshes) {
    double edgeLen = 0.0;
    int numEdges = 0;
    for (int ti = 0; ti < m.numTriangles(); ti++) {
      for (int j = 0; j < 3; j++) {
        const ES::V3d &p0 = m.pos(ti, j);
        const ES::V3d &p1 = m.pos(ti, (j + 1) % 3);
        edgeLen += (p1 - p0).norm();
        numEdges += 1;
      }
    }

    averageEdgeLength.emplace_back(edgeLen / numEdges);
  }
}

pgo::Mesh::TriMeshGeo MedialAxisRepresentation::TemplateSpheres::createIcosahedron(double radius) const
{
  pgo::Mesh::TriMeshGeo mesh;

  double phi = (1.0 + std::sqrt(5.0)) * 0.5;  // golden ratio
  double a = 1.0;
  double b = 1.0 / phi;

  // add vertices
  mesh.addPos(ES::V3d(0, b, -a));
  mesh.addPos(ES::V3d(b, a, 0));
  mesh.addPos(ES::V3d(-b, a, 0));
  mesh.addPos(ES::V3d(0, b, a));
  mesh.addPos(ES::V3d(0, -b, a));
  mesh.addPos(ES::V3d(-a, 0, b));
  mesh.addPos(ES::V3d(0, -b, -a));
  mesh.addPos(ES::V3d(a, 0, -b));
  mesh.addPos(ES::V3d(a, 0, b));
  mesh.addPos(ES::V3d(-a, 0, -b));
  mesh.addPos(ES::V3d(b, -a, 0));
  mesh.addPos(ES::V3d(-b, -a, 0));

  // add triangles
  mesh.addTri(ES::V3i(2, 1, 0));
  mesh.addTri(ES::V3i(1, 2, 3));
  mesh.addTri(ES::V3i(5, 4, 3));
  mesh.addTri(ES::V3i(4, 8, 3));
  mesh.addTri(ES::V3i(7, 6, 0));
  mesh.addTri(ES::V3i(6, 9, 0));
  mesh.addTri(ES::V3i(11, 10, 4));
  mesh.addTri(ES::V3i(10, 11, 6));
  mesh.addTri(ES::V3i(9, 5, 2));
  mesh.addTri(ES::V3i(5, 9, 11));
  mesh.addTri(ES::V3i(8, 7, 1));
  mesh.addTri(ES::V3i(7, 8, 10));
  mesh.addTri(ES::V3i(2, 5, 3));
  mesh.addTri(ES::V3i(8, 1, 3));
  mesh.addTri(ES::V3i(9, 2, 0));
  mesh.addTri(ES::V3i(1, 7, 0));
  mesh.addTri(ES::V3i(11, 9, 6));
  mesh.addTri(ES::V3i(7, 10, 6));
  mesh.addTri(ES::V3i(5, 11, 4));
  mesh.addTri(ES::V3i(10, 8, 4));

  double s = radius / std::sqrt(a * a + b * b);

  for (int i = 0; i < mesh.numVertices(); i++) {
    mesh.pos(i) *= s;
  }

  return mesh;
}

void MedialAxisRepresentation::TemplateSpheres::generateIcosahedronSphereMeshes(int numVerticesMin, int numVerticesMax)
{
  pgo::Mesh::TriMeshGeo ih = createIcosahedron(1.0);

  std::vector<pgo::Mesh::TriMeshGeo> meshes;

  while (ih.numVertices() < numVerticesMax) {
    std::map<std::pair<int, int>, int> edgeVertexIDs;

    pgo::Mesh::TriMeshGeo newMesh;
    for (int vi = 0; vi < ih.numVertices(); vi++) {
      newMesh.addPos(ih.pos(vi));
    }

    std::pair<int, int> triVertexIDs[4][3] = {
      { std::pair<int, int>(0, -1),
        std::pair<int, int>(0, 1),
        std::pair<int, int>(0, 2) },

      { std::pair<int, int>(1, -1),
        std::pair<int, int>(1, 2),
        std::pair<int, int>(1, 0) },

      { std::pair<int, int>(2, -1),
        std::pair<int, int>(2, 0),
        std::pair<int, int>(2, 1) },

      { std::pair<int, int>(0, 1),
        std::pair<int, int>(1, 2),
        std::pair<int, int>(2, 0) },
    };

    for (int trii = 0; trii < ih.numTriangles(); trii++) {
      // for each sub triangle
      for (int casei = 0; casei < 4; casei++) {
        int vids[3] = { -1, -1, -1 };
        for (int j = 0; j < 3; j++) {
          std::pair<int, int> ids = triVertexIDs[casei][j];
          // if it is vertex
          if (ids.second < 0) {
            vids[j] = ih.triVtxID(trii, ids.first);
          }
          // if it is a edge
          else {
            std::pair<int, int> edgeID(ih.triVtxID(trii, ids.first), ih.triVtxID(trii, ids.second));
            if (edgeID.first > edgeID.second)
              std::swap(edgeID.first, edgeID.second);

            auto iter = edgeVertexIDs.find(edgeID);
            if (iter != edgeVertexIDs.end()) {
              vids[j] = iter->second;
            }
            else {
              ES::V3d p = (ih.pos(trii, ids.first) + ih.pos(trii, ids.second)) * 0.5;
              newMesh.addPos(p);

              vids[j] = newMesh.numVertices() - 1;
              edgeVertexIDs.emplace(edgeID, vids[j]);
            }
          }
        }

        newMesh.addTri(ES::V3i(vids[0], vids[1], vids[2]));
      }
    }

    for (int vi = 0; vi < newMesh.numVertices(); vi++) {
      ES::V3d p = newMesh.pos(vi).normalized();
      newMesh.pos(vi) = p;
    }

    meshes.emplace_back(newMesh);

    ih = newMesh;
  }

  int i = 0;
  for (; i < (int)meshes.size(); i++) {
    if (meshes[i].numVertices() >= numVerticesMin)
      break;
  }

  for (; i < (int)meshes.size(); i++) {
    sphereMeshes.emplace_back(meshes[i]);
    pgo::Mesh::TriMeshGeo mesh1 = pgo::CGALInterface::smoothMesh(sphereMeshes.back(), 10, 180.0);
    sphereMeshes.back() = mesh1;
  }
}

void MedialAxisRepresentation::TemplateSpheres::generateIsotropicSphereMeshes(int numVerticesMin, int numVerticesMax)
{
  int nv = numVerticesMin;
  while (nv < numVerticesMax) {
    pgo::Mesh::TriMeshGeo m = generataeSphereMesh(nv);
    sphereMeshes.emplace_back(m);
    nv *= 2;
  }
}

pgo::Mesh::TriMeshGeo MedialAxisRepresentation::TemplateSpheres::generataeSphereMesh(int numVertices)
{
  int approxNumTriangles = numVertices * 2;
  pgo::Mesh::TriMeshGeo mesh = pgo::Mesh::createSphereMesh(1.0, 20, 10);
  pgo::Mesh::TriMeshGeo outMesh;

  double triArea = 4 * M_PI / approxNumTriangles;
  double triEdgeLength = std::sqrt(triArea * 2 / std::sqrt(3));

  pgo::Mesh::TriMeshGeo meshTemp = pgo::CGALInterface::isotropicRemeshing(mesh, triEdgeLength, 5, 180.0);
  double left, right;
  if (numVertices < meshTemp.numVertices()) {
    right = triEdgeLength;
    left = triEdgeLength * 2;

    while (1) {
      meshTemp = pgo::CGALInterface::isotropicRemeshing(mesh, left, 5, 180.0);
      if (meshTemp.numVertices() == numVertices) {
        outMesh = meshTemp;
        for (int vi = 0; vi < outMesh.numVertices(); vi++) {
          outMesh.pos(vi) /= outMesh.pos(vi).norm();
        }

        return outMesh;
      }
      else if (meshTemp.numVertices() < numVertices) {
        break;
      }
      else {
        left *= 2.0;
      }
    }
  }
  else if (numVertices > meshTemp.numVertices()) {
    left = triEdgeLength;
    right = triEdgeLength * 0.5;

    while (1) {
      meshTemp = pgo::CGALInterface::isotropicRemeshing(mesh, right, 5, 180.0);
      if (meshTemp.numVertices() == numVertices) {
        outMesh = meshTemp;
        for (int vi = 0; vi < outMesh.numVertices(); vi++) {
          outMesh.pos(vi) /= outMesh.pos(vi).norm();
        }

        return outMesh;
      }
      else if (meshTemp.numVertices() > numVertices) {
        break;
      }
      else {
        right *= 0.5;
      }
    }
  }
  else {
    outMesh = meshTemp;
    for (int vi = 0; vi < outMesh.numVertices(); vi++) {
      outMesh.pos(vi) /= outMesh.pos(vi).norm();
    }

    return outMesh;
  }

  int numIter = 5;
  while (left > right && numIter >= 0) {
    double mid = (left + right) * 0.5;
    meshTemp = pgo::CGALInterface::isotropicRemeshing(mesh, mid, 5, 180.0);
    if (meshTemp.numVertices() > numVertices) {
      right = mid;
    }
    else if (meshTemp.numVertices() < numVertices) {
      left = mid;
    }
    else {
      outMesh = meshTemp;
      for (int vi = 0; vi < outMesh.numVertices(); vi++) {
        outMesh.pos(vi) /= outMesh.pos(vi).norm();
      }

      return outMesh;
    }

    numIter--;
  }

  outMesh = meshTemp;
  for (int vi = 0; vi < outMesh.numVertices(); vi++) {
    outMesh.pos(vi) /= outMesh.pos(vi).norm();
  }

  return outMesh;
}
