#pragma once

#include "triMeshGeo.h"

#include <vector>

namespace MedialAxisRepresentation
{
enum class SphereRefinementMethod
{
  ISOTROPIC,
  ICOSAHEDRON,
};

class TemplateSpheres
{
public:
  TemplateSpheres(int numVerticesMin, int numVerticesMax, SphereRefinementMethod srm, const char *cacheFolder = nullptr);

  int getNumTemplateSpheres() const { return (int)sphereMeshes.size(); }
  double getSphereMeshAverageEdgeLength(int id) const { return averageEdgeLength[id]; }
  int getSphereMeshIDFromEdgeLength(double edgeLen) const;
  const pgo::Mesh::TriMeshGeo &getSphereMesh(int id) const { return sphereMeshes[id]; }
  const pgo::Mesh::TriMeshGeo &getSphereMeshFromEdgeLength(double edgeLen) const { return sphereMeshes[getSphereMeshIDFromEdgeLength(edgeLen)]; }

  static pgo::Mesh::TriMeshGeo generataeSphereMesh(int numVertices);

  pgo::Mesh::TriMeshGeo createIcosahedron(double radius) const;
  void generateIcosahedronSphereMeshes(int numVerticesMin, int numVerticesMax);
  void generateIsotropicSphereMeshes(int numVerticesMin, int numVerticesMax);

  std::vector<pgo::Mesh::TriMeshGeo> sphereMeshes;
  std::vector<double> averageEdgeLength;
};

inline int MedialAxisRepresentation::TemplateSpheres::getSphereMeshIDFromEdgeLength(double edgeLen) const
{
  auto it = averageEdgeLength.begin();
  while (it != averageEdgeLength.end()) {
    if (*it <= edgeLen) {
      break;
    }
    ++it;
  }

  if (it == averageEdgeLength.end())
    --it;

  return int(it - averageEdgeLength.begin());
}

}  // namespace MedialAxisRepresentation
