#pragma once

#include "triMeshGeo.h"
#include "boundingVolumeTree.h"
#include "triMeshPseudoNormal.h"

namespace MedialAxisRepresentation
{
int corefineSphereMeshWithTarget(const pgo::Mesh::TriMeshGeo &sphereMesh, const double center[3],
  const pgo::Mesh::TriMeshGeo &targetMesh, const pgo::Mesh::TriMeshBVTree &targetMeshBVTree, const pgo::Mesh::TriMeshPseudoNormal &targetMeshNormals,
  double maxConfidentDist, double maxAllowedRadiusDist, pgo::Mesh::TriMeshGeo &meshOut, const char *filename = nullptr);

int remeshSphereMeshWithTarget(const pgo::Mesh::TriMeshGeo &sphereMesh, const double center[3],
  const pgo::Mesh::TriMeshGeo &targetMesh, const pgo::Mesh::TriMeshBVTree &targetMeshBVTree, const pgo::Mesh::TriMeshPseudoNormal &targetMeshNormals,
  double maxConfidentDist, double maxAllowedRadiusDist, pgo::Mesh::TriMeshGeo &meshOut);
}  // namespace MedialAxisRepresentation