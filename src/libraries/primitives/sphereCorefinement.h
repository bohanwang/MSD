#pragma once

#include "triMeshGeo.h"
#include "boundingVolumeTree.h"
#include "triMeshPseudoNormal.h"

namespace MedialAxisRepresentation
{
int corefineSphereMeshWithTarget(const TriMeshGeo &sphereMesh, const double center[3],
  const TriMeshGeo &targetMesh, const TriMeshBVTree &targetMeshBVTree, const TriMeshPseudoNormal &targetMeshNormals,
  double maxConfidentDist, double maxAllowedRadiusDist, TriMeshGeo &meshOut, const char *filename = nullptr);

int remeshSphereMeshWithTarget(const TriMeshGeo &sphereMesh, const double center[3],
  const TriMeshGeo &targetMesh, const TriMeshBVTree &targetMeshBVTree, const TriMeshPseudoNormal &targetMeshNormals,
  double maxConfidentDist, double maxAllowedRadiusDist, TriMeshGeo &meshOut);
}  // namespace MedialAxisRepresentation