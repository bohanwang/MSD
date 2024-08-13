#pragma once

#include "triMeshGeo.h"
#include "boundingVolumeTree.h"
#include "triMeshPseudoNormal.h"

namespace MedialAxisRepresentation
{
int corefinePrismMeshWithTarget(const TriMeshGeo &prismMesh, const Vec3d centers[3],
  const TriMeshGeo &targetMesh, const TriMeshBVTree &targetMeshBVTree, const TriMeshPseudoNormal &targetMeshNormals,
  double maxConfidentDist, double maxAllowedRadiusDist, TriMeshGeo &meshOut, const char *filename = nullptr);

void fillPrism(const TriMeshGeo &prismMesh, const Vec3d centers[3], TriMeshGeo &meshOut);
}  // namespace MedialAxisRepresentation