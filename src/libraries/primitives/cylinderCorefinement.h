#pragma once

#include "triMeshGeo.h"
#include "boundingVolumeTree.h"
#include "triMeshPseudoNormal.h"

namespace MedialAxisRepresentation
{
int corefineCylinderMeshWithTarget(const TriMeshGeo &sphereMesh, const Vec3d centers[2],
  const TriMeshGeo &targetMesh, const TriMeshBVTree &targetMeshBVTree, const TriMeshPseudoNormal &targetMeshNormals,
  double maxConfidentDist, double maxAllowedRadiusDist, TriMeshGeo &meshOut, const char *filename = nullptr);


int fillCylinder(const TriMeshGeo &sphereMesh, const Vec3d centers[2], TriMeshGeo &meshOut);
}  // namespace MedialAxisRepresentation