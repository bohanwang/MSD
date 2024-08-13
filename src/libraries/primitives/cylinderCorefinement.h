#pragma once

#include "triMeshGeo.h"
#include "boundingVolumeTree.h"
#include "triMeshPseudoNormal.h"

namespace MedialAxisRepresentation
{
int corefineCylinderMeshWithTarget(const pgo::Mesh::TriMeshGeo &sphereMesh, const pgo::EigenSupport::V3d centers[2],
  const pgo::Mesh::TriMeshGeo &targetMesh, const pgo::Mesh::TriMeshBVTree &targetMeshBVTree, const pgo::Mesh::TriMeshPseudoNormal &targetMeshNormals,
  double maxConfidentDist, double maxAllowedRadiusDist, pgo::Mesh::TriMeshGeo &meshOut, const char *filename = nullptr);

int fillCylinder(const pgo::Mesh::TriMeshGeo &sphereMesh, const pgo::EigenSupport::V3d centers[2], pgo::Mesh::TriMeshGeo &meshOut);
}  // namespace MedialAxisRepresentation