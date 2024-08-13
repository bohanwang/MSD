#pragma once

#include "triMeshGeo.h"
#include "boundingVolumeTree.h"
#include "triMeshPseudoNormal.h"

namespace MedialAxisRepresentation
{
int corefinePrismMeshWithTarget(const pgo::Mesh::TriMeshGeo &prismMesh, const pgo::EigenSupport::V3d centers[3],
  const pgo::Mesh::TriMeshGeo &targetMesh, const pgo::Mesh::TriMeshBVTree &targetMeshBVTree, const pgo::Mesh::TriMeshPseudoNormal &targetMeshNormals,
  double maxConfidentDist, double maxAllowedRadiusDist, pgo::Mesh::TriMeshGeo &meshOut, const char *filename = nullptr);

void fillPrism(const pgo::Mesh::TriMeshGeo &prismMesh, const pgo::EigenSupport::V3d centers[3], pgo::Mesh::TriMeshGeo &meshOut);
}  // namespace MedialAxisRepresentation