#pragma once

#include "skeletonExact.h"

#include "EigenDef.h"
#include "triMeshGeo.h"
#include "boundingVolumeTree.h"
#include "triMeshPseudoNormal.h"

namespace MedialAxisRepresentation
{
void solveSkeleton(const pgo::Mesh::TriMeshGeo &targetMesh, const pgo::Mesh::TriMeshBVTree &targetMeshBVTree, const pgo::Mesh::TriMeshPseudoNormal &targetNormals,
  const pgo::Mesh::TriMeshGeo &targetMeshSmall, const pgo::Mesh::TriMeshBVTree &targetMeshBVTreeSmall, const pgo::Mesh::TriMeshPseudoNormal &targetNormalsSmall,
  const std::string &maInitFilename, int nPt, int nIt, int numAddedPt,
  std::vector<pgo::EigenSupport::V3d> &finalSkeletonPoints,
  std::vector<std::pair<int, int>> &finalSkeletonEdges,
  std::vector<std::tuple<int, int, int>> &finalSkeletonTriangles,
  std::vector<pgo::Mesh::TriMeshGeo> *finalFitMeshes = nullptr);
}  // namespace MedialAxisRepresentation
