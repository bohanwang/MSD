#pragma once

#include "skeletonExact.h"

#include "triMeshGeo.h"
#include "boundingVolumeTree.h"
#include "EigenDef.h"

namespace MedialAxisRepresentation
{
void computeRVD_internal_Dijkstras_withTri_refine_withEdge(
  const pgo::Mesh::TriMeshGeo &targetMesh, const pgo::Mesh::TriMeshBVTree &targetMeshBVTree,
  const pgo::Mesh::TriMeshGeo &targetMeshSmall, const pgo::Mesh::TriMeshBVTree &targetMeshBVTreeSmall,
  const std::vector<pgo::EigenSupport::V3d> &samplePoints,
  std::vector<EK::Point_3> &medialAxisVertices, std::vector<std::vector<int>> &medialAxisFacets,
  std::vector<pgo::EigenSupport::V3d> &finalSkeletonPoints, std::vector<std::pair<int, int>> &finalSkeletonEdges, std::vector<std::tuple<int, int, int>> &finalSkeletonTriangles,
  std::vector<pgo::EigenSupport::V3d> *cellCenters, const std::string &saveExactResults = "");

void preprocessingMA(const std::vector<EK::Point_3> &medialAxisVertices, const std::vector<std::vector<int>> &medialAxisFacets,
  const pgo::Mesh::TriMeshGeo &mesh, const pgo::Mesh::TriMeshBVTree &meshBVTree,
  std::vector<EK::Point_3> &finalMedialAxisVertices, std::vector<std::vector<int>> &finalMedialAxisFacets);
}  // namespace MedialAxisRepresentation
