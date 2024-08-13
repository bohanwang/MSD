#pragma once

#include "EigenDef.h"

#include <string>

class TriMeshGeo;

namespace MedialAxisRepresentation
{
struct VoronoiDiagram
{
  // vertices
  VegaFEM::EigenSupport::M3Xd vertexPositions;
  // edges
  std::vector<std::tuple<int, int>> edges;
  VegaFEM::EigenSupport::M3Xd edgeDirs;
  // facets
  std::vector<std::vector<int>> facets;
  std::vector<std::tuple<int, int>> facetCells;

  std::vector<std::vector<int>> cells;
};

namespace TetGenUtilities
{
int computeTetMesh(const TriMeshGeo &mesh, const std::string &switcher, VegaFEM::EigenSupport::MXd &vtx, VegaFEM::EigenSupport::MXi &tet);
}  // namespace TetGenUtilities
}  // namespace MedialAxisRepresentation