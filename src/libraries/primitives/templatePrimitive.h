#pragma once

#include "triMeshGeo.h"
#include "triMeshNeighbor.h"
#include "triMeshPseudoNormal.h"

namespace MedialAxisRepresentation
{
void rotAlignMat(const pgo::EigenSupport::V3d &v0, const pgo::EigenSupport::V3d &v1, pgo::EigenSupport::M3d &rot);

pgo::Mesh::TriMeshGeo createPrismWallMeshForFigures(
  const pgo::EigenSupport::V3d &e1, const pgo::EigenSupport::V3d &e2, const pgo::EigenSupport::V3d &e3, const double thickness, const double approxTargetRadius,
  std::vector<pgo::EigenSupport::V3d> &prismVtxDri, std::vector<pgo::EigenSupport::V3d> &baryCentricWeights, std::vector<int> &rotAxisID);

class TemplatePrimitive
{
public:
  TemplatePrimitive(){};
  TemplatePrimitive(double radius_, int numCylinderSlices_, int primitiveType_):
    radius(radius_), numCylinderSlices(numCylinderSlices_), primitiveType(primitiveType_){};
  ~TemplatePrimitive(){};

  virtual void init(const pgo::EigenSupport::MXd &centers, const pgo::EigenSupport::VXd &centerRadii, const double useHardCodedTargetEdgeLenForTri = -1);
  virtual void init(const pgo::EigenSupport::MXd &centers, const pgo::EigenSupport::VXd &centerRadii, const int level);
  virtual void initForFigures(const pgo::EigenSupport::MXd &centers, const pgo::EigenSupport::VXd &centerRadii);
  void update(const pgo::EigenSupport::MXd &centers, const pgo::EigenSupport::VXd &centerRadii);

  double radius;
  int numCylinderSlices;
  int primitiveType;                  // 1 for sphere, 2 for cylinder, 3 for prism
  pgo::EigenSupport::MXd rayDirInit;  // normalized direction of the ray
  pgo::EigenSupport::MXd rayDir;      // normalized direction of the ray
  pgo::EigenSupport::MXd rayStartPosInit;
  pgo::EigenSupport::MXd rayStartPos;
  pgo::EigenSupport::VXd rayInitialLength;
  pgo::EigenSupport::VXd rayCurrLength;
  pgo::EigenSupport::MXd centersInit;
  pgo::EigenSupport::MXd centersCurr;
  // only used for prism
  pgo::EigenSupport::VXi rotAxisID;
  pgo::EigenSupport::MXd rayDirLocalFrameWeights;
  //
  pgo::EigenSupport::MXd baryCentricWeights;
  pgo::Mesh::TriMeshGeo primitiveTemplateMesh;
  pgo::Mesh::TriMeshPseudoNormal primitiveTemplateMeshNormal;
  std::vector<std::vector<int>> primitiveTemplateMeshVertexNeighboringVertices;

  int verbose = 0;
};
}  // namespace MedialAxisRepresentation
