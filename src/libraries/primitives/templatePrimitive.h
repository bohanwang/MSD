#pragma once

#include "CGALUtilities.h"
#include "libiglInterface.h"
#include "createTriMesh.h"
#include "triMeshGeo.h"
#include "triMeshNeighbor.h"
#include "boundingVolumeTree.h"
#include "triMeshPseudoNormal.h"

namespace MedialAxisRepresentation
{
void rotAlignMat(const ES::V3d &v0, const ES::V3d &v1, ES::M3d &rot);
TriMeshGeo createPrismWallMeshForFigures(const ES::V3d &e1, const ES::V3d &e2, const ES::V3d &e3, const double thickness, const double approxTargetRadius,
  std::vector<ES::V3d> &prismVtxDri, std::vector<ES::V3d> &baryCentricWeights, std::vector<int> &rotAxisID);
class TemplatePrimitive
{
public:
  TemplatePrimitive(){};
  TemplatePrimitive(double radius_, int numCylinderSlices_, int primitiveType_):
    radius(radius_), numCylinderSlices(numCylinderSlices_), primitiveType(primitiveType_){};
  ~TemplatePrimitive(){};

  virtual void init(const ES::MXd &centers, const ES::VXd &centerRadii, const double useHardCodedTargetEdgeLenForTri = -1);
  virtual void init(const ES::MXd &centers, const ES::VXd &centerRadii, const int level);
  virtual void initForFigures(const ES::MXd &centers, const ES::VXd &centerRadii);
  void update(const ES::MXd &centers, const ES::VXd &centerRadii);

  double radius;
  int numCylinderSlices;
  int primitiveType;   // 1 for sphere, 2 for cylinder, 3 for prism
  ES::MXd rayDirInit;  // normalized direction of the ray
  ES::MXd rayDir;      // normalized direction of the ray
  ES::MXd rayStartPosInit;
  ES::MXd rayStartPos;
  ES::VXd rayInitialLength;
  ES::VXd rayCurrLength;
  ES::MXd centersInit;
  ES::MXd centersCurr;
  // only used for prism
  ES::VXi rotAxisID;
  ES::MXd rayDirLocalFrameWeights;
  //
  ES::MXd baryCentricWeights;
  TriMeshGeo primitiveTemplateMesh;
  TriMeshPseudoNormal primitiveTemplateMeshNormal;
  std::vector<std::vector<int>> primitiveTemplateMeshVertexNeighboringVertices;

  int verbose = 0;
};
}  // namespace MedialAxisRepresentation
