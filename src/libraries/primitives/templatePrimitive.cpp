#include "templatePrimitive.h"
#include "templateSpheres.h"

#include "pgoLogging.h"
#include "createTriMesh.h"
#include "cgalInterface.h"
#include "libiglInterface.h"

#include <fmt/format.h>

namespace MedialAxisRepresentation
{
namespace ES = pgo::EigenSupport;

void rotatePoint(const ES::V3d &p, const ES::V3d &origin, const ES::V3d &rotAxism, double angle, ES::V3d &pRotated);

pgo::Mesh::TriMeshGeo createPrismWallMesh(const ES::V3d &e1, const ES::V3d &e2, const ES::V3d &e3, const double thickness, const double approxTargetRadius,
  std::vector<ES::V3d> &prismVtxDri, std::vector<ES::V3d> &baryCentricWeights, std::vector<int> &rotAxisID, double useHardCodedTargetEdgeLenForTri = -1);

// v1 = rot * v0
void rotAlignMat(const ES::V3d &v0, const ES::V3d &v1, ES::M3d &rot)
{
  PGO_ALOG(std::abs(v0.norm() - 1) < 1e-6);
  PGO_ALOG(std::abs(v1.norm() - 1) < 1e-6);

  // calculate rotation matrix
  ES::V3d v = v0.cross(v1);
  double s = v.norm();
  double c = v0.dot(v1);
  ES::M3d vx;
  vx << 0, -v(2), v(1),
    v(2), 0, -v(0),
    -v(1), v(0), 0;

  rot = ES::M3d::Identity();
  if (s > 1e-6) {
    rot += vx + vx * vx * (1 - c) / (s * s);
  }
  else {
    rot += vx;
  }
}

void getTemplateMeshForLevel(int level, int primitiveType, pgo::Mesh::TriMeshGeo &outMesh)
{
  if (primitiveType == 1) {
    if (level == -1) {
      TemplateSpheres templateSpheres(200, 3800, MedialAxisRepresentation::SphereRefinementMethod::ISOTROPIC, "template-spheres");
      outMesh = templateSpheres.sphereMeshes[3];  // radius 1.0, center (0, 0, 0)
    }
  }
  else if (primitiveType == 2) {
    if (level == -1) {
    }
  }
}

void TemplatePrimitive::init(const ES::MXd &centers, const ES::VXd &centerRadii, const int level)
{
  PGO_ALOG(centers.rows() == primitiveType);
  PGO_ALOG(centerRadii.size() == primitiveType);

  if (primitiveType == 1) {
    // sphere
    if (level == 0) {
      TemplateSpheres templateSpheres(500, 500, MedialAxisRepresentation::SphereRefinementMethod::ISOTROPIC, "template-spheres");
      primitiveTemplateMesh = templateSpheres.sphereMeshes[2];  // radius 1.0, center (0, 0, 0)
    }
    else if (level == 1) {
      TemplateSpheres templateSpheres(200, 200, MedialAxisRepresentation::SphereRefinementMethod::ISOTROPIC, "template-spheres");
      primitiveTemplateMesh = templateSpheres.sphereMeshes[1];  // radius 1.0, center (0, 0, 0)
    }
    else if (level == 2) {
      TemplateSpheres templateSpheres(50, 50, MedialAxisRepresentation::SphereRefinementMethod::ISOTROPIC, "template-spheres");
      primitiveTemplateMesh = templateSpheres.sphereMeshes[0];  // radius 1.0, center (0, 0, 0)
    }
    else {
      PGO_ALOG(false);
    }
    // TemplateSpheres templateSpheres(5, 6, MedialAxisRepresentation::SphereRefinementMethod::ISOTROPIC, "debug-spheres");
    // fmt::print("templateSpheres.sphereMeshes.size() = {}\n", templateSpheres.sphereMeshes.size());
    // primitiveTemplateMesh = templateSpheres.sphereMeshes[0];  // radius 1.0, center (0, 0, 0)

    int numVtx = primitiveTemplateMesh.numVertices();
    ES::V3d center = centers.row(0);

    rayDir.resize(numVtx, 3);
    rayStartPos.resize(numVtx, 3);
    for (int i = 0; i < numVtx; i++) {
      ES::V3d p = primitiveTemplateMesh.pos(i);
      rayDir.row(i) = p / p.norm();
      rayStartPos.row(i) = center;
    }

    rayInitialLength.resize(numVtx);
    // rayInitialLength.setConstant(centerRadii(0));
    rayInitialLength.setConstant(radius);
    for (int i = 0; i < numVtx; i++) {
      // primitiveTemplateMesh.pos(i) *= centerRadii(0);
      primitiveTemplateMesh.pos(i) *= radius;
      primitiveTemplateMesh.pos(i) += center;
    }

    baryCentricWeights.resize(numVtx, 1);
    baryCentricWeights.setConstant(1.0);
  }
  else if (primitiveType == 2) {
    // cylinder
    ES::V3d c1 = centers.row(0);
    ES::V3d c2 = centers.row(1);
    // double approxTargetRadius = (centerRadii(0) + centerRadii(1)) / 2;
    double approxTargetRadius = std::max(centerRadii(0), centerRadii(1));
    double height = (c1 - c2).norm();

    if (false) {
      fmt::print("approxTargetRadius = {}, height = {}, initRadius = {}\n", approxTargetRadius, height, radius);
    }

    // int numCirclePts = (int)(approxTargetRadius / 0.2 * numCylinderSlices);            // 10000
    // int numHeightPts = (int)(height / (2 * M_PI * approxTargetRadius) * numCirclePts);  // 1000

    int totalPt = -1;
    if (level == 0) {
      totalPt = 500;
    }
    else if (level == 1) {
      totalPt = 200;
    }
    else if (level == 2) {
      totalPt = 50;
    }
    else {
      PGO_ALOG(false);
    }

    double ratio = height / (2 * M_PI * approxTargetRadius);
    int numCirclePts = (int)std::sqrt(totalPt / ratio);
    int numHeightPts = (int)(ratio * numCirclePts);

    // debug
    // numCirclePts = 20;
    // numHeightPts = 30;

    primitiveTemplateMesh = pgo::Mesh::createCylinderWallMesh(radius, height, numCirclePts, numHeightPts);
    int numVtx = primitiveTemplateMesh.numVertices();
    fmt::print("numVtx = {}\n", numVtx);

    // transform cylinder to the given axis
    ES::V3d c1Hat(0, height / 2, 0), c2Hat(0, -height / 2, 0);
    ES::V3d translate = (c1 + c2) / 2;  // - origin
    c1Hat += translate;
    c2Hat += translate;
    ES::V3d axis = (c2 - c1).cross(c2Hat - c1Hat);
    axis = axis / axis.norm();
    double angle = acos((c2Hat - c1Hat).dot(c1 - c2) / ((c2Hat - c1Hat).norm() * (c1 - c2).norm()));

    // Use Rodrigues rotation formula
    ES::M3d W;
    W << 0, -axis(2), axis(1),
      axis(2), 0, -axis(0),
      -axis(1), axis(0), 0;
    ES::M3d R = ES::M3d::Identity() + sin(angle) * W + (1 - cos(angle)) * W * W;

    for (int i = 0; i < numVtx; i++) {
      ES::V3d v = primitiveTemplateMesh.pos(i);
      v = R * v + translate;
      primitiveTemplateMesh.pos(i) = v;
    }

    rayDir.resize(numVtx, 3);
    rayStartPos.resize(numVtx, 3);
    rayInitialLength.resize(numVtx);
    baryCentricWeights.resize(numVtx, 2);

    // ES::V3d cyAlis = centers.row(1) - centers.row(0);
    // cyAlis /= cyAlis.norm();

    for (int i = 0; i < numVtx; i++) {
      ES::V3d seg = c2 - c1;
      ES::V3d pos = primitiveTemplateMesh.pos(i);
      ES::V3d vec = pos - c1;
      double vTs_sTs = vec.dot(seg) / seg.dot(seg);
      ES::V3d dir = vec - vTs_sTs * seg;

      rayInitialLength(i) = dir.norm();
      rayDir.row(i) = dir / rayInitialLength(i);

      baryCentricWeights(i, 0) = 1 - vTs_sTs;
      baryCentricWeights(i, 1) = vTs_sTs;

      rayStartPos.row(i) = c1 + vTs_sTs * seg;
    }
  }
  else {
    // prism
    PGO_ALOG(primitiveType == 3);
    ES::V3d c0 = centers.row(0);
    ES::V3d c1 = centers.row(1);
    ES::V3d c2 = centers.row(2);
    double approxTargetRadius = (centerRadii(0) + centerRadii(1) + centerRadii(2)) / 2;
    double thickness = radius;

    std::vector<ES::V3d> prismVtxDri;
    std::vector<ES::V3d> baryCentricWeightsVec;
    std::vector<int> rotAxisIDVec;

    double useHardCodedTargetEdgeLenForTri = 0.1;
    if (level == 0) {
      useHardCodedTargetEdgeLenForTri = 0.05;
    }
    else if (level == 1) {
      useHardCodedTargetEdgeLenForTri = 0.1;
    }
    else if (level == 2) {
      useHardCodedTargetEdgeLenForTri = 0.5;
    }
    else {
      PGO_ALOG(false);
    }

    primitiveTemplateMesh = createPrismWallMesh(c0, c1, c2, thickness, approxTargetRadius, prismVtxDri, baryCentricWeightsVec, rotAxisIDVec, useHardCodedTargetEdgeLenForTri);

    PGO_ALOG(baryCentricWeightsVec.size() == rotAxisIDVec.size());

    int numVtx = prismVtxDri.size();
    fmt::print("numVtx = {}\n", numVtx);

    rayDir.resize(numVtx, 3);
    baryCentricWeights.resize(numVtx, 3);
    rayStartPos.resize(numVtx, 3);

    rayInitialLength.resize(numVtx);
    rayInitialLength.setConstant(thickness);

    rotAxisID.resize(numVtx);
    rotAxisID.setConstant(-1);

    rayDirLocalFrameWeights.resize(numVtx, 2);
    rayDirLocalFrameWeights.setZero();

    for (int i = 0; i < numVtx; i++) {
      ES::V3d rayDir_i = prismVtxDri[i];
      rayDir.row(i) = rayDir_i;
      PGO_ALOG(std::abs(rayDir_i.norm() - 1) < 1e-6);
      baryCentricWeights.row(i) = baryCentricWeightsVec[i];

      rayStartPos.row(i) = c0 * baryCentricWeightsVec[i](0) + c1 * baryCentricWeightsVec[i](1) + c2 * baryCentricWeightsVec[i](2);

      std::array<ES::V3d, 3> rayDirRotAxis = { ES::V3d(c1 - c0), ES::V3d(c2 - c1), ES::V3d(c0 - c2) };
      ES::V3d triPlaneNormal = rayDirRotAxis[0].cross(rayDirRotAxis[1]);
      triPlaneNormal /= triPlaneNormal.norm();

      // bool isOnEdge = false;
      // int rotAxisID_i = -1;
      // // 0 - 12, 1 - 20, 2 - 01; 01 - 0, 12 - 1, 20 - 2;
      // // 0 - 1, 1 - 2, 2 - 0
      // for (int k = 0; k < 3; k++) {
      //   if (isOnEdge) {
      //     PGO_ALOG(false);
      //   }
      //   if (baryCentricWeightsVec[i](k) < 1e-6) {
      //     isOnEdge = true;
      //     rotAxisID_i = (k + 1) % 3;
      //     break;
      //   }
      // }
      int rotAxisID_i = rotAxisIDVec[i];
      rotAxisID(i) = rotAxisID_i;

      if (rotAxisID_i == -1) {
        if (triPlaneNormal.dot(rayDir_i) < 0) {
          rayDirLocalFrameWeights(i, 0) = -1;
        }
        else {
          rayDirLocalFrameWeights(i, 0) = 1;
        }
      }
      else {
        ES::V3d rayDirRotAxis_i = rayDirRotAxis[rotAxisID_i];
        rayDirRotAxis_i /= rayDirRotAxis_i.norm();
        ES::V3d axisZ_i = triPlaneNormal.cross(rayDirRotAxis_i);
        axisZ_i /= axisZ_i.norm();

        double debugDot = rayDir_i.dot(rayDirRotAxis_i);
        PGO_ALOG(std::abs(rayDir_i.dot(rayDirRotAxis_i)) < 1e-6);

        rayDirLocalFrameWeights(i, 0) = rayDir_i.dot(triPlaneNormal);
        rayDirLocalFrameWeights(i, 1) = rayDir_i.dot(axisZ_i);

        ES::V3d debugDir = rayDirLocalFrameWeights(i, 0) * triPlaneNormal + rayDirLocalFrameWeights(i, 1) * axisZ_i;
        PGO_ALOG((debugDir - rayDir_i).norm() < 1e-6);
      }
    }
  }
  centersInit = centers;
  rayDirInit = rayDir;
  rayStartPosInit = rayStartPos;
  rayCurrLength = rayInitialLength;

  primitiveTemplateMeshNormal.buildPseudoNormals(primitiveTemplateMesh);
  // primitiveTemplateMesh.save("t.obj");

  primitiveTemplateMeshNormal.updateVertexPositions(primitiveTemplateMesh);

  pgo::Mesh::TriMeshNeighbor meshNeighbor(primitiveTemplateMesh);
  primitiveTemplateMeshVertexNeighboringVertices.resize(primitiveTemplateMesh.numVertices());
  for (int vi = 0; vi < (int)primitiveTemplateMeshVertexNeighboringVertices.size(); vi++) {
    primitiveTemplateMeshVertexNeighboringVertices[vi] = meshNeighbor.getVtxNearbyVertices(vi, primitiveTemplateMesh);
  }
}

void TemplatePrimitive::init(const ES::MXd &centers, const ES::VXd &centerRadii, const double useHardCodedTargetEdgeLenForTri)
{
  PGO_ALOG(centers.rows() == primitiveType);
  PGO_ALOG(centerRadii.size() == primitiveType);

  if (primitiveType == 1) {
    // sphere
    TemplateSpheres templateSpheres(200, 3800, MedialAxisRepresentation::SphereRefinementMethod::ISOTROPIC, "template-spheres");
    // // PGO_ALOG(templateSpheres.sphereMeshes.size() == 4);
    primitiveTemplateMesh = templateSpheres.sphereMeshes[3];  // radius 1.0, center (0, 0, 0)
    // TemplateSpheres templateSpheres(5, 6, MedialAxisRepresentation::SphereRefinementMethod::ISOTROPIC, "debug-spheres");
    fmt::print("templateSpheres.sphereMeshes.size() = {}\n", templateSpheres.sphereMeshes.size());
    // primitiveTemplateMesh = templateSpheres.sphereMeshes[0];  // radius 1.0, center (0, 0, 0)

    int numVtx = primitiveTemplateMesh.numVertices();
    ES::V3d center = centers.row(0);

    rayDir.resize(numVtx, 3);
    rayStartPos.resize(numVtx, 3);
    for (int i = 0; i < numVtx; i++) {
      ES::V3d p = primitiveTemplateMesh.pos(i);
      rayDir.row(i) = p / p.norm();
      rayStartPos.row(i) = center;
    }

    rayInitialLength.resize(numVtx);
    // rayInitialLength.setConstant(centerRadii(0));
    rayInitialLength.setConstant(radius);
    for (int i = 0; i < numVtx; i++) {
      // primitiveTemplateMesh.pos(i) *= centerRadii(0);
      primitiveTemplateMesh.pos(i) *= radius;
      primitiveTemplateMesh.pos(i) += center;
    }

    baryCentricWeights.resize(numVtx, 1);
    baryCentricWeights.setConstant(1.0);
  }
  else if (primitiveType == 2) {
    // cylinder
    ES::V3d c1 = centers.row(0);
    ES::V3d c2 = centers.row(1);
    // double approxTargetRadius = (centerRadii(0) + centerRadii(1)) / 2;
    double approxTargetRadius = std::max(centerRadii(0), centerRadii(1));
    double height = (c1 - c2).norm();

    if (false) {
      fmt::print("approxTargetRadius = {}, height = {}, initRadius = {}\n", approxTargetRadius, height, radius);
    }

    // int numCirclePts = (int)(approxTargetRadius / 0.2 * numCylinderSlices);            // 10000
    // int numHeightPts = (int)(height / (2 * M_PI * approxTargetRadius) * numCirclePts);  // 1000

    double ratio = height / (2 * M_PI * approxTargetRadius);
    int totalPt = 1500;
    int numCirclePts = (int)std::sqrt(totalPt / ratio);
    int numHeightPts = (int)(ratio * numCirclePts);

    // debug
    // numCirclePts = 20;
    // numHeightPts = 30;

    primitiveTemplateMesh = pgo::Mesh::createCylinderWallMesh(radius, height, numCirclePts, numHeightPts);
    int numVtx = primitiveTemplateMesh.numVertices();
    fmt::print("numVtx = {}\n", numVtx);

    // transform cylinder to the given axis
    ES::V3d c1Hat(0, height / 2, 0), c2Hat(0, -height / 2, 0);
    ES::V3d translate = (c1 + c2) / 2;  // - origin
    c1Hat += translate;
    c2Hat += translate;
    ES::V3d axis = (c2 - c1).cross(c2Hat - c1Hat);
    axis = axis / axis.norm();
    double angle = acos((c2Hat - c1Hat).dot(c1 - c2) / ((c2Hat - c1Hat).norm() * (c1 - c2).norm()));
    // Use Rodrigues rotation formula
    ES::M3d W;
    W << 0, -axis(2), axis(1),
      axis(2), 0, -axis(0),
      -axis(1), axis(0), 0;
    ES::M3d R = ES::M3d::Identity() + sin(angle) * W + (1 - cos(angle)) * W * W;

    for (int i = 0; i < numVtx; i++) {
      ES::V3d v = primitiveTemplateMesh.pos(i);
      v = R * v + translate;
      primitiveTemplateMesh.pos(i) = v;
    }

    rayDir.resize(numVtx, 3);
    rayStartPos.resize(numVtx, 3);
    rayInitialLength.resize(numVtx);
    baryCentricWeights.resize(numVtx, 2);

    // ES::V3d cyAlis = centers.row(1) - centers.row(0);
    // cyAlis /= cyAlis.norm();

    for (int i = 0; i < numVtx; i++) {
      ES::V3d seg = c2 - c1;
      ES::V3d pos = primitiveTemplateMesh.pos(i);
      ES::V3d vec = pos - c1;
      double vTs_sTs = vec.dot(seg) / seg.dot(seg);
      ES::V3d dir = vec - vTs_sTs * seg;

      rayInitialLength(i) = dir.norm();
      rayDir.row(i) = dir / rayInitialLength(i);

      baryCentricWeights(i, 0) = 1 - vTs_sTs;
      baryCentricWeights(i, 1) = vTs_sTs;

      rayStartPos.row(i) = c1 + vTs_sTs * seg;
    }
  }
  else {
    // prism
    PGO_ALOG(primitiveType == 3);
    ES::V3d c0 = centers.row(0);
    ES::V3d c1 = centers.row(1);
    ES::V3d c2 = centers.row(2);
    double approxTargetRadius = (centerRadii(0) + centerRadii(1) + centerRadii(2)) / 2;
    double thickness = radius;

    std::vector<ES::V3d> prismVtxDri;
    std::vector<ES::V3d> baryCentricWeightsVec;
    std::vector<int> rotAxisIDVec;
    primitiveTemplateMesh = createPrismWallMesh(c0, c1, c2, thickness, approxTargetRadius, prismVtxDri, baryCentricWeightsVec, rotAxisIDVec, useHardCodedTargetEdgeLenForTri);

    PGO_ALOG(baryCentricWeightsVec.size() == rotAxisIDVec.size());

    int numVtx = prismVtxDri.size();

    rayDir.resize(numVtx, 3);
    baryCentricWeights.resize(numVtx, 3);
    rayStartPos.resize(numVtx, 3);

    rayInitialLength.resize(numVtx);
    rayInitialLength.setConstant(thickness);

    rotAxisID.resize(numVtx);
    rotAxisID.setConstant(-1);

    rayDirLocalFrameWeights.resize(numVtx, 2);
    rayDirLocalFrameWeights.setZero();

    for (int i = 0; i < numVtx; i++) {
      ES::V3d rayDir_i = prismVtxDri[i];
      rayDir.row(i) = rayDir_i;
      PGO_ALOG(std::abs(rayDir_i.norm() - 1) < 1e-6);
      baryCentricWeights.row(i) = baryCentricWeightsVec[i];

      rayStartPos.row(i) = c0 * baryCentricWeightsVec[i](0) + c1 * baryCentricWeightsVec[i](1) + c2 * baryCentricWeightsVec[i](2);

      std::array<ES::V3d, 3> rayDirRotAxis = { ES::V3d(c1 - c0), ES::V3d(c2 - c1), ES::V3d(c0 - c2) };
      ES::V3d triPlaneNormal = rayDirRotAxis[0].cross(rayDirRotAxis[1]);
      triPlaneNormal /= triPlaneNormal.norm();

      // bool isOnEdge = false;
      // int rotAxisID_i = -1;
      // // 0 - 12, 1 - 20, 2 - 01; 01 - 0, 12 - 1, 20 - 2;
      // // 0 - 1, 1 - 2, 2 - 0
      // for (int k = 0; k < 3; k++) {
      //   if (isOnEdge) {
      //     PGO_ALOG(false);
      //   }
      //   if (baryCentricWeightsVec[i](k) < 1e-6) {
      //     isOnEdge = true;
      //     rotAxisID_i = (k + 1) % 3;
      //     break;
      //   }
      // }
      int rotAxisID_i = rotAxisIDVec[i];
      rotAxisID(i) = rotAxisID_i;

      if (rotAxisID_i == -1) {
        if (triPlaneNormal.dot(rayDir_i) < 0) {
          rayDirLocalFrameWeights(i, 0) = -1;
        }
        else {
          rayDirLocalFrameWeights(i, 0) = 1;
        }
      }
      else {
        ES::V3d rayDirRotAxis_i = rayDirRotAxis[rotAxisID_i];
        rayDirRotAxis_i /= rayDirRotAxis_i.norm();
        ES::V3d axisZ_i = triPlaneNormal.cross(rayDirRotAxis_i);
        axisZ_i /= axisZ_i.norm();

        double debugDot = rayDir_i.dot(rayDirRotAxis_i);
        PGO_ALOG(std::abs(rayDir_i.dot(rayDirRotAxis_i)) < 1e-6);

        rayDirLocalFrameWeights(i, 0) = rayDir_i.dot(triPlaneNormal);
        rayDirLocalFrameWeights(i, 1) = rayDir_i.dot(axisZ_i);

        ES::V3d debugDir = rayDirLocalFrameWeights(i, 0) * triPlaneNormal + rayDirLocalFrameWeights(i, 1) * axisZ_i;
        PGO_ALOG((debugDir - rayDir_i).norm() < 1e-6);
      }
    }
  }
  centersInit = centers;
  rayDirInit = rayDir;
  rayStartPosInit = rayStartPos;
  rayCurrLength = rayInitialLength;

  primitiveTemplateMeshNormal.buildPseudoNormals(primitiveTemplateMesh);
  // primitiveTemplateMesh.save("t.obj");

  primitiveTemplateMeshNormal.updateVertexPositions(primitiveTemplateMesh);

  pgo::Mesh::TriMeshNeighbor meshNeighbor(primitiveTemplateMesh);
  primitiveTemplateMeshVertexNeighboringVertices.resize(primitiveTemplateMesh.numVertices());
  for (int vi = 0; vi < (int)primitiveTemplateMeshVertexNeighboringVertices.size(); vi++) {
    primitiveTemplateMeshVertexNeighboringVertices[vi] = meshNeighbor.getVtxNearbyVertices(vi, primitiveTemplateMesh);
  }
}

void TemplatePrimitive::update(const ES::MXd &centers, const ES::VXd &rayRadii)
{
  int numVtx = primitiveTemplateMesh.numVertices();
  PGO_ALOG(rayRadii.size() == numVtx);
  PGO_ALOG(centers.rows() == primitiveType);

  centersCurr = centers;
  rayCurrLength = rayRadii;

  // update rayDir for cylinders and prisims
  // update Mesh
  if (primitiveType == 1) {
    ES::V3d maVtxPos = centers.row(0);
    for (int i = 0; i < numVtx; i++) {
      ES::V3d rayDir_i = rayDir.row(i);
      ES::V3d posNew = maVtxPos + rayDir_i * rayRadii(i);
      primitiveTemplateMesh.pos(i) = posNew;
      rayStartPos.row(i) = maVtxPos;
    }
  }
  else if (primitiveType == 2) {
    ES::V3d cylAxisInit = centersInit.row(1) - centersInit.row(0);
    ES::V3d cylAxisInitNorm = cylAxisInit / cylAxisInit.norm();

    ES::V3d cylAxis = centers.row(1) - centers.row(0);
    ES::V3d cylAxisNorm = cylAxis / cylAxis.norm();

    ES::M3d rotMat;
    rotAlignMat(cylAxisInitNorm, cylAxisNorm, rotMat);

    for (int i = 0; i < numVtx; i++) {
      ES::V3d rayDirInit_i = rayDirInit.row(i);
      ES::V3d rayDir_i = rotMat * rayDirInit_i;
      PGO_ALOG(std::abs(rayDir_i.norm() - 1) < 1e-6);
      // rayDir_i /= rayDir_i.norm();
      rayDir.row(i) = rayDir_i;

      double w0 = baryCentricWeights(i, 0), w1 = baryCentricWeights(i, 1);
      ES::V3d centerNew = w0 * centers.row(0) + w1 * centers.row(1);
      ES::V3d posNew = centerNew + rayDir_i * rayRadii(i);
      primitiveTemplateMesh.pos(i) = posNew;
      rayStartPos.row(i) = centerNew;
    }
  }
  else if (primitiveType == 3) {
    ES::V3d maVtxPos0 = centers.row(0), maVtxPos1 = centers.row(1), maVtxPos2 = centers.row(2);
    std::array<ES::V3d, 3> rayDirRotAxis = { ES::V3d(maVtxPos1 - maVtxPos0), ES::V3d(maVtxPos2 - maVtxPos1), ES::V3d(maVtxPos0 - maVtxPos2) };
    ES::V3d triPlaneNormal = rayDirRotAxis[0].cross(rayDirRotAxis[1]);
    triPlaneNormal /= triPlaneNormal.norm();

    for (int i = 0; i < numVtx; i++) {
      int rotAxisID_i = rotAxisID(i);

      ES::V3d rayDir_i;
      ES::V2d rayDirLocalFrameWeights_i = rayDirLocalFrameWeights.row(i);
      if (rotAxisID_i == -1) {
        rayDir_i = triPlaneNormal * rayDirLocalFrameWeights_i(0);
      }
      else {
        ES::V3d rayDirRotAxis_i = rayDirRotAxis[rotAxisID_i];
        rayDirRotAxis_i /= rayDirRotAxis_i.norm();
        ES::V3d axisZ_i = triPlaneNormal.cross(rayDirRotAxis_i);
        axisZ_i /= axisZ_i.norm();

        rayDir_i = triPlaneNormal * rayDirLocalFrameWeights_i(0) + axisZ_i * rayDirLocalFrameWeights_i(1);
      }

      PGO_ALOG(std::abs(rayDir_i.norm() - 1) < 1e-6);
      rayDir.row(i) = rayDir_i;

      double w0 = baryCentricWeights(i, 0), w1 = baryCentricWeights(i, 1), w2 = baryCentricWeights(i, 2);
      ES::V3d centerNew = w0 * centers.row(0) + w1 * centers.row(1) + w2 * centers.row(2);
      ES::V3d posNew = centerNew + rayDir_i * rayRadii(i);
      primitiveTemplateMesh.pos(i) = posNew;
      rayStartPos.row(i) = centerNew;
    }
  }
}

pgo::Mesh::TriMeshGeo createPrismWallMesh(const ES::V3d &e1, const ES::V3d &e2, const ES::V3d &e3, const double thickness, const double approxTargetRadius,
  std::vector<ES::V3d> &prismVtxDri, std::vector<ES::V3d> &baryCentricWeights, std::vector<int> &rotAxisID, double useHardCodedTargetEdgeLenForTri)
{
  pgo::Mesh::TriMeshGeo prismMesh;  // final output
  const std::string prefix = "./";

  pgo::Mesh::TriMeshGeo triPlane;
  triPlane.addPos(e1);
  triPlane.addPos(e2);
  triPlane.addPos(e3);
  triPlane.addTri(ES::V3i(0, 1, 2));

  double minEdgeLength = std::min({ (e1 - e2).norm(), (e2 - e3).norm(), (e3 - e1).norm() });
  double maxEdgeLength = std::max({ (e1 - e2).norm(), (e2 - e3).norm(), (e3 - e1).norm() });
  double minmaxRatio = maxEdgeLength / minEdgeLength;
  double targetEdgeSegs = std::max(std::sqrt(200 / minmaxRatio), 2.0);
  double targetEdgeLen = minEdgeLength / targetEdgeSegs;

  if (useHardCodedTargetEdgeLenForTri > 0) {
    targetEdgeLen = useHardCodedTargetEdgeLenForTri;
  }

  pgo::Mesh::TriMeshGeo triPlaneRemesh = pgo::CGALInterface::isotropicRemeshing(triPlane, targetEdgeLen, 3, 90);
  // triPlaneRemesh.save(fmt::format("{}/triPlaneRemesh.obj", prefix));

  int numHeight = std::max(int(approxTargetRadius * M_PI / targetEdgeLen), 10);  // 0.01

  if (false) {
    fmt::print("numHeight = {}\n", numHeight);
  }

  double da = M_PI / numHeight;

  int numPlanePts = triPlaneRemesh.numVertices();
  std::vector<int> boundaryVtxID;
  pgo::libiglInterface::boundaryLoop(triPlaneRemesh, boundaryVtxID);

  // upper and lower plane
  ES::V3d upperDir = (e2 - e1).cross(e3 - e1);
  upperDir = upperDir / upperDir.norm();

  pgo::Mesh::TriMeshGeo upperPlane = triPlaneRemesh;
  pgo::Mesh::TriMeshGeo lowerPlane = triPlaneRemesh;
  std::vector<ES::V3d> baryCentricWeightsPlane;
  for (int i = 0; i < upperPlane.numVertices(); i++) {
    ES::V3d vPos = upperPlane.pos(i);
    ES::V3d vUp = vPos + upperDir * thickness;
    ES::V3d vLow = vPos - upperDir * thickness;
    upperPlane.pos(i) = vUp;
    lowerPlane.pos(i) = vLow;

    // calculate barycentric weights
    ES::V3d v0 = e2 - e1, v1 = e3 - e1, v2 = vPos - e1;
    double d00 = v0.dot(v0), d01 = v0.dot(v1), d11 = v1.dot(v1), d20 = v2.dot(v0), d21 = v2.dot(v1);
    double denom = d00 * d11 - d01 * d01;
    double v = (d11 * d20 - d01 * d21) / denom;
    double w = (d00 * d21 - d01 * d20) / denom;
    double u = 1 - v - w;
    baryCentricWeightsPlane.push_back(ES::V3d(u, v, w));
  }

  // reverse the triangle order
  for (int t = 0; t < lowerPlane.numTriangles(); t++) {
    ES::V3i tri = lowerPlane.tri(t);
    lowerPlane.tri(t) = ES::V3i(tri[2], tri[1], tri[0]);
  }
  // upperPlane.save(fmt::format("{}/upperPlane.obj", prefix));
  // lowerPlane.save(fmt::format("{}/lowerPlane.obj", prefix));

  prismMesh = upperPlane;
  prismMesh.addMesh(lowerPlane);

  // side faces
  int e1ID = -1, e2ID = -1, e3ID = -1;
  // align to e1
  for (int i = 0; i < (int)boundaryVtxID.size(); i++) {
    int id = boundaryVtxID[i];
    if ((triPlaneRemesh.pos(id) - e1).norm() < 1e-6) {
      e1ID = id;
    }
    if ((triPlaneRemesh.pos(id) - e2).norm() < 1e-6) {
      e2ID = id;
    }
    if ((triPlaneRemesh.pos(id) - e3).norm() < 1e-6) {
      e3ID = id;
    }
  }
  PGO_ALOG(e1ID != -1 && e2ID != -1 && e3ID != -1);
  PGO_ALOG(numPlanePts * 2 == prismMesh.numVertices());

  for (int i = 0; i < numPlanePts; i++) {
    prismVtxDri.push_back(upperDir);
  }
  for (int i = 0; i < numPlanePts; i++) {
    prismVtxDri.push_back(-1 * upperDir);
  }

  prismMesh.save("r.obj");

  baryCentricWeights.insert(baryCentricWeights.end(), baryCentricWeightsPlane.begin(), baryCentricWeightsPlane.end());
  baryCentricWeights.insert(baryCentricWeights.end(), baryCentricWeightsPlane.begin(), baryCentricWeightsPlane.end());

  rotAxisID.insert(rotAxisID.end(), numPlanePts, -1);
  rotAxisID.insert(rotAxisID.end(), numPlanePts, -1);

  // rotate boundaryVtxID to align to e1
  std::rotate(boundaryVtxID.begin(), std::find(boundaryVtxID.begin(), boundaryVtxID.end(), e1ID), boundaryVtxID.end());

  std::vector<int> startIDs = { e1ID, e2ID, e3ID };
  std::vector<ES::V3d> upperTriVtx = { e1, e2, e3 };
  std::vector<ES::V3d> triVtx = { e1, e2, e3 };
  for (int i = 0; i < 3; i++) {
    upperTriVtx[i] += upperDir * thickness;
  }
  int startBndID = 0;
  for (int id_ = 0; id_ < 3; id_++) {
    int startID = startIDs[id_];
    ES::V3d triVtxPrev = upperTriVtx[(id_ + 2) % 3], triVtxCurr = upperTriVtx[id_], triVtxNext = upperTriVtx[(id_ + 1) % 3];
    std::vector<ES::V3d> rot2Axis = { (triVtxCurr - triVtxPrev) / (triVtxCurr - triVtxPrev).norm(), (triVtxNext - triVtxCurr) / (triVtxNext - triVtxCurr).norm() };

    // Generate cylinder wall
    int startBndIDForTri = startBndID;
    int numCyWallVtx = 0;
    int startCyWallVtxID = prismMesh.numVertices();
    for (int bndId_ = startBndID; bndId_ <= (int)boundaryVtxID.size(); bndId_++) {
      int bndID = boundaryVtxID[bndId_ % boundaryVtxID.size()];
      int nextID = startIDs[(id_ + 1) % 3];

      int upperID = bndID;
      ES::V3d upperStart = upperPlane.pos(upperID);
      ES::V3d upperNext;
      ES::V3d rotAxis = rot2Axis[1];
      for (int i = 0; i < numHeight - 1; i++) {
        ES::V3d rotPt;
        ES::V3d startPt = upperStart;
        ES::V3d center = triPlaneRemesh.pos(upperID);
        rotatePoint(startPt, center, rotAxis, da, rotPt);
        upperNext = rotPt;
        prismMesh.addPos(upperNext);

        ES::V3d dirNewPos = upperNext - center;
        prismVtxDri.push_back(dirNewPos / dirNewPos.norm());

        ES::V3d baryCentricWeightCurr;
        baryCentricWeightCurr[(id_ + 2) % 3] = 0;
        baryCentricWeightCurr[id_] = 1 - (center - triVtx[id_]).norm() / (triVtx[(id_ + 1) % 3] - triVtx[id_]).norm();
        baryCentricWeightCurr[(id_ + 1) % 3] = 1 - baryCentricWeightCurr[id_];
        baryCentricWeights.push_back(baryCentricWeightCurr);
        rotAxisID.push_back(((id_ + 2) % 3 + 1) % 3);

        upperStart = upperNext;
      }
      numCyWallVtx++;

      if (bndID == nextID) {
        startBndID = bndId_;
        break;
      }
    }
    // std::cout << "numCyWallVtx: " << numCyWallVtx << std::endl;

    for (int n = 0; n < numCyWallVtx - 1; n++) {
      int upperID0 = boundaryVtxID[(startBndIDForTri + n) % boundaryVtxID.size()];
      int upperID1 = boundaryVtxID[(startBndIDForTri + n + 1) % boundaryVtxID.size()];
      // upper part
      int ID0 = startCyWallVtxID + n * (numHeight - 1);
      int ID1 = startCyWallVtxID + (n + 1) * (numHeight - 1);
      prismMesh.addTri(ES::V3i(upperID0, ID0, upperID1));
      prismMesh.addTri(ES::V3i(upperID1, ID0, ID1));

      // middle
      for (int i = 0; i < numHeight - 2; i++) {
        int IDi0 = startCyWallVtxID + n * (numHeight - 1) + i;
        int IDi1 = startCyWallVtxID + n * (numHeight - 1) + i + 1;
        int IDi2 = startCyWallVtxID + (n + 1) * (numHeight - 1) + i;
        int IDi3 = startCyWallVtxID + (n + 1) * (numHeight - 1) + i + 1;
        prismMesh.addTri(ES::V3i(IDi0, IDi1, IDi2));
        prismMesh.addTri(ES::V3i(IDi1, IDi3, IDi2));
      }

      // lower part
      int lowerID0 = upperID0 + numPlanePts;
      int lowerID1 = upperID1 + numPlanePts;
      ID0 = startCyWallVtxID + n * (numHeight - 1) + numHeight - 2;
      ID1 = startCyWallVtxID + (n + 1) * (numHeight - 1) + numHeight - 2;
      prismMesh.addTri(ES::V3i(ID0, lowerID0, ID1));
      prismMesh.addTri(ES::V3i(lowerID0, lowerID1, ID1));
    }
    // prismMesh.save(fmt::format("{}/prismMesh.obj", prefix));
  }

  PGO_ALOG(prismMesh.numVertices() == (int)prismVtxDri.size());
  PGO_ALOG(prismMesh.numVertices() == (int)baryCentricWeights.size());
  return prismMesh;
}

void rotatePoint(const ES::V3d &p, const ES::V3d &origin, const ES::V3d &rotAxism, double angle, ES::V3d &pRotated)
{
  PGO_ALOG(std::abs(rotAxism.norm() - 1) < 1e-6);
  Eigen::Quaternion<double> q(std::cos(angle / 2), std::sin(angle / 2) * rotAxism[0], std::sin(angle / 2) * rotAxism[1], std::sin(angle / 2) * rotAxism[2]);
  Eigen::Quaterniond pQuatern;
  pQuatern.w() = 0;
  pQuatern.vec() = p - origin;
  Eigen::Quaterniond pRotQuatern = q * pQuatern * q.inverse();
  pRotated = pRotQuatern.vec() + origin;
}

}  // namespace MedialAxisRepresentation
