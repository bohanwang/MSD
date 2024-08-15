#pragma once

#include "templatePrimitive.h"

#include "triMeshGeo.h"
#include "boundingVolumeTree.h"
#include "triMeshPseudoNormal.h"

#include <memory>

namespace MedialAxisRepresentation
{
struct PrimitiveICPSolverData;

enum class SolverType : int
{
  LINEARSOLVER = 0,
  PDSOLVER = 1
};

struct PrimitiveICPSolverParameters
{
  int maxNumIter = 20;
  int maxPDNumIter = 5;

  double sphereRadius = 0.01;  // change w/ minRadius of medial axis
  double expansionRadius = 0.1;
  // double smoothnessCoeff = 1.0; //1e4;
  double smoothnessCoeff = 1.0;  // 1e4;
  // double smoothnessCoeff = 1e4; //1e4;
  double smoothnessDecay = 1.0;
  // double expansionCoeff = 0.5;
  double expansionCoeff = 1.0;
  double contactCoeff = 1e3;
  // cylinder
  // int numCylinderSlices = 40; // for radius 0.01
  int numCylinderSlices = 10;  // for radius 0.01
  // SolverType solType = LINEARSOLVER;
  SolverType solType = SolverType::PDSOLVER;
  bool useCorefine = true;
  bool useEdgeOnlyMRF = true;
  bool verbose = false;
};

class PrimitiveICPSolver
{
public:
  PrimitiveICPSolver(const pgo::EigenSupport::MXd &centers, const pgo::EigenSupport::VXd &centerRadii, const PrimitiveICPSolverParameters &param_, int primitiveType_);
  ~PrimitiveICPSolver();

  int sim(const pgo::EigenSupport::MXd &centers, const pgo::Mesh::TriMeshGeo &targetMesh,
    const pgo::Mesh::TriMeshBVTree &targetMeshBVTree,
    const pgo::Mesh::TriMeshPseudoNormal &targetMeshNormals,
    const std::string &prefix, pgo::Mesh::TriMeshGeo &meshOut);

  // protected:
  PrimitiveICPSolverParameters params;
  int primitiveType;  // 1 for sphere, 2 for cylinder, 3 for prism
  std::shared_ptr<TemplatePrimitive> primitiveTemp;
  PrimitiveICPSolverData *pd;
};
}  // namespace MedialAxisRepresentation
