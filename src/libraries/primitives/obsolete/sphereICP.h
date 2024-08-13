#pragma once

#include "triMeshGeo.h"
#include "boundingVolumeTree.h"
#include "triMeshPseudoNormal.h"

namespace MedialAxisRepresentation
{
enum class CenterOptimizationFlag : int
{
  FIXED = 0,
  OPTIM = 1,
  FOLLOW = 2
};

enum class ExpansionConstraintType : int
{
  LINEAR = 0,
  QUADRATIC = 1,
};

enum class NumPositionDOFs : int
{
  THREE = 3,
  TWELVE = 12,
};
struct SphereICPSolverParameters
{
  // smoothness
  double smooth_w_init = 10;
  double smooth_w_final = 1e-4;
  // icp/contact
  double icp_w = 1e2;
  // inertia
  double inertial_w = 0.0;
  // expansion
  double expansion_w = 0.1;
  // barrier
  double barrier_w = 1.0;
  // constrain rigid motion (no translation/rotation)
  double fixTranslation_w = 1e7;

  // max allow triangle scale
  double triangleScale = 5;
  // expansion step size
  double expansionRelStepSize = 0.1;
  double expansionRelMaxStepSize = 0.1;
  // timestep if doing dynamic (-1 means no dynamics)
  double timestep = -1.0;

  // # iterations
  int numIter = 20; 
  // # iterations that weight will change
  int numWeightIter = 12;

  // inner solver config
  int solverVerbose = 0;
  int solverNumIter = 30;

  // pnorm = 2: quadratic energy for distance/icp constraints
  double pnorm = 2.0;

  // use hard inequality constraints
  bool addConstraints = false;

  // whether resample triangles =1: no; >1: yes
  int sampleCount = 1;

  // use position as DOF to solve the problem or use affine transformation
  NumPositionDOFs numPositionDOFs = NumPositionDOFs::THREE;
  // whether the center is fixed or optimized
  CenterOptimizationFlag centerFlag = CenterOptimizationFlag::FIXED;
  // use quadratic expansion term
  ExpansionConstraintType expansionConstraintType = ExpansionConstraintType::QUADRATIC;
};

struct SphereICPSolverData;

class SphereICPSolver
{
public:
  SphereICPSolver(const TriMeshGeo &sphereRestMesh, const SphereICPSolverParameters &params);
  ~SphereICPSolver();

  int icp(const double center[3], double r,
    const TriMeshGeo &targetMesh, const TriMeshBVTree &targetMeshBVTree, const TriMeshPseudoNormal &targetMeshNormals,
    const std::string &prefix, TriMeshGeo &meshOut);

  int sim(const double center[3], double r,
    const TriMeshGeo &targetMesh, const TriMeshBVTree &targetMeshBVTree, const TriMeshPseudoNormal &targetMeshNormals,
    const std::string &prefix, TriMeshGeo &meshOut, Vec3d *updatedCenter = nullptr);

protected:
  TriMeshGeo sphereTemplateMesh;
  SphereICPSolverParameters params;
  SphereICPSolverData *sd;
};

int sphereICP(const TriMeshGeo &sphereMesh,
  const TriMeshGeo &targetMesh, const TriMeshBVTree &targetMeshBVTree, const TriMeshPseudoNormal &targetMeshNormals,
  double w_init, double w_final, int numIter, double timestep, const std::string &prefix, TriMeshGeo &meshOut);

}  // namespace MedialAxisRepresentation