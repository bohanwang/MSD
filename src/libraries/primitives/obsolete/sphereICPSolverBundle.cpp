#include "sphereICPSolverBundle.h"

MedialAxisRepresentation::SphereICPSolverBundle::SphereICPSolverBundle(const TemplateSpheres &ts, const SphereICPSolverParameters &params):
  templateSpheres(ts)
{
  for (int si = 0; si < templateSpheres.getNumTemplateSpheres(); si++) {
    icpSolvers.emplace_back();
  }

  int sampleCount = 1;
  for (int si = (int)templateSpheres.getNumTemplateSpheres() - 1; si >= 0; si--) {
    SphereICPSolverParameters p = params;
    p.sampleCount = sampleCount;

    icpSolvers[si] = std::make_shared<SphereICPSolver>(templateSpheres.getSphereMesh(si), p);
    sampleCount = std::min(sampleCount + 1, 5);
  }
}