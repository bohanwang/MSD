#pragma once

#include "sphericalFunction/templateSpheres.h"

#include <memory>

namespace MedialAxisRepresentation
{
class SphereICPSolverBundle
{
public:
  SphereICPSolverBundle(const TemplateSpheres &templateSpheres, const SphereICPSolverParameters &params);

  std::shared_ptr<SphereICPSolver> getICPSolver(int sid) const { return icpSolvers[sid]; }
  std::shared_ptr<SphereICPSolver> getICPSolver(double edgeLength) const;
  const TemplateSpheres &getTemplateSpheres() const { return templateSpheres; }

protected:
  const TemplateSpheres &templateSpheres;
  std::vector<std::shared_ptr<SphereICPSolver>> icpSolvers;
};

inline std::shared_ptr<SphereICPSolver> SphereICPSolverBundle::getICPSolver(double edgeLength) const
{
  int sid = templateSpheres.getSphereMeshIDFromEdgeLength(edgeLength);
  return icpSolvers[sid];
}

}  // namespace MedialAxisRepresentation