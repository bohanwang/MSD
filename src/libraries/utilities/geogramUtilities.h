#pragma once

#include "triMeshGeo.h"

#include <string>

namespace MedialAxisRepresentation
{
namespace GeogramUtilities
{
void initGEO();
int remesh(const std::string &meshFilename, int targetNumPoints, double sizeFactor, double anisotropy, TriMeshGeo &meshOut);
}  // namespace GeogramUtilities
}  // namespace MedialAxisRepresentation