#include "geogramUtilities.h"

#include "logger.h"

#include <geogram/mesh/mesh.h>
#include <geogram/mesh/mesh_io.h>
#include <geogram/mesh/mesh_geometry.h>
#include <geogram/mesh/mesh_preprocessing.h>
#include <geogram/mesh/mesh_reorder.h>
#include <geogram/mesh/mesh_frame_field.h>
#include <geogram/mesh/mesh_tetrahedralize.h>
#include <geogram/mesh/mesh_decimate.h>
#include <geogram/mesh/mesh_remesh.h>
#include <geogram/basic/command_line_args.h>

#include <mutex>

namespace MedialAxisRepresentation
{
namespace GeogramUtilities
{
std::once_flag geo_init;
}
}  // namespace MedialAxisRepresentation

void MedialAxisRepresentation::GeogramUtilities::initGEO()
{
  std::call_once(geo_init, []() {
    GEO::initialize();
    GEO::CmdLine::import_arg_group("standard");
    GEO::CmdLine::import_arg_group("remesh");
    GEO::CmdLine::import_arg_group("algo");
  });
}

int MedialAxisRepresentation::GeogramUtilities::remesh(const std::string &meshFilename, int targetNumPoints, double sizeFactor, double anisotropy, TriMeshGeo &meshOut)
{
  GEO::Mesh inputMesh;
  GEO::Mesh remeshedMesh;
  // GEO::mesh_load(meshFilename, inputMesh);

  if (!GEO::mesh_load(meshFilename, inputMesh)) {
    SPDLOG_LOGGER_ERROR(Logger::lgr(), "Cannot open file {}", meshFilename);
    return 1;
  }

  if (sizeFactor > 0) {
    GEO::compute_sizing_field(inputMesh, sizeFactor);
  }

  if (anisotropy > 0) {
    GEO::compute_normals(inputMesh);
    GEO::simple_Laplacian_smooth(inputMesh, 3, true);  // true: smooth normals
    GEO::set_anisotropy(inputMesh, anisotropy * 0.02);
  }

  GEO::remesh_smooth(inputMesh, remeshedMesh, targetNumPoints);

  meshOut.clear();
  for (GEO::index_t vi = 0; vi < remeshedMesh.vertices.nb(); ++vi) {
    const double *p = remeshedMesh.vertices.point_ptr(vi);
    Vec3d pos(p);
    meshOut.addPos(pos);
  }

  for (GEO::index_t fi = 0; fi < remeshedMesh.facets.nb(); ++fi) {
    GEO::index_t nvtx = remeshedMesh.facets.nb_vertices(fi);
    ALOG(nvtx == 3);

    Vec3i tri(remeshedMesh.facets.vertex(fi, 0),
      remeshedMesh.facets.vertex(fi, 1),
      remeshedMesh.facets.vertex(fi, 2));

    meshOut.addTri(tri);
  }

  return 0;
}