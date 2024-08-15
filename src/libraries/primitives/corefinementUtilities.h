#pragma once

#include "cgalBasic.h"
#include "cgalTemplateUtilities.h"

#include "EigenSupport.h"
#include "createTriMesh.h"

#include <CGAL/Polygon_mesh_processing/orient_polygon_soup.h>
#include <CGAL/Polygon_mesh_processing/triangulate_faces.h>
#include <CGAL/Polygon_mesh_processing/corefinement.h>
#include <CGAL/Polygon_mesh_processing/remesh.h>

#include <CGAL/AABB_triangle_primitive.h>
#include <CGAL/Barycentric_coordinates_2/triangle_coordinates_2.h>
#include <CGAL/Polygon_mesh_processing/orient_polygon_soup_extension.h>

#include <CGAL/Cartesian_converter.h>
#include <CGAL/gmpxx.h>

namespace MedialAxisRepresentation
{
namespace ES = pgo::EigenSupport;

// using K = CGALUtilities::KernelExact;
using K = CGAL::Simple_cartesian<mpq_class>;
using Poly = pgo::CGALInterface::Polyhedron<K>;

using Triangle = K::Triangle_3;
using Iterator = std::vector<Triangle>::iterator;
using TrianaglePrimitive = CGAL::AABB_triangle_primitive<K, Iterator>;
using FaceGraphPrimitive = CGAL::AABB_face_graph_triangle_primitive<Poly>;
using AABB_triangle_traits = CGAL::AABB_traits<K, TrianaglePrimitive>;
using AABB_face_graph_traits = CGAL::AABB_traits<K, FaceGraphPrimitive>;
using TriangleTree = CGAL::AABB_tree<AABB_triangle_traits>;
using FaceGraphTree = CGAL::AABB_tree<AABB_face_graph_traits>;

using HClock = std::chrono::high_resolution_clock;
using HClockPt = HClock::time_point;

using real = CGAL::Gmpfr;
using V2 = Eigen::Matrix<real, 2, 1>;
using V3 = Eigen::Matrix<real, 3, 1>;
using IK = CGAL::Simple_cartesian<real>;

template<typename T>
concept SupportedRealType = std::is_same<std::decay_t<T>, double>::value ||
  std::is_same<std::decay_t<T>, CGAL::Gmpfr>::value;

using K_to_IK = CGAL::Cartesian_converter<K, IK>;
using IK_to_K = CGAL::Cartesian_converter<IK, K>;

template<typename T>
  requires SupportedRealType<T>
inline T EK_to_IK(const K::FT &v)
{
  if constexpr (std::is_same<T, double>::value) {
    return CGAL::to_double(v);
  }
  else if constexpr (std::is_same<T, CGAL::Gmpfr>::value) {
    T vIK;
    mpfr_set_q(vIK.fr(), v.get_mpq_t(), MPFR_RNDF);
    return vIK;
  }
}

template<typename T>
  requires SupportedRealType<T>
inline K::FT IK_to_EK(const T &v)
{
  if constexpr (std::is_same<T, double>::value) {
    return K::FT(v);
  }
  else if constexpr (std::is_same<T, CGAL::Gmpfr>::value) {
    K::FT vK;
    mpfr_get_q(vK.get_mpq_t(), v.fr());
    return vK;
  }
}

template<typename T>
  requires SupportedRealType<T> || std::is_same<T, K::FT>::value
inline double tod(const T &v)
{
  if constexpr (std::is_same<T, double>::value) {
    return v;
  }
  else if constexpr (std::is_same<T, CGAL::Gmpfr>::value) {
    return v.to_double();
  }
  else if constexpr (std::is_same<T, K::FT>::value) {
    return CGAL::to_double(v);
  }
}

inline double duraSecond(const HClockPt &t0, const HClockPt &t1)
{
  return double(std::chrono::duration_cast<std::chrono::microseconds>(t1 - t0).count()) * 1e-6;
}

inline ES::V3d toVec3(const K::Point_3 &pt)
{
  return ES::V3d(CGAL::to_double(pt.x()), CGAL::to_double(pt.y()), CGAL::to_double(pt.z()));
}

inline ES::V3d toVec3(const K::Vector_3 &pt)
{
  return ES::V3d(CGAL::to_double(pt[0]), CGAL::to_double(pt[1]), CGAL::to_double(pt[2]));
}

inline ES::V3d toVec3(const IK::Point_3 &pt)
{
  return ES::V3d(pt[0].to_double(), pt[1].to_double(), pt[2].to_double());
}

inline ES::V3d toVec3(const V3 &pt)
{
  return ES::V3d(pt[0].to_double(), pt[1].to_double(), pt[2].to_double());
}

inline K::Point_3 toP3_EK(const ES::V3d &pt)
{
  return K::Point_3(pt[0], pt[1], pt[2]);
}

inline K::Point_3 toP3_EK(const V3 &pt)
{
  return K::Point_3(IK_to_EK(pt[0]), IK_to_EK(pt[1]), IK_to_EK(pt[2]));
}

inline K::Point_3 toP3_EK(const IK::Point_3 &pt)
{
  return K::Point_3(IK_to_EK(pt[0]), IK_to_EK(pt[1]), IK_to_EK(pt[2]));
}

inline V3 toV3(const ES::V3d &pt)
{
  V3 ret;
  ret << pt[0], pt[1], pt[2];
  return ret;
}

inline V3 toV3(const K::Point_3 &pt)
{
  V3 ret;
  ret << EK_to_IK<real>(pt[0]), EK_to_IK<real>(pt[1]), EK_to_IK<real>(pt[2]);
  return ret;
}

inline IK::Point_3 toP3_IK(const K::Point_3 &pt)
{
  return IK::Point_3(EK_to_IK<real>(pt[0]), EK_to_IK<real>(pt[1]), EK_to_IK<real>(pt[2]));
}

#define BASIC_MATH_FUNC(_func)                                              \
  template<typename T>                                                      \
    requires SupportedRealType<T>                                           \
  T _func(const T &src)                                                     \
  {                                                                         \
    if constexpr (std::is_same<std::decay_t<T>, double>::value) {           \
      return std::_func(src);                                               \
    }                                                                       \
    else if constexpr (std::is_same<std::decay_t<T>, CGAL::Gmpfr>::value) { \
      T ret;                                                                \
      mpfr_##_func(ret.fr(), src.fr(), MPFR_RNDF);                          \
      return ret;                                                           \
    }                                                                       \
  }

BASIC_MATH_FUNC(sqrt)
BASIC_MATH_FUNC(acos)
BASIC_MATH_FUNC(abs)
BASIC_MATH_FUNC(sin)
BASIC_MATH_FUNC(cos)

inline void minFaceSize(const Poly &finalMesh)
{
  int min_c = 10;
  for (Poly::Halfedge_const_handle h = finalMesh.halfedges_begin(); h != finalMesh.halfedges_end(); ++h) {
    Poly::Halfedge_const_handle h1 = h;
    int count = 0;
    do {
      h1 = h1->next();
      count++;
    } while (h1 != h);

    min_c = std::min(min_c, count);
  }
  std::cout << "min c: " << min_c << std::endl;
}

inline void faceSize(Poly::Facet_const_handle f)
{
  Poly::Halfedge_const_handle h = f->halfedge();
  Poly::Halfedge_const_handle h1 = h;
  int count = 0;
  do {
    h1 = h1->next();
    count++;
  } while (h1 != h);

  std::cout << "f size: " << count << std::endl;
}

inline pgo::Mesh::TriMeshGeo dumpTri(const K::Triangle_2 &tri)
{
  ES::V3d p[3];
  for (int j = 0; j < 3; j++) {
    p[j][0] = tod(tri[j][0]);
    p[j][1] = tod(tri[j][1]);
    p[j][2] = 0;
  }

  return pgo::Mesh::createSingleTriangleMesh(p[0], p[1], p[2]);
}

inline pgo::Mesh::TriMeshGeo dumpTri(const K::Triangle_3 &tri)
{
  ES::V3d p[3];
  for (int j = 0; j < 3; j++) {
    p[j][0] = tod(tri[j][0]);
    p[j][1] = tod(tri[j][1]);
    p[j][2] = tod(tri[j][2]);
  }

  return pgo::Mesh::createSingleTriangleMesh(p[0], p[1], p[2]);
}

inline pgo::Mesh::TriMeshGeo dumpSeg(const K::Segment_2 &seg)
{
  ES::V3d p[2];
  for (int j = 0; j < 2; j++) {
    p[j][0] = tod(seg[j][0]);
    p[j][1] = tod(seg[j][1]);
    p[j][2] = 0;
  }

  return pgo::Mesh::createSingleTriangleMesh(p[0], p[1], p[0] + pgo::asVec3d(1e-4));
}

inline pgo::Mesh::TriMeshGeo dumpSeg(const K::Segment_3 &seg)
{
  ES::V3d p[2];
  for (int j = 0; j < 2; j++) {
    p[j][0] = tod(seg[j][0]);
    p[j][1] = tod(seg[j][1]);
    p[j][2] = tod(seg[j][2]);
  }

  return pgo::Mesh::createSingleTriangleMesh(p[0], p[1], p[0] + pgo::asVec3d(1e-4));
}

struct Patch2D
{
  // for bfs cutting
  std::list<K::Triangle_2> rawTriangles;
  // final vertex positions
  std::vector<K::Point_2> vertices;
  // final triangles, the final result will still be a single connected components
  std::vector<std::array<int, 3>> triangles;
  // the IDs of vertices on each edge
  std::vector<int> edgeVertices[3];
  // local triangle ID, j-th edge, target edge ID
  std::vector<std::tuple<int, int, int>> triangleEdgeIDs;
  // is vertex on target mesh vertex, or target mesh edge,
  // if on vertex, record target mesh vertex ID, 0
  // if on edge, record target mesh edge ID, 1
  std::vector<std::array<int, 2>> vertexIsOnTargetMesh;
  // corner vertex IDs
  int cornerVertices[3];
};

struct Patch3D
{
  std::list<K::Triangle_3> rawTriangles;
  std::vector<K::Point_3> vertices;
  std::vector<std::array<int, 3>> triangles;
  std::vector<int> edgeVertices[3];
  std::vector<std::tuple<int, int, int>> triangleEdgeIDs;
  // is vertex on target mesh vertex, or target mesh edge,
  // if on vertex, record target mesh vertex ID, 0
  // if on edge, record target mesh edge ID, 1
  std::vector<std::array<int, 2>> vertexIsOnTargetMesh;
  int cornerVertices[3];
};

enum class VertexAttachmentType : int
{
  DETACHED = -1,        // -1, detached
  ATTACHED = 0,         // 0: attached
  TARGET_VERTEX = 1,    // 1: target vertex
  TARGET_EDGE = 2,      // 2: target edge
  SKELETON_VERTEX = 4,  // 4: center
};

struct VertexProperty
{
  VertexAttachmentType vaType = VertexAttachmentType::DETACHED;
  int isBorder = 0;
  std::array<K::Point_3, 3> attachedTargets;
  std::array<int, 3> targetIDs;
};

using EdgeMap = std::unordered_set<std::pair<int, int>, pgo::EigenSupport::IntPairHash, pgo::EigenSupport::IntPairEqual>;

bool collapseEdge(Poly &finalMesh, Poly::Halfedge_handle h,
  const std::vector<VertexProperty> &vertexTypes, const EdgeMap &targetEdges);
bool collapseSmallTriangleDouble(Poly &finalMesh, const double &eps2,
  std::vector<VertexProperty> &vertexTypes, const EdgeMap &targetEdges);
bool collapseSmallEdge(Poly &finalMesh, const K::FT &eps2,
  const std::vector<VertexProperty> &vertexTypes, const EdgeMap &targetEdges);
bool collapseColinearVertices(Poly &finalMesh, int &iter, std::set<Poly::Halfedge_handle> &skippedHE,
  const std::vector<VertexProperty> &vertexTypes, const EdgeMap &targetEdges);
bool smoothVertices(Poly &finalMesh, std::vector<int> &isVertexChanged,
  const std::vector<VertexProperty> &vertexTypes, const EdgeMap &targetEdges);
void cleanUpMesh(const Poly &targetMeshPoly, const FaceGraphTree &targetMeshBVTreeExact,
  Poly &finalMesh, double edgeLengthThreshold,
  std::vector<VertexProperty> &vertexTypes, const EdgeMap &targetEdges);

void orientMesh(const std::vector<K::Point_3> &vertices, std::vector<std::vector<int>> &faces, const K::Point_3 centers[3], int type);
void orientMesh(const std::vector<K::Point_3> &vertices, std::vector<std::vector<int>> &faces, Poly &finalMesh);
bool triangulateMesh(Poly &finalMesh);

void saveExact(Poly &finalMesh, const char *filename);
}  // namespace MedialAxisRepresentation