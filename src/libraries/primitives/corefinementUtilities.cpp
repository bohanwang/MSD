#include "corefinementUtilities.h"

#include "createTriMesh.h"
#include "EigenSupport.h"
#include "pgoLogging.h"

#include <fmt/format.h>

#include <tbb/parallel_for.h>
#include <tbb/concurrent_hash_map.h>

namespace MedialAxisRepresentation
{
inline bool isFixed(VertexAttachmentType vaType)
{
  return vaType == VertexAttachmentType::TARGET_EDGE ||
    vaType == VertexAttachmentType::TARGET_VERTEX;
}

inline bool shareEdge(int vid0, int vid1,
  const std::vector<VertexProperty> &vertexTypes)
{
  if (vertexTypes[vid0].vaType == VertexAttachmentType::TARGET_EDGE && vertexTypes[vid1].vaType == VertexAttachmentType::TARGET_EDGE) {
    std::pair edgeID0{ vertexTypes[vid0].targetIDs[0], vertexTypes[vid0].targetIDs[1] };
    std::pair edgeID1{ vertexTypes[vid1].targetIDs[0], vertexTypes[vid1].targetIDs[1] };

    if (edgeID0.first > edgeID0.second)
      std::swap(edgeID0.first, edgeID0.second);

    if (edgeID1.first > edgeID1.second)
      std::swap(edgeID1.first, edgeID1.second);

    return edgeID0 == edgeID1;
  }
  else if (vertexTypes[vid0].vaType == VertexAttachmentType::TARGET_EDGE && vertexTypes[vid1].vaType == VertexAttachmentType::TARGET_VERTEX) {
    std::pair edgeID0{ vertexTypes[vid0].targetIDs[0], vertexTypes[vid0].targetIDs[1] };
    int tgtVID = vertexTypes[vid1].targetIDs[0];
    return (tgtVID == edgeID0.first || tgtVID == edgeID0.second);
  }
  else if (vertexTypes[vid0].vaType == VertexAttachmentType::TARGET_VERTEX && vertexTypes[vid1].vaType == VertexAttachmentType::TARGET_EDGE) {
    std::pair edgeID1{ vertexTypes[vid1].targetIDs[0], vertexTypes[vid1].targetIDs[1] };
    int tgtVID = vertexTypes[vid0].targetIDs[0];
    return (tgtVID == edgeID1.first || tgtVID == edgeID1.second);
  }
  else if (vertexTypes[vid0].isBorder == vertexTypes[vid1].isBorder) {
    return true;
  }

  return false;
}

inline bool coEdge(Poly::Halfedge_handle h0,
  Poly::Halfedge_handle g0, Poly::Halfedge_handle g1,
  const std::vector<VertexProperty> &vertexTypes)
{
  Poly::Vertex_handle v0 = h0->vertex();
  Poly::Vertex_handle v_l = g0->vertex();
  Poly::Vertex_handle v_r = g1->vertex();

  int vid0 = (int)v0->id();
  int vid_l = (int)v_l->id();
  int vid_r = (int)v_r->id();

  if (vertexTypes[vid0].vaType == VertexAttachmentType::TARGET_VERTEX)
    return false;

  return shareEdge(vid0, vid_l, vertexTypes) && shareEdge(vid0, vid_r, vertexTypes);
}

inline std::tuple<bool, K::Point_3, int> computeCollapsingPt(Poly::Halfedge_handle h,
  const std::vector<VertexProperty> &vertexTypes)
{
  Poly::Vertex_handle v0 = h->vertex();
  Poly::Vertex_handle v1 = h->opposite()->vertex();

  int vid0 = (int)v0->id();
  int vid1 = (int)v1->id();

  VertexAttachmentType vaType0 = vertexTypes[vid0].vaType;
  VertexAttachmentType vaType1 = vertexTypes[vid1].vaType;
  bool isBorder0 = vertexTypes[vid0].isBorder ? true : false;
  bool isBorder1 = vertexTypes[vid1].isBorder ? true : false;

  K::Point_3 pt(0, 0, 0);
  int vid = 0;
  bool canCollapse = false;

  if (vaType0 == VertexAttachmentType::SKELETON_VERTEX &&
    vaType1 == VertexAttachmentType::SKELETON_VERTEX) {
    throw std::runtime_error("Impossible");
  }

  // if either one of them is skeleton vertex, return false
  if (vaType0 == VertexAttachmentType::SKELETON_VERTEX ||
    vaType1 == VertexAttachmentType::SKELETON_VERTEX) {
    return std::tuple(false, pt, -1);
  }

  if (isBorder0 && isBorder1) {
    if (isFixed(vaType0) && isFixed(vaType1)) {
      canCollapse = false;
    }
    else if (isFixed(vaType0) && !isFixed(vaType1)) {
      pt = v0->point();
      vid = vid0;
      canCollapse = true;
    }
    else if (!isFixed(vaType0) && isFixed(vaType1)) {
      pt = v1->point();
      vid = vid1;
      canCollapse = true;
    }
    else {
      pt = CGAL::midpoint(v0->point(), v1->point());
      vid = vid0;
      canCollapse = true;
    }
  }
  else if (isBorder0 && !isBorder1) {
    if (isFixed(vaType0) && isFixed(vaType1)) {
      canCollapse = false;
    }
    else if (isFixed(vaType0) && !isFixed(vaType1)) {
      pt = v0->point();
      vid = vid0;
      canCollapse = true;
    }
    else if (!isFixed(vaType0) && isFixed(vaType1)) {
      canCollapse = false;
    }
    else {
      pt = v0->point();
      vid = vid0;
      canCollapse = true;
    }
  }
  else if (!isBorder0 && isBorder1) {
    if (isFixed(vaType0) && isFixed(vaType1)) {
      canCollapse = false;
    }
    else if (isFixed(vaType0) && !isFixed(vaType1)) {
      canCollapse = false;
    }
    else if (!isFixed(vaType0) && isFixed(vaType1)) {
      pt = v1->point();
      vid = vid1;
      canCollapse = true;
    }
    else {
      pt = v1->point();
      vid = vid1;
      canCollapse = true;
    }
  }
  else {
    if (isFixed(vaType0) && isFixed(vaType1)) {
      canCollapse = false;
    }
    else if (isFixed(vaType0) && !isFixed(vaType1)) {
      pt = v0->point();
      vid = vid0;
      canCollapse = true;
    }
    else if (!isFixed(vaType0) && isFixed(vaType1)) {
      pt = v1->point();
      vid = vid1;
      canCollapse = true;
    }
    else {
      pt = CGAL::midpoint(v0->point(), v1->point());
      vid = vid0;
      canCollapse = true;
    }
  }

  if (vaType0 == VertexAttachmentType::TARGET_EDGE && vaType1 == VertexAttachmentType::TARGET_EDGE &&
    isBorder0 == isBorder1) {
    std::pair edgeID0{ vertexTypes[vid0].targetIDs[0], vertexTypes[vid0].targetIDs[1] };
    std::pair edgeID1{ vertexTypes[vid1].targetIDs[0], vertexTypes[vid1].targetIDs[1] };

    if (edgeID0.first > edgeID0.second)
      std::swap(edgeID0.first, edgeID0.second);

    if (edgeID1.first > edgeID1.second)
      std::swap(edgeID1.first, edgeID1.second);

    // if both point on the same edge
    if (edgeID0 == edgeID1) {
      canCollapse = true;
      pt = CGAL::midpoint(v0->point(), v1->point());
      vid = vid0;
    }
  }
  else if (vaType0 == VertexAttachmentType::TARGET_EDGE && vaType1 == VertexAttachmentType::TARGET_VERTEX &&
    isBorder0 == isBorder1) {
    std::pair edgeID0{ vertexTypes[vid0].targetIDs[0], vertexTypes[vid0].targetIDs[1] };
    int tgtVID = vertexTypes[vid1].targetIDs[0];
    if (tgtVID == edgeID0.first || tgtVID == edgeID0.second) {
      canCollapse = true;
      pt = v1->point();
      vid = vid1;
    }
  }
  else if (vaType0 == VertexAttachmentType::TARGET_VERTEX && vaType1 == VertexAttachmentType::TARGET_EDGE &&
    isBorder0 == isBorder1) {
    std::pair edgeID1{ vertexTypes[vid1].targetIDs[0], vertexTypes[vid1].targetIDs[1] };
    int tgtVID = vertexTypes[vid0].targetIDs[0];
    if (tgtVID == edgeID1.first || tgtVID == edgeID1.second) {
      canCollapse = true;
      pt = v0->point();
      vid = vid0;
    }
  }

  return std::tuple(canCollapse, pt, vid);
}

bool checkEdgeCollapsingFlip(Poly &finalMesh, Poly::Halfedge_handle h, const K::Point_3 oldPos[2], const K::Point_3 &newPos)
{
  Poly::Halfedge_handle initial_h[2] = { h, h->opposite() };

  for (int j = 0; j < 2; j++) {
    Poly::Halfedge_handle h0 = initial_h[j];
    bool flipped = false;
    do {
      const K::Point_3 &p1 = h0->next()->vertex()->point();
      const K::Point_3 &p2 = h0->next()->next()->vertex()->point();

      if (oldPos[j] == p1 || p1 == p2 || oldPos[j] == p2) {
      }
      else {
        K::Vector_3 e01 = p1 - oldPos[j];
        K::Vector_3 e02 = p2 - oldPos[j];
        K::Vector_3 n1 = CGAL::cross_product(e01, e02);

        K::Vector_3 e01_new = p1 - newPos;
        K::Vector_3 e02_new = p2 - newPos;
        K::Vector_3 n2 = CGAL::cross_product(e01_new, e02_new);

        if (CGAL::scalar_product(n1, n2) < 0) {
          flipped = true;
          break;
        }
      }

      h0 = h0->next()->opposite();
    } while (h0 != initial_h[j]);

    if (flipped) {
      return false;
    }
  }

  return true;
}
}  // namespace MedialAxisRepresentation

bool MedialAxisRepresentation::collapseEdge(Poly &finalMesh, Poly::Halfedge_handle h,
  const std::vector<VertexProperty> &vertexTypes, const EdgeMap &targetEdges)
{
  constexpr int debug = 0;

  auto collapseRet = computeCollapsingPt(h, vertexTypes);

  // if (h->vertex()->id() && h->opposite()->vertex()->id())
  //   return false;

  if (std::get<0>(collapseRet) == false)
    return false;

  // check if the edge can be collapsed
  if (CGAL::Euler::does_satisfy_link_condition(CGAL::edge(h, finalMesh), finalMesh)) {
    Poly::Vertex_handle v = CGAL::Euler::collapse_edge(CGAL::edge(h, finalMesh), finalMesh);
    v->point() = std::get<1>(collapseRet);

    return true;
  }
  // if no
  else {
#if 0
    // Poly finalMeshCopy = finalMesh;
    bool h_face_exist = true;
    bool h_o_face_exist = true;
    Poly::Halfedge_handle g = h->opposite();

    // if two triangles share one edge
    auto shareEdge = [&](Poly::Halfedge_handle h1, Poly::Halfedge_handle h2) -> bool {
      if (h1->is_border() || h2->is_border()) {
        return false;
      }

      if (h_face_exist &&
        (h1->facet() == h->facet() || h2->facet() == h->facet())) {
        return false;
      }

      if (h_o_face_exist &&
        (h1->facet() == g->facet() ||
          h2->facet() == g->facet())) {
        return false;
      }

      Poly::Halfedge_handle hh = h1;
      for (int i = 0; i < 3; i++, hh = hh->next()) {
        Poly::Halfedge_handle gg = h2;
        for (int j = 0; j < 3; j++, gg = gg->next()) {
          if (hh->opposite() == gg) {
            if ((hh->vertex() == h1->vertex() && gg->vertex() == h2->vertex()) ||
              (hh->vertex() == h2->vertex() && gg->vertex() == h1->vertex())) {
              return false;
            }
            else {
              return true;
            }
          }
        }
      }

      return false;
    };

    auto sameFace = [](Poly::Halfedge_handle h1, Poly::Halfedge_handle h2) -> bool {
      if (h1->is_border() || h2->is_border()) {
        return false;
      }

      Poly::Halfedge_handle g = h1;
      do {
        if (g == h2) {
          return true;
        }
        g = g->next();
      } while (g != h1);

      return false;
    };

    Poly::Halfedge_handle v0v1 = h;
    Poly::Halfedge_handle v1v0 = h->opposite();
    Poly::Halfedge_handle v0v1_1 = v0v1;
    std::ofstream("aa.off") << finalMesh;
    while (1) {
      bool modified = false;
      Poly::Halfedge_handle v1v0_1 = v1v0;

      TriMeshGeo mam;
      mam.addMesh(createSingleTriangleMesh(toVec3(v0v1_1->vertex()->point()),
        toVec3(v0v1_1->opposite()->vertex()->point()) + pgo::asVec3d(5e-4),
        toVec3(v0v1_1->opposite()->vertex()->point()) - pgo::asVec3d(5e-4)));
      mam.save("aa1.obj");

      mam.clear();
      mam.addMesh(createSingleTriangleMesh(toVec3(v1v0_1->vertex()->point()),
        toVec3(v1v0_1->opposite()->vertex()->point()) + pgo::asVec3d(5e-4),
        toVec3(v1v0_1->opposite()->vertex()->point()) - pgo::asVec3d(5e-4)));
      mam.save("aa2.obj");

      while (1) {
        // if two halfedges from both point are in the same triangle
        if (sameFace(v0v1_1, v1v0_1)) {
          TriMeshGeo mam;
          Poly::Halfedge_handle zzz = v0v1_1;
          do {
            mam.addMesh(createSingleTriangleMesh(toVec3(zzz->vertex()->point()),
              toVec3(zzz->opposite()->vertex()->point()) + pgo::asVec3d(1e-4),
              toVec3(zzz->opposite()->vertex()->point()) - pgo::asVec3d(1e-4)));

            zzz = zzz->next();
          } while (zzz != v0v1_1);
          mam.save("aaz.obj");

          Poly::Halfedge_handle next01 = v0v1_1->next()->opposite();
          Poly::Halfedge_handle next10 = v1v0_1->next()->opposite();

          if (next01->is_border())
            next01 = next01->next()->opposite();

          if (next10->is_border())
            next10 = next10->next()->opposite();

          std::cout << next01->is_border() << ',' << next10->is_border() << std::endl;
          std::cout << (next01->facet() == v0v1_1->facet()) << std::endl;
          std::cout << (next10->facet() == v0v1_1->facet()) << std::endl;

          if (h_face_exist && sameFace(h, v0v1_1))
            h_face_exist = false;

          if (h_o_face_exist && sameFace(g, v0v1_1))
            h_o_face_exist = false;

          CGAL::Euler::remove_face(v0v1_1, finalMesh);

          v0v1_1 = next01;
          v0v1 = next01;

          v1v0 = next10;
          v1v0_1 = next10;

          modified = true;
          TriMeshGeo m;
          pgo::CGALInterface::polyhedron2TriangleMesh(finalMesh, m);
          m.save("aa.obj");
        }
        // if two halfedges pointing triangles share the same edge
        else if (shareEdge(v0v1_1, v1v0_1)) {
          // remove first face
          Poly::Halfedge_handle next01 = v0v1_1->next()->opposite();
          Poly::Halfedge_handle next10 = v1v0_1->next()->opposite();

          if (sameFace(next01, v1v0_1)) {
            next01 = next01->next()->opposite();
          }

          if (next01->is_border())
            next01 = next01->next()->opposite();

          if (sameFace(next10, v0v1_1)) {
            next10 = next10->next()->opposite();
          }

          if (next10->is_border())
            next10 = next10->next()->opposite();

          CGAL::Euler::remove_face(v0v1_1, finalMesh);
          CGAL::Euler::remove_face(v1v0_1, finalMesh);

          v0v1_1 = next01;
          v0v1 = next01;

          v1v0 = next10;
          v1v0_1 = next10;

          modified = true;
          TriMeshGeo m;
          pgo::CGALInterface::polyhedron2TriangleMesh(finalMesh, m);
          m.save("aa.obj");
        }
        else {
          v1v0_1 = v1v0_1->next()->opposite();

          if (v1v0_1 == v1v0) {
            break;
          }
        }
      }

      if (!modified) {
        v0v1_1 = v0v1_1->next()->opposite();
        if (v0v1_1 == v0v1) {
          break;
        }
      }
    }

    // todo
    // reconstruct holes
    while (v0v1->is_border() == false)
      v0v1 = v0v1->next()->opposite();

    while (v1v0->is_border() == false)
      v1v0 = v1v0->next()->opposite();

    Poly::Halfedge_handle hh = finalMesh.add_facet_to_border(v0v1, v1v0);
    hh = finalMesh.fill_hole(hh->opposite());

    bool linkCond = CGAL::Euler::does_satisfy_link_condition(CGAL::edge(hh, finalMesh), finalMesh);
    std::cout << "link cond: " << linkCond << std::endl;
    std::cout << "is closed: " << finalMesh.is_closed() << std::endl;
    if (linkCond) {
      CGAL::Euler::collapse_edge(CGAL::edge(hh, finalMesh), finalMesh);
      std::ofstream("aa.off") << finalMesh;

      return true;
    }
    else {
      std::cin.get();
      return false;
      // finalMesh = finalMeshCopy;
    }
#else
    Poly finalMeshCopy = finalMesh;
    std::set<Poly::Halfedge_handle> borderHE;
    for (auto hit = finalMesh.halfedges_begin(); hit != finalMesh.halfedges_end(); ++hit) {
      if (hit->is_border()) {
        borderHE.emplace(hit);
      }
    }

    Poly::Halfedge_handle v0v1 = h;
    Poly::Halfedge_handle v1v0 = h->opposite();
    Poly::Halfedge_handle v0v1_1 = v0v1;

    Poly::Vertex_handle v1 = v0v1->vertex();
    Poly::Vertex_handle v0 = v1v0->vertex();

    if constexpr (debug) {
      pgo::Mesh::TriMeshGeo mm;
      pgo::CGALInterface::polyhedron2TriangleMesh(finalMesh, mm);
      mm.save("aa.obj");

      pgo::Mesh::TriMeshGeo aaaa;
      aaaa.addMesh(pgo::Mesh::createSingleTriangleMesh(toVec3(v0->point()), toVec3(v1->point()), toVec3(v0->point()) + pgo::asVec3d(1e-6)));
      aaaa.save("cc.obj");
    }

    std::set<Poly::Facet_handle> facesToDelete;
    std::set<Poly::Halfedge_handle> halfEdgeVisited;

    std::unordered_set<Poly::Facet_handle> faceSet;
    std::vector<std::tuple<Poly::Halfedge_handle, int>> Q;

    faceSet.reserve(finalMesh.size_of_facets());
    Q.reserve(finalMesh.size_of_facets() * 2);

    while (1) {
      if (halfEdgeVisited.find(v0v1_1) != halfEdgeVisited.end()) {
        v0v1_1 = v0v1_1->next()->opposite();
        if (v0v1_1 == v0v1)
          break;
        else
          continue;
      }

      Poly::Halfedge_handle v1v0_1 = v1v0;

      if constexpr (debug) {
        pgo::Mesh::TriMeshGeo mam;
        mam.addMesh(pgo::Mesh::createSingleTriangleMesh(toVec3(v0v1_1->vertex()->point()),
          toVec3(v0v1_1->opposite()->vertex()->point()) + pgo::asVec3d(1e-4),
          toVec3(v0v1_1->opposite()->vertex()->point()) - pgo::asVec3d(1e-4)));
        mam.save("aa1.obj");
      }

      while (1) {
        if (halfEdgeVisited.find(v1v0_1) != halfEdgeVisited.end()) {
          v1v0_1 = v1v0_1->next()->opposite();
          if (v1v0_1 == v1v0)
            break;
          else
            continue;
        }

        if constexpr (debug) {
          pgo::Mesh::TriMeshGeo mam;
          mam.addMesh(pgo::Mesh::createSingleTriangleMesh(toVec3(v1v0_1->vertex()->point()),
            toVec3(v1v0_1->opposite()->vertex()->point()) + pgo::asVec3d(1e-4),
            toVec3(v1v0_1->opposite()->vertex()->point()) - pgo::asVec3d(1e-4)));
          mam.save("aa2.obj");
        }

        // if there is a case where two edges form a triangle
        if (v0v1_1->opposite()->vertex() == v1v0_1->opposite()->vertex()) {
          faceSet.clear();
          Q.clear();

          Q.emplace_back(h, 0);
          Q.emplace_back(v0v1_1->opposite(), 0);
          Q.emplace_back(v1v0_1, 0);

          Q.emplace_back(h->opposite(), 1);
          Q.emplace_back(v0v1_1, 1);
          Q.emplace_back(v1v0_1->opposite(), 1);

          for (auto h : Q) {
            faceSet.emplace(std::get<0>(h)->facet());
          }

          size_t start = 0;
          while (start < Q.size()) {
            Poly::Halfedge_handle h_cur = std::get<0>(Q[start]);
            int setID_cur = std::get<1>(Q[start]);
            start++;

            Poly::Halfedge_handle h1 = h_cur->next();
            Poly::Halfedge_handle h2 = h1->next();

            if (faceSet.find(h1->opposite()->facet()) == faceSet.end()) {
              Q.emplace_back(h1->opposite(), setID_cur);
              faceSet.emplace(h1->opposite()->facet());
            }

            if (faceSet.find(h2->opposite()->facet()) == faceSet.end()) {
              Q.emplace_back(h2->opposite(), setID_cur);
              faceSet.emplace(h2->opposite()->facet());
            }
          }

          int counter[2] = { 0, 0 };
          for (int i = 0; i < (int)Q.size(); i++) {
            counter[std::get<1>(Q[i])]++;
          }

          if (counter[0] < counter[1]) {
            for (int i = 0; i < (int)Q.size(); i++) {
              if (std::get<1>(Q[i]) == 0) {
                facesToDelete.emplace(std::get<0>(Q[i])->facet());
              }
            }
          }
          else {
            for (int i = 0; i < (int)Q.size(); i++) {
              if (std::get<1>(Q[i]) == 1) {
                facesToDelete.emplace(std::get<0>(Q[i])->facet());
              }
            }
          }

          if constexpr (debug) {
            pgo::Mesh::TriMeshGeo m1;
            std::map<Poly::Vertex_handle, int> vids;
            for (auto fh : facesToDelete) {
              Poly::Halfedge_handle h0 = fh->halfedge();
              Poly::Halfedge_handle h1 = h0->next();
              Poly::Halfedge_handle h2 = h1->next();

              vids.emplace(h0->vertex(), 0);
              vids.emplace(h1->vertex(), 0);
              vids.emplace(h2->vertex(), 0);
            }

            for (int inc = 0; auto &pr : vids) {
              ES::V3d v0 = toVec3(pr.first->point());
              m1.addPos(v0);
              pr.second = inc++;
            }

            for (auto fh : facesToDelete) {
              Poly::Halfedge_handle h0 = fh->halfedge();
              Poly::Halfedge_handle h1 = h0->next();
              Poly::Halfedge_handle h2 = h1->next();

              m1.addTri(ES::V3i(vids[h0->vertex()], vids[h1->vertex()], vids[h2->vertex()]));
            }
            m1.save("bb.obj");
          }

          halfEdgeVisited.clear();
          for (auto fh : facesToDelete) {
            halfEdgeVisited.emplace(fh->halfedge());
            halfEdgeVisited.emplace(fh->halfedge()->next());
            halfEdgeVisited.emplace(fh->halfedge()->next()->next());
          }
        }

        v1v0_1 = v1v0_1->next()->opposite();
        if (v1v0_1 == v1v0) {
          break;
        }
      }

      v0v1_1 = v0v1_1->next()->opposite();
      if (v0v1_1 == v0v1) {
        break;
      }
    }

    bool canCollapse = true;
    for (auto fh : facesToDelete) {
      Poly::Halfedge_handle h0 = fh->halfedge();
      Poly::Halfedge_handle h1 = h0;
      do {
        if (h1->vertex()->id()) {
          canCollapse = false;
          break;
        }
        h1 = h1->next();
      } while (h1 != h0);
    }

    if constexpr (debug) {
      pgo::Mesh::TriMeshGeo m1;
      std::map<Poly::Vertex_handle, int> vids;
      for (auto fh : facesToDelete) {
        Poly::Halfedge_handle h0 = fh->halfedge();
        Poly::Halfedge_handle h1 = h0->next();
        Poly::Halfedge_handle h2 = h1->next();

        vids.emplace(h0->vertex(), 0);
        vids.emplace(h1->vertex(), 0);
        vids.emplace(h2->vertex(), 0);
      }

      for (int inc = 0; auto &pr : vids) {
        ES::V3d v0 = toVec3(pr.first->point());
        m1.addPos(v0);
        pr.second = inc++;
      }

      for (auto fh : facesToDelete) {
        Poly::Halfedge_handle h0 = fh->halfedge();
        Poly::Halfedge_handle h1 = h0->next();
        Poly::Halfedge_handle h2 = h1->next();
        m1.addTri(ES::V3i(vids[h0->vertex()], vids[h1->vertex()], vids[h2->vertex()]));
      }
      m1.save("bb.obj");
    }

    for (auto fh : facesToDelete) {
      finalMesh.erase_facet(fh->halfedge());
    }

    if constexpr (0) {
      pgo::Mesh::TriMeshGeo mm;
      pgo::CGALInterface::polyhedron2TriangleMesh(finalMesh, mm);
      mm.save("aa.obj");
    }

    Poly::Halfedge_handle border_h;
    for (border_h = finalMesh.halfedges_begin(); border_h != finalMesh.halfedges_end(); ++border_h) {
      if (borderHE.find(border_h) != borderHE.end())
        continue;

      if (border_h->is_border())
        break;
    }

    Poly::Halfedge_handle border_h1 = border_h;
    do {
      if (border_h1->vertex() == v1)
        break;

      border_h1 = border_h1->next();
    } while (border_h1 != border_h);

    PGO_ALOG(border_h1->vertex() == v1);

    border_h = border_h1;
    do {
      if (border_h->vertex() == v0)
        break;

      border_h = border_h->next();
    } while (border_h1 != border_h);

    PGO_ALOG(border_h->vertex() == v0);

    Poly::Halfedge_handle hh = finalMesh.add_facet_to_border(border_h, border_h1);
    hh = finalMesh.fill_hole(hh->opposite());

    // std::vector<Poly::Facet_handle> newFaces;
    // CGAL::Polygon_mesh_processing::triangulate_hole(finalMesh, hh->opposite(), std::back_inserter(newFaces));
    // std::cout << "new faces: " << newFaces.size() << std::endl;
    std::cout << "is triangle: " << finalMesh.is_pure_triangle() << std::endl;

    if constexpr (debug) {
      pgo::Mesh::TriMeshGeo mm;
      pgo::CGALInterface::polyhedron2TriangleMesh(finalMesh, mm);
      mm.save("aa.obj");
    }

    bool linkCond = CGAL::Euler::does_satisfy_link_condition(CGAL::edge(hh, finalMesh), finalMesh);
    std::cout << "link cond: " << linkCond << std::endl;
    std::cout << "is closed: " << finalMesh.is_closed() << std::endl;
    std::cout << "is triangle: " << finalMesh.is_pure_triangle() << std::endl;

    if (linkCond) {
      Poly::Vertex_handle v = CGAL::Euler::collapse_edge(CGAL::edge(hh, finalMesh), finalMesh);
      v->point() = std::get<1>(collapseRet);

      std::cout << "after collapsing\n";
      std::cout << "is closed: " << finalMesh.is_closed() << std::endl;
      std::cout << "is triangle: " << finalMesh.is_pure_triangle() << std::endl;

      if constexpr (debug) {
        pgo::Mesh::TriMeshGeo mm;
        pgo::CGALInterface::polyhedron2TriangleMesh(finalMesh, mm);
        mm.save("aa.obj");
      }

      return true;
    }
    else {
      std::cin.get();
      finalMesh = finalMeshCopy;
      return false;
    }
#endif
  }
}

bool MedialAxisRepresentation::collapseSmallTriangleDouble(Poly &finalMesh, const double &eps2,
  std::vector<VertexProperty> &vertexTypes, const EdgeMap &targetEdges)
{
  using IK = pgo::CGALInterface::KernelSimple;
  using EK_to_IK = CGAL::Cartesian_converter<K, IK>;
  EK_to_IK to_inexact;

  using IK_to_EK = CGAL::Cartesian_converter<IK, K>;
  IK_to_EK to_exact;

  constexpr int debug = 0;
  static int iii = 0;

  for (auto it : finalMesh.facet_handles()) {
    if (it->id() & 1) {
      continue;
    }

    Poly::Halfedge_handle h[3] = { it->halfedge() };
    h[1] = h[0]->next();
    h[2] = h[1]->next();

    bool isModified = false;

    for (int j = 0; j < 3; j++) {
      IK::Point_3 vtx_pos[3] = {
        to_inexact(h[j % 3]->vertex()->point()),
        to_inexact(h[(j + 1) % 3]->vertex()->point()),
        to_inexact(h[(j + 2) % 3]->vertex()->point())
      };

      // check intersection point location, if it is on the segment or not.
      // we only want the ones with it on
      // (t (p1 - p0) + p0 - x)T (p1 - p0)
      // t (p1- p0)^T(p1 - p0
      IK::Vector_3 dir = vtx_pos[2] - vtx_pos[1];
      IK::Vector_3 diff = vtx_pos[1] - vtx_pos[0];

      IK::FT a = dir.squared_length();
      if (CGAL::to_double(CGAL::abs(a)) < 1e-16) {
        std::cout << "a too smal a=" << CGAL::to_double(a) << std::endl;
        continue;
      }

      IK::Line_3 line(vtx_pos[1], vtx_pos[2]);
      IK::FT d2 = CGAL::squared_distance(vtx_pos[0], line);
      if (d2 > eps2) {
        continue;
      }

      std::cout << "Find one bad triangle: " << CGAL::to_double(d2) << std::endl;

      IK::FT b = CGAL::scalar_product(dir, diff);
      IK::FT t = -b / a;

      if (t <= 0 || t >= 1) {
        std::cout << "t not in [0, 1]; t=" << CGAL::to_double(t) << std::endl;
        continue;
      }
      else {
        std::cout << "t in [0, 1]; t=" << CGAL::to_double(t) << std::endl;
      }

      if constexpr (0) {
        static int jjj = 0;

        Poly::Halfedge_handle h_oppo = h[(j + 2) % 3]->opposite()->next();

        pgo::Mesh::TriMeshGeo aar;
        ES::V3d p0 = toVec3(h[j]->vertex()->point());
        ES::V3d p1 = toVec3(h[(j + 1) % 3]->vertex()->point());
        ES::V3d p2 = toVec3(h[(j + 2) % 3]->vertex()->point());
        ES::V3d p3 = toVec3(h_oppo->vertex()->point());

        aar.addPos(p0);
        aar.addPos(p1);
        aar.addPos(p2);
        aar.addPos(p3);

        aar.addTri(ES::V3i(0, 1, 2));
        aar.addTri(ES::V3i(2, 1, 3));

        aar.save(fmt::format("aa{}.obj", jjj));

        jjj = (jjj + 1) % 10;
      }

      // new vtx pos
      IK::Point_3 newPtIK = vtx_pos[1] + dir * t;
      bool hasFixedVertices = false;
      if (isFixed(vertexTypes[h[j]->vertex()->id()].vaType)) {
        if (vertexTypes[h[(j + 1) % 3]->vertex()->id()].isBorder &&
          vertexTypes[h[(j + 2) % 3]->vertex()->id()].isBorder) {
          std::cout << "Middle edge is a border edge. Skip." << std::endl;
          continue;
        }

        if (vertexTypes[h[(j + 1) % 3]->vertex()->id()].vaType == VertexAttachmentType::TARGET_VERTEX &&
          vertexTypes[h[(j + 2) % 3]->vertex()->id()].vaType == VertexAttachmentType::TARGET_VERTEX) {
          // get target mesh vid
          int vid0 = vertexTypes[h[(j + 1) % 3]->vertex()->id()].targetIDs[0];
          int vid1 = vertexTypes[h[(j + 2) % 3]->vertex()->id()].targetIDs[0];

          std::pair<int, int> edgeID{ vid0, vid1 };
          if (edgeID.first > edgeID.second)
            std::swap(edgeID.first, edgeID.second);

          if (targetEdges.find(edgeID) != targetEdges.end()) {
            std::cout << "middle edge is a target edge. skip." << std::endl;
            continue;
          }
        }

        hasFixedVertices = true;
        newPtIK = to_inexact(h[j]->vertex()->point());
      }
      else if (vertexTypes[h[j]->vertex()->id()].isBorder) {
        if (isFixed(vertexTypes[h[(j + 1) % 3]->vertex()->id()].vaType) &&
          isFixed(vertexTypes[h[(j + 2) % 3]->vertex()->id()].vaType)) {
          if (shareEdge((int)h[(j + 1) % 3]->vertex()->id(), (int)h[(j + 2) % 3]->vertex()->id(), vertexTypes)) {
            std::cout << "border + fixed, nothing can move" << std::endl;
            continue;
          }
        }
        /*
        * &&
        isFixed(vertexTypes[h[(j + 1) % 3]->vertex()->id()].vaType) &&
        isFixed(vertexTypes[h[(j + 2) % 3]->vertex()->id()].vaType)
        */
        // std::cout << "border" << std::endl;
        // continue;
        hasFixedVertices = true;
        newPtIK = to_inexact(h[j]->vertex()->point());
      }
      else {
        hasFixedVertices = false;
      }

      bool isMovingValid = true;

      // check if the current vertex is moved to the target location,
      // all related triangles has small area
      Poly ::Halfedge_handle hh = h[j]->next()->opposite();
      do {
        Poly::Halfedge_handle hv0 = hh;
        Poly::Halfedge_handle hv1 = hv0->next();
        Poly::Halfedge_handle hv2 = hv1->next();

        IK::Point_3 this_tri_vtx[3] = {
          to_inexact(hv0->vertex()->point()),
          to_inexact(hv1->vertex()->point()),
          to_inexact(hv2->vertex()->point()),
        };

        if (hv0->is_border()) {
        }
        else {
          IK::Point_3 v_new[3] = { newPtIK, this_tri_vtx[1], this_tri_vtx[2] };
          IK::Point_3 v_old[3] = { this_tri_vtx[0], this_tri_vtx[1], this_tri_vtx[2] };
          IK::Triangle_3 tri_new(v_new[0], v_new[1], v_new[2]);
          IK::FT area2_new = tri_new.squared_area() * 4;

          IK::Triangle_3 tri_old(v_old[0], v_old[1], v_old[2]);
          IK::FT area2_old = tri_old.squared_area() * 4;

          // check old
          bool isOldSmall = false;
          for (int k = 0; k < 3; k++) {
            IK::FT d2 = CGAL::squared_distance(v_old[k], v_old[(k + 1) % 3]);
            if (d2 < 1e-16) {
              // std::cout << "old" << CGAL::to_double(d2) << ',' << CGAL::to_double(area2_old) << std::endl;
              isOldSmall = true;
              break;
            }
            IK::FT h2 = area2_old / d2;

            if (h2 < eps2) {
              isOldSmall = true;
              break;
            }
          }

          // check new
          bool isNewSmall = false;
          for (int k = 0; k < 3; k++) {
            IK::FT d2 = CGAL::squared_distance(v_new[k], v_new[(k + 1) % 3]);
            if (d2 < 1e-16) {
              isNewSmall = true;
              break;
            }
            IK::FT h2 = area2_new / d2;

            if (h2 < eps2) {
              isNewSmall = true;
              break;
            }
          }
          // if the old triangle is not small but the new one is small
          if (isOldSmall == false && isNewSmall) {
            isMovingValid = false;
            break;
          }
        }

        hh = hh->next()->opposite();
      } while (hh != h[j]);

      if (isMovingValid == false) {
        continue;
      }

      // if the opposite is a border
      if (h[(j + 2) % 3]->opposite()->is_border()) {
        throw std::runtime_error("Impossible");
#if 0
        // create a point in the middle of the segment
        Poly::Halfedge_handle hnew = finalMesh.split_edge(h[(j + 2) % 3]);
        if (hasFixedVertices) {
          hnew->vertex()->point() = h[j]->vertex()->point();
        }
        else {
          K::Vector_3 dirEK = h[(j + 2) % 3]->vertex()->point() - h[(j + 1) % 3]->vertex()->point();
          K::Point_3 newPtEK = h[(j + 1) % 3]->vertex()->point() + dirEK * t;
          hnew->vertex()->point() = newPtEK;
        }

        // create two faces
        Poly::Halfedge_handle h1 = finalMesh.split_facet(hnew, h[j]);

        // mark dirty
        Poly ::Halfedge_handle hh = h[j];
        do {
          if (hh->is_border() == false)
            hh->facet()->id() = 0;

          hh = hh->next()->opposite();
        } while (hh != h[j]);

        if (collapseEdge(finalMesh, h1, vertexTypes, targetEdges)) {
          // std::cout << "After splitting: is triangle mesh: " << finalMesh.is_pure_triangle() << std::endl;
          std::cout << "Triangle collapsed. " << std::endl;
          isModified = true;
          break;
        }
        else {
          std::cout << "Failed to collapse." << std::endl;
        }
#endif
      }
      // otherwise
      else {
        // check opposite triangle, if the area is small too, then we skip this temporarily
        Poly::Halfedge_handle h_oppo = h[(j + 2) % 3]->opposite()->next();

        IK::Point_3 vtx_oppo = to_inexact(h_oppo->vertex()->point());
        IK::FT oppo_area[2] = {
          CGAL::squared_area(newPtIK, vtx_pos[1], vtx_oppo) * 4.0,
          CGAL::squared_area(newPtIK, vtx_pos[2], vtx_oppo) * 4.0,
        };

        // check if all new opposite triangles are big enough
        IK::Vector_3 segs[6] = {
          newPtIK - vtx_pos[1],
          vtx_pos[2] - newPtIK,
          vtx_oppo - newPtIK,
          vtx_oppo - newPtIK,
          vtx_oppo - vtx_pos[1],
          vtx_oppo - vtx_pos[2]
        };
        int segTriID[6] = { 0, 1, 0, 1, 0, 1 };

        bool isOppoValid = true;
        for (int i = 0; i < 6; i++) {
          IK::FT d2 = segs[i].squared_length();
          const IK::FT &a2 = oppo_area[segTriID[i]];
          if (d2 < 1e-16) {
            std::cout << "Opposite triangle seg is too small: d2=" << CGAL::to_double(d2) << std::endl;
            isOppoValid = false;
            break;
          }

          IK::FT h2 = a2 / d2;
          if (h2 < eps2) {
            std::cout << "Opposite triangle is too small: h2=" << CGAL::to_double(h2) << std::endl;
            isOppoValid = false;
            break;
          }
          else {
            std::cout << "Opposite triangle h2=" << CGAL::to_double(h2) << std::endl;
          }
        }

        if (isOppoValid == false) {
          continue;
        }

        // create a point in the middle of the segment
        Poly::Halfedge_handle hnew = finalMesh.split_edge(h[(j + 2) % 3]);
        Poly::Halfedge_handle hnew_oppo = hnew->next()->opposite();

        K::Vector_3 dirEK = h[(j + 2) % 3]->vertex()->point() - h[(j + 1) % 3]->vertex()->point();
        K::Point_3 newPtEK = h[(j + 1) % 3]->vertex()->point() + dirEK * t;

        K::Point_3 oldPositions[2] = {
          newPtEK,
          h[j]->vertex()->point(),
        };

        if (hasFixedVertices) {
          hnew->vertex()->point() = h[j]->vertex()->point();
          hnew->vertex()->id() = h[j]->vertex()->id();
        }
        else {
          hnew->vertex()->point() = newPtEK;

          hnew->vertex()->id() = vertexTypes.size();
          vertexTypes.emplace_back();

          if (vertexTypes[h[(j + 1) % 3]->vertex()->id()].isBorder &&
            vertexTypes[h[(j + 2) % 3]->vertex()->id()].isBorder) {
            vertexTypes.back().isBorder = 1;
          }

          if (vertexTypes[h[(j + 1) % 3]->vertex()->id()].vaType == VertexAttachmentType::TARGET_VERTEX &&
            vertexTypes[h[(j + 2) % 3]->vertex()->id()].vaType == VertexAttachmentType::TARGET_VERTEX) {
            // get target mesh vid
            int vid0 = vertexTypes[h[(j + 1) % 3]->vertex()->id()].targetIDs[0];
            int vid1 = vertexTypes[h[(j + 2) % 3]->vertex()->id()].targetIDs[0];

            std::pair<int, int> edgeID{ vid0, vid1 };
            if (edgeID.first > edgeID.second)
              std::swap(edgeID.first, edgeID.second);

            if (targetEdges.find(edgeID) != targetEdges.end()) {
              vertexTypes.back().vaType = VertexAttachmentType::TARGET_EDGE;
              vertexTypes.back().targetIDs[0] = edgeID.first;
              vertexTypes.back().targetIDs[1] = edgeID.second;
            }
          }
          else if (isFixed(vertexTypes[h[(j + 1) % 3]->vertex()->id()].vaType) &&
            isFixed(vertexTypes[h[(j + 2) % 3]->vertex()->id()].vaType)) {
            if (shareEdge((int)h[(j + 1) % 3]->vertex()->id(), (int)h[(j + 2) % 3]->vertex()->id(), vertexTypes)) {
              vertexTypes.back().vaType = VertexAttachmentType::TARGET_EDGE;

              int vid0 = vertexTypes[h[(j + 1) % 3]->vertex()->id()].targetIDs[0];
              int vid1 = vertexTypes[h[(j + 2) % 3]->vertex()->id()].targetIDs[0];

              std::pair<int, int> edgeID{ vid0, vid1 };
              if (edgeID.first > edgeID.second)
                std::swap(edgeID.first, edgeID.second);

              vertexTypes.back().targetIDs[0] = edgeID.first;
              vertexTypes.back().targetIDs[1] = edgeID.second;
            }
          }
        }

        // Poly finalMeshCopy = finalMesh;

        Poly::Halfedge_handle h1 = finalMesh.split_facet(hnew, h[j]);
        h[j]->vertex()->point() = hnew->vertex()->point();

        // mark dirty
        Poly ::Halfedge_handle hh = h[j];
        do {
          hh->facet()->id() = 0;
          hh = hh->next()->opposite();
        } while (hh != h[j]);

        Poly::Halfedge_handle g = hnew_oppo->next()->next();
        Poly::Halfedge_handle h2 = finalMesh.split_facet(hnew_oppo, g);

        if (checkEdgeCollapsingFlip(finalMesh, h1, oldPositions, hnew->vertex()->point()) == false) {
          std::cout << "is tri: " << finalMesh.is_pure_triangle() << std::endl;
          continue;
        }

        if constexpr (debug) {
          Poly pp;
          pp.make_triangle(h[0]->vertex()->point(), h[1]->vertex()->point(), h[2]->vertex()->point());
          pp.make_triangle(h[0]->vertex()->point(), h[1]->vertex()->point(), h[2]->vertex()->point());
          std::ofstream("aa1.off") << pp;
        }

        // mark dirty
        h2->facet()->id() = 0;
        h2->opposite()->facet()->id() = 0;

        if (collapseEdge(finalMesh, h1, vertexTypes, targetEdges)) {
          // std::cout << "After splitting: is triangle mesh: " << finalMesh.is_pure_triangle() << std::endl;
          std::cout << "Triangle collapsed. " << std::endl;

          // TriMeshGeo mm;
          // pgo::CGALInterface::polyhedron2TriangleMesh(finalMesh, mm);
          // iii++;
          // if (iii > 10) {
          //   mm.save("rr.obj");
          //   iii = 0;
          // }

          isModified = true;
          break;
        }
        else {
          std::cout << "Failed to collapse." << std::endl;
        }
      }  //
    }    // j = 0 to 3

    if (isModified)
      return true;
    else {
      it->id() = 1;
    }
  }

  return false;
}

bool MedialAxisRepresentation::collapseSmallEdge(Poly &finalMesh, const K::FT &eps2,
  const std::vector<VertexProperty> &vertexTypes, const EdgeMap &targetEdges)
{
  for (Poly::Halfedge_handle hit = finalMesh.halfedges_begin(); hit != finalMesh.halfedges_end(); ++hit) {
    if ((hit->id() & 1) || (hit->opposite()->id() & 1))
      continue;

    if (hit->opposite()->is_border() || hit->is_border())
      continue;

    auto ret = computeCollapsingPt(hit, vertexTypes);
    if (std::get<0>(ret) == false)
      continue;

    K::FT d2 = CGAL::squared_distance(hit->vertex()->point(), hit->opposite()->vertex()->point());
    if (d2 < eps2) {
      // mark dirty
      Poly::Halfedge_handle h = hit;
      do {
        h->id() = 0;
        h->opposite()->id() = 0;

        h = h->next()->opposite();
      } while (h != hit);

      h = hit->opposite();
      do {
        h->id() = 0;
        h->opposite()->id() = 0;

        h = h->next()->opposite();
      } while (h != hit->opposite());

      K::Point_3 oldPositions[2] = {
        hit->vertex()->point(),
        hit->opposite()->vertex()->point()
      };

      hit->vertex()->point() = std::get<1>(ret);
      hit->opposite()->vertex()->point() = std::get<1>(ret);

      if (checkEdgeCollapsingFlip(finalMesh, hit, oldPositions, std::get<1>(ret)) == false) {
        continue;
      }

      hit->vertex()->id() = std::get<2>(ret);
      hit->opposite()->vertex()->id() = std::get<2>(ret);

      // perform collapsing
      if (collapseEdge(finalMesh, hit, vertexTypes, targetEdges)) {
        std::cout << "Collapse one edge: d2=" << CGAL::to_double(d2) << std::endl;
        return true;
      }
      else {
        std::cout << "Failed to collapse." << std::endl;
      }
    }
    else {
      hit->id() = 1;
      hit->opposite()->id() = 1;
    }
  }

  return false;
}

bool MedialAxisRepresentation::collapseColinearVertices(Poly &finalMesh, int &iter, std::set<Poly::Halfedge_handle> &skippedHE,
  const std::vector<VertexProperty> &vertexTypes,
  const EdgeMap &targetEdges)
{
  constexpr int debug = 0;

  for (auto hit = finalMesh.halfedges_begin(); hit != finalMesh.halfedges_end(); ++hit) {
    Poly::Halfedge_handle h0 = hit, h1 = hit;
    Poly::Halfedge_handle g0 = h0->opposite();
    Poly::Halfedge_handle g1 = g0;

    if (skippedHE.find(hit) != skippedHE.end())
      continue;

    if (hit->id() & 2)
      continue;

    // if both are constrained, we skip
    VertexAttachmentType vaType0 = vertexTypes[h0->vertex()->id()].vaType;
    VertexAttachmentType vaType1 = vertexTypes[g0->vertex()->id()].vaType;

    if (vaType0 == VertexAttachmentType::TARGET_VERTEX || vaType1 == VertexAttachmentType::TARGET_VERTEX)
      continue;

    K::Vector_3 dir_gh = h0->vertex()->point() - g0->vertex()->point();
    K::Ray_3 ray_gh(h0->vertex()->point(), h0->vertex()->point() + dir_gh);

    K::Vector_3 dir_hg = g0->vertex()->point() - h0->vertex()->point();
    K::Ray_3 ray_hg(g0->vertex()->point(), g0->vertex()->point() + dir_hg);

    bool found_h = false, found_g = false;
    do {
      h1 = h1->next();
      if (h1->opposite() != h0 && ray_gh.has_on(h1->vertex()->point())) {
        found_h = true;
        break;
      }

      h1 = h1->opposite();
    } while (h1 != h0);

    do {
      g1 = g1->next();
      if (g1->opposite() != g0 && ray_hg.has_on(g1->vertex()->point())) {
        found_g = true;
        break;
      }

      g1 = g1->opposite();
    } while (g1 != g0);

    if (found_g && found_h) {
      ES::V3d p0 = toVec3(h0->vertex()->point());
      ES::V3d p1 = toVec3(g0->vertex()->point());
      ES::V3d p2 = toVec3(h1->vertex()->point());
      ES::V3d p3 = toVec3(g1->vertex()->point());

      if constexpr (debug) {
        pgo::Mesh::TriMeshGeo mm;
        mm.addPos(p0);
        mm.addPos(p1);
        mm.addPos(p2);
        mm.addPos(p3);

        mm.save(fmt::format("temp/aa{:05d}.1.obj", iter));

        Poly::Halfedge_handle pq = h0;
        Poly::Halfedge_handle qp = pq->opposite();
        Poly::Halfedge_handle pt = pq->prev()->opposite();
        Poly::Halfedge_handle qb = qp->prev()->opposite();

        std::cout << "pt: " << pt->vertex()->degree() << std::endl;
        std::cout << "qb: " << qb->vertex()->degree() << std::endl;

        ES::V3d q = toVec3(pq->vertex()->point());
        ES::V3d p = toVec3(qp->vertex()->point());
        ES::V3d t = toVec3(pt->vertex()->point());
        ES::V3d b = toVec3(qb->vertex()->point());

        mm.addMesh(pgo::Mesh::createSingleTriangleMesh(p, q, p + pgo::asVec3d(1e-8)));
        mm.addMesh(pgo::Mesh::createSingleTriangleMesh(p, t, p + pgo::asVec3d(1e-8)));
        mm.addMesh(pgo::Mesh::createSingleTriangleMesh(q, b, q + pgo::asVec3d(1e-8)));
        mm.save(fmt::format("temp/aa{:05d}.2.obj", iter));
      }

      Poly::Halfedge_handle hs[4] = { h0, h1, g0, g1 };
      for (int i = 0; i < 4; i++) {
        Poly::Halfedge_handle hh = hs[i];
        do {
          hh->id() = 0;
          hh->opposite()->id() = 0;
          hh = hh->next()->opposite();
        } while (hh != hs[i]);
      }

      K::Point_3 mid_pt = CGAL::midpoint(h0->vertex()->point(), h0->opposite()->vertex()->point());
      K::Point_3 oldPositions[2] = {
        h0->vertex()->point(),
        h0->opposite()->vertex()->point()
      };

      h0->vertex()->point() = mid_pt;
      h0->opposite()->vertex()->point() = mid_pt;

      if (checkEdgeCollapsingFlip(finalMesh, h0, oldPositions, mid_pt) == false) {
        continue;
      }

      hit->opposite()->vertex()->id() = hit->vertex()->id();

      if (collapseEdge(finalMesh, h0, vertexTypes, targetEdges)) {
        // std::ofstream(fmt::format("temp/aa{:05d}.off", iter)) << finalMesh;
        // ++iter;
        // std::cout << "After Collapse: is triangle mesh: " << finalMesh.is_pure_triangle() << std::endl;
        std::cout << "Collapse one edge." << std::endl;

        return true;
      }
      else {
        std::cout << "Failed to collapse one edge." << std::endl;
        skippedHE.insert(h0);
      }
    }
    else {
      hit->id() |= 2;
    }
  }

  return false;
}

bool MedialAxisRepresentation::smoothVertices(Poly &finalMesh, std::vector<int> &isVertexChanged,
  const std::vector<VertexProperty> &vertexTypes, const EdgeMap &targetEdges)
{
// std::cout << finalMesh.size_of_vertices() << std::endl;
#if 0
  int modifiedCount = 0;
  for (auto hit = finalMesh.halfedges_begin(); hit != finalMesh.halfedges_end(); ++hit) {
    Poly::Halfedge_handle h0 = hit, h1 = hit;
    Poly::Halfedge_handle g0 = h0->opposite();

    // std::cout << h0->vertex()->id() << std::endl;

    if (isVertexChanged[h0->vertex()->id()])
      continue;

    if (vertexTypes[h0->vertex()->id()].vaType == VertexAttachmentType::TARGET_VERTEX)
      continue;

    if (vertexTypes[h0->vertex()->id()].isBorder)
      continue;

    // pgo::asVec3d z0(-0.115436, 0.225465, 0.0327274);
    // pgo::asVec3d z1 = toVec3(h0->vertex()->point());
    // if (len(z0 - z1) < 1e-3) {
    //   std::cout << "catched." << std::endl;
    // }

    K::Line_3 base_line(h0->vertex()->point(), g0->vertex()->point());
    bool found_h = false;
    do {
      h1 = h1->next();
      if (h1->opposite() != h0 && base_line.has_on(h1->vertex()->point())) {
        found_h = true;
        break;
      }

      h1 = h1->opposite();
    } while (h1 != h0);

    if (found_h) {
      // if (len(z0 - z1) < 1e-3) {
      //   std::cout << "catched." << std::endl;
      //   std::cout << toVec3(h0->vertex()->point()) << std::endl;
      //   std::cout << toVec3(g0->vertex()->point()) << std::endl;
      //   std::cout << toVec3(h1->vertex()->point()) << std::endl;
      // }

      if (coEdge(h0, g0, h1, vertexTypes) == false) {
        continue;
      }

      h0->vertex()->point() = CGAL::midpoint(g0->vertex()->point(), h1->vertex()->point());
      // isVertexChanged[h0->vertex()->id()] = 1;

      std::cout << "move one to midpoint of the edge." << std::endl;
      modifiedCount++;
    }
  }

  if constexpr (1) {
    TriMeshGeo tempM;
    pgo::CGALInterface::polyhedron2TriangleMesh(finalMesh, tempM);
    tempM.save("tc1.obj");
  }
#endif

  for (auto vit = finalMesh.vertices_begin(); vit != finalMesh.vertices_end(); ++vit) {
    int vid = (int)vit->id();

    if (isVertexChanged[vid])
      continue;

    if (isFixed(vertexTypes[vid].vaType))
      continue;

    if (vertexTypes[vid].isBorder == 0)
      continue;

    Poly::Halfedge_handle h = vit->halfedge();
    Poly::Halfedge_handle g0, g1;
    Poly::Halfedge_handle h1 = h;

    bool found = false;
    do {
      int vid1 = h1->opposite()->vertex()->id();
      if (vertexTypes[vid1].isBorder) {
        g0 = h1->opposite();
        found = true;
        break;
      }

      h1 = h1->next()->opposite();
    } while (h1 != h);

    PGO_ALOG(found == true);
    h1 = h1->next()->opposite();

    found = false;
    do {
      int vid1 = h1->opposite()->vertex()->id();
      if (vertexTypes[vid1].isBorder) {
        g1 = h1->opposite();
        found = true;
        break;
      }

      h1 = h1->next()->opposite();
    } while (h1 != h);

    if (found == false) {
      Poly p;
      p.make_triangle(g0->vertex()->point(), h->vertex()->point(), g0->vertex()->point());
      std::ofstream("aa.off") << p;
    }
    else {
      vit->point() = CGAL::midpoint(g0->vertex()->point(), g1->vertex()->point());
    }
  }

  if constexpr (0) {
    pgo::Mesh::TriMeshGeo tempM;
    pgo::CGALInterface::polyhedron2TriangleMesh(finalMesh, tempM);
    tempM.save("tc2.obj");
  }

  return true;
}

void MedialAxisRepresentation::cleanUpMesh(const Poly &targetMeshPoly, const FaceGraphTree &targetMeshBVTreeExact,
  Poly &finalMesh, double edgeLengthThreshold,
  std::vector<VertexProperty> &vertexTypes,
  const EdgeMap &targetEdges)
{
  constexpr int debug = 0;

  edgeLengthThreshold = 1e-3;

  // ========================================
  // mark feature vertices
  K::FT eps = 1e-10;
  K::FT eps2 = eps * eps;

#if 0
  // vtx handles
  tbb::concurrent_hash_map<Poly::Vertex_handle, Poly::Vertex_handle> thisMeshToTargetMesh;  
  std::vector<Poly::Vertex_handle> vtxHandles;
  vtxHandles.reserve(finalMesh.size_of_vertices());
  for (auto vit = finalMesh.vertices_begin(); vit != finalMesh.vertices_end(); ++vit) {
    vit->id() = 0;
    vtxHandles.emplace_back(vit);
  }

  tbb::parallel_for(0, (int)vtxHandles.size(), [&](int vi) {
    Poly::Vertex_handle vit = vtxHandles[vi];

    // query the closest point on the target mesh
    auto ret = targetMeshBVTreeExact.closest_point_and_primitive(vit->point());
    K::FT min_d2(1e100);

    Poly::Halfedge_handle h0 = ret.second->halfedge();
    Poly::Halfedge_handle h1 = h0, hmin;
    do {
      K::FT d2 = CGAL::squared_distance(h1->vertex()->point(), vit->point());
      if (d2 < min_d2) {
        min_d2 = d2;
        hmin = h1;
      }
      h1 = h1->next();
    } while (h1 != h0);

    // if the closest point is a vertex, we don't change it
    if (min_d2 < eps2) {
      std::cout << "g0 is a feature. d2=" << CGAL::to_double(min_d2) << std::endl;
      vit->id() = 1;
      thisMeshToTargetMesh.emplace(vit, hmin->vertex());
    }
  });
#endif

  // ========================================
  // remove small pieces
  eps = edgeLengthThreshold;
  eps2 = eps * eps;

  // int r = SMS::edge_collapse(finalMesh,
  //   CGAL::Surface_mesh_simplification::Edge_length_stop_predicate<K::FT>(edgeLengthThreshold),
  //   CGAL::parameters::get_cost(SMS::Edge_length_cost<Poly>())
  //     .get_placement(SMS::Midpoint_placement<Poly>()));

  // std::cout << "#collapsed edges: " << r << std::endl;

  for (Poly::Halfedge_handle hit = finalMesh.halfedges_begin(); hit != finalMesh.halfedges_end(); ++hit) {
    hit->id() = 0;
  }

  for (Poly::Facet_handle fit = finalMesh.facets_begin(); fit != finalMesh.facets_end(); ++fit) {
    fit->id() = 0;
  }

  int iter = 0;
  std::set<Poly::Halfedge_handle> skippedHE;
  while (1) {
    // if constexpr (1) {
    //   TriMeshGeo tempM;
    //   pgo::CGALInterface::polyhedron2TriangleMesh(finalMesh, tempM);
    //   tempM.save("tc.obj");
    // }

    // std::cout << "Is triangle mesh: " << finalMesh.is_pure_triangle() << std::endl;

    // perform one edge collapsing
    if (collapseSmallEdge(finalMesh, eps2, vertexTypes, targetEdges)) {
#if 0
      for (Poly::Facet_handle fi = finalMesh.facets_begin(); fi != finalMesh.facets_end(); ++fi) {
        Poly::Halfedge_handle h0 = fi->halfedge();
        Poly::Halfedge_handle h1 = h0->next();
        Poly::Halfedge_handle h2 = h1->next();

        K::Vector_3 e01 = h1->vertex()->point() - h0->vertex()->point();
        K::Vector_3 e02 = h2->vertex()->point() - h0->vertex()->point();
        K::Vector_3 n = CGAL::cross_product(e01, e02);

        Poly::Halfedge_handle g[3] = {
          h0->opposite(), h1->opposite(), h2->opposite()
        };

        K::Vector_3 n1[3];
        for (int j = 0; j < 3; j++) {
          Poly::Halfedge_handle r0 = g[j];
          Poly::Halfedge_handle r1 = r0->next();
          Poly::Halfedge_handle r2 = r1->next();

          K::Vector_3 e01 = r1->vertex()->point() - r0->vertex()->point();
          K::Vector_3 e02 = r2->vertex()->point() - r0->vertex()->point();
          n1[j] = CGAL::cross_product(e01, e02);
        }

        pgo::asVec3d n_d = toVec3(n);
        pgo::asVec3d n_d0 = toVec3(n1[0]);
        pgo::asVec3d n_d1 = toVec3(n1[1]);
        pgo::asVec3d n_d2 = toVec3(n1[2]);
        double l = len(n_d);
        double l1[3] = { len(n_d0), len(n_d1), len(n_d2) };

        if (l < 1e-16 || l1[0] < 1e-6 || l1[1] < 1e-6 || l1[2] < 1e-6)
          continue;

        n_d /= l;
        n_d0 /= l1[0];
        n_d1 /= l1[1];
        n_d2 /= l1[2];

        if (dot(n_d, n_d0) < -0.8 ||
          dot(n_d, n_d1) < -0.8||
          dot(n_d, n_d2) < -0.8) {
          Poly p;
          p.make_triangle(h0->vertex()->point(), h1->vertex()->point(), h2->vertex()->point());
          std::ofstream("zzzz.off") << p;
          std::cout << "zzz";
        }
      }
#endif

      continue;
    }

    // if constexpr (1) {
    //   TriMeshGeo tempM;
    //   pgo::CGALInterface::polyhedron2TriangleMesh(finalMesh, tempM);
    //   tempM.save("tc.obj");
    // }

    // std::cout << "Is triangle mesh: " << finalMesh.is_pure_triangle() << std::endl;

    // perform one triangle collapsing
    if (collapseSmallTriangleDouble(finalMesh, edgeLengthThreshold * edgeLengthThreshold, vertexTypes, targetEdges))
      continue;

    // if constexpr (1) {
    //   TriMeshGeo tempM;
    //   pgo::CGALInterface::polyhedron2TriangleMesh(finalMesh, tempM);
    //   tempM.save("tc.obj");
    // }

    // int z = 0;
    // for (Poly::Halfedge_handle hit = finalMesh.halfedges_begin(); hit != finalMesh.halfedges_end(); ++hit) {
    //   if (hit->id() == 1)
    //     z++;
    // }
    // std::cout << z << ',' << finalMesh.size_of_halfedges() << std::endl;
    // std::cout << "Is triangle mesh: " << finalMesh.is_pure_triangle() << std::endl;

    // perform one colinear collapsing
    if (collapseColinearVertices(finalMesh, iter, skippedHE, vertexTypes, targetEdges)) {
      continue;
    }

    break;
  }

  std::cout << "Done." << std::endl;

  if constexpr (0) {
    pgo::Mesh::TriMeshGeo tempM;
    pgo::CGALInterface::polyhedron2TriangleMesh(finalMesh, tempM);
    tempM.save("tc.obj");
  }

  int maxvid = 0;
  for (auto vit = finalMesh.vertices_begin(); vit != finalMesh.vertices_end(); ++vit) {
    maxvid = std::max(maxvid, (int)vit->id());
  }

  std::vector<int> isVtxChanged(maxvid + 1, 0);
  for (int i = 0; i < 2; i++) {
    smoothVertices(finalMesh, isVtxChanged, vertexTypes, targetEdges);
  }

  if constexpr (0) {
    pgo::Mesh::TriMeshGeo tempM;
    pgo::CGALInterface::polyhedron2TriangleMesh(finalMesh, tempM);
    tempM.save("tc3.obj");
  }

  //  for (auto vit = finalMesh.vertices_begin(); vit != finalMesh.vertices_end(); ++vit) {
  //    vit->id() = 0;
  //
  //    // query the closest point on the target mesh
  //    auto ret = targetMeshBVTreeExact.closest_point_and_primitive(vit->point());
  //    K::FT min_d2 = CGAL::squared_distance(ret.first, vit->point());
  //
  //    // if the closest point is a vertex, we don't change it
  //    if (min_d2 < eps2) {
  //      std::cout << "g0 is a feature: " << CGAL::to_double(min_d2) << std::endl;
  //      vit->id() = 1;
  //      continue;
  //    }
  //  }
}

void MedialAxisRepresentation::orientMesh(const std::vector<K::Point_3> &vertices, std::vector<std::vector<int>> &faces, const K::Point_3 centers[3], int type)
{
  if (type == 0) {
  }
  else if (type == 1) {
    K::Line_3 line(centers[0], centers[1]);

    for (int ti = 0; ti < (int)faces.size(); ti++) {
      std::vector<int> &f = faces[ti];
      K::FT val[3] = { 0, 0, 0 };
      for (int j = 0; j < (int)f.size(); j++) {
        val[0] += vertices[f[j]][0];
        val[1] += vertices[f[j]][1];
        val[2] += vertices[f[j]][2];
      }
      K::Point_3 triCenter(val[0] / double(f.size()), val[1] / double(f.size()), val[2] / double(f.size()));
      K::Vector_3 dir = centers[1] - centers[0];
      K::Vector_3 diff = triCenter - centers[0];
      // dT (t d + p0 - x) = 0
      K::FT t = CGAL::scalar_product(diff, dir) / dir.squared_length();
      K::Point_3 cpt = centers[0] + dir * t;
      dir = triCenter - cpt;

      // n
      K::Vector_3 e0 = vertices[f[1]] - vertices[f[0]];
      K::Vector_3 e1 = vertices[f[2]] - vertices[f[0]];
      K::Line_3 e0Line(vertices[f[0]], vertices[f[1]]);
      bool found = false;
      for (int j = 2; j < (int)faces.size(); j++) {
        if (e0Line.has_on(vertices[f[j]]) == false) {
          e1 = vertices[f[j]] - vertices[f[0]];
          found = true;
          break;
        }
      }
      PGO_ALOG(found == true);

      K::Vector_3 n = CGAL::cross_product(e0, e1);
      if (CGAL::scalar_product(n, dir) < 0) {
        std::reverse(f.begin(), f.end());
      }
    }
  }
  else if (type == 2) {
  }
  else {
    throw std::invalid_argument("unknown type");
  }

  if constexpr (0) {
    std::ofstream outfile("mm.obj");

    for (const auto &p : vertices) {
      ES::V3d a = ES::Mp<const ES::V3d>(pgo::CGALInterface::toDoublePt3(p).data());
      outfile << "v " << a[0] << ' ' << a[1] << ' ' << a[2] << std::endl;
    }

    for (int i = 0; i < (int)faces.size(); i++) {
      outfile << "f";
      for (int j = 0; j < (int)faces[i].size(); j++) {
        outfile << " " << faces[i][j] + 1;
      }
      outfile << std::endl;
    }
    outfile.close();
  }

  exit(1);

  /*
  std::unordered_map<std::pair<int, int>, std::array<int, 5>, VegaFEM::EigenSupport::IntPairHash, VegaFEM::EigenSupport::IntPairEqual> edgeTriangles;
  for (int ti = 0; ti < (int)triangles.size(); ti++) {
    for (int j = 0; j < 3; j++) {
      std::pair<int, int> edgeID{ triangles[ti][j], triangles[ti][(j + 1) % 3] };
      if (edgeID.first > edgeID.second)
        std::swap(edgeID.first, edgeID.second);

      auto it = edgeTriangles.find(edgeID);
      if (it != edgeTriangles.end()) {
        int &count = it->second.back();
        PGO_ALOG(count < 2);

        it->second[count++] = ti;
      }
      else {
        std::array<int, 5> triIDs = { 0, 0, 0, 0, 0 };
        triIDs[0] = ti;
        triIDs.back() = 1;

        edgeTriangles.emplace(edgeID, triIDs);
      }
    }
  }

  std::vector<std::array<int, 5>> triangleNeighbors(triangles.size(), std::array<int, 5>{ 0, 0, 0, 0, 0 });
  for (auto it = edgeTriangles.begin(); it != edgeTriangles.end(); ++it) {
    PGO_ALOG(it->second.back() == 2);

    // t0
    int &count = triangleNeighbors[it->second[0]].back();
    PGO_ALOG(count < 3);
    triangleNeighbors[it->second[0]][count++] = it->second[1];

    // t1
    int &count1 = triangleNeighbors[it->second[1]].back();
    PGO_ALOG(count1 < 3);
    triangleNeighbors[it->second[1]][count1++] = it->second[0];
  }

  std::vector<int> isTriVisited(triangles.size(), 0);
  std::vector<int> triQ;
  size_t start = 0;

  triQ.emplace_back(0);
  while (triQ.size() > start) {
    int curTriID = triQ[start++];

    for (int j = 0; j < triangleNeighbors[curTriID].back(); j++) {
      int nTriID = triangleNeighbors[curTriID][j];
      if (isTriVisited[nTriID])
        continue;

      std::array<int, 3> reverseIdx = { -1, -1, -1 };
      for (int k = 0; k < 3; k++) {
        for (int r = 0; r < 3; r++) {
          if (triangles[curTriID][k] == triangles[nTriID][r]) {
            reverseIdx[k] = r;
            break;
          }
        }
      }

      for (int k = 0; k < 3; k++) {
        if (reverseIdx[k] >= 0 &&
          reverseIdx[(k + 1) % 3] >= 0) {
          if (reverseIdx[k] > reverseIdx[(k + 1) % 3]) {
            std::cout << "????" << std::endl;
          }
          else {
            std::swap(triangles[nTriID][1], triangles[nTriID][2]);
          }
        }
      }

      triQ.emplace_back(nTriID);
    }
  }

  if constexpr (1) {
    TriMeshGeo zza;
    for (const auto &p : vertices) {
      zza.addPos(pgo::asVec3d(pgo::CGALInterface::toDoublePt3(p).data()));
    }
    for (int i = 0; i < (int)triangles.size(); i++) {
      if (1) {
        zza.addTri(pgo::asVec3i(triangles[i].data()));
      }
    }
    zza.save("mm.obj");
  }
  exit(1);
  */
}

namespace MedialAxisRepresentation
{
class PolyhedronBuilderFace : public CGAL::Modifier_base<Poly::HalfedgeDS>
{
protected:
  const std::vector<K::Point_3> &vertices;
  std::vector<std::vector<int>> &faces;

public:
  PolyhedronBuilderFace(const std::vector<K::Point_3> &vtx, std::vector<std::vector<int>> &faces_):
    vertices(vtx), faces(faces_) {}

  using HalfedgeDS = Poly::HalfedgeDS;
  void operator()(HalfedgeDS &hds)
  {
    typedef typename HalfedgeDS::Vertex Vertex;
    typedef typename Vertex::Point Point;

    std::set<std::pair<int, int>> visitedEdges;
    std::vector<int> isFaceVisited(faces.size(), 0);

    // create a cgal incremental builder
    CGAL::Polyhedron_incremental_builder_3<HalfedgeDS> B(hds, true);
    B.begin_surface(vertices.size(), faces.size());

    // add the polyhedron vertices
    for (int vi = 0; vi < (int)vertices.size(); vi++) {
      HalfedgeDS::Vertex_handle vit = B.add_vertex(vertices[vi]);
      vit->id() = vi;
    }

    while (1) {
      int count = 0;

      for (int fi = 0; fi < (int)faces.size(); fi++) {
        // if face is added
        if (isFaceVisited[fi])
          continue;

        // if the face is not shared edges
        // we temporarily skip
        int visitedEdgeCount = 0;
        for (int j = 0; j < (int)faces[fi].size(); j++) {
          std::pair edgeID{ faces[fi][j], faces[fi][(j + 1) % faces[fi].size()] };
          if (edgeID.first > edgeID.second)
            std::swap(edgeID.first, edgeID.second);

          if (visitedEdges.find(edgeID) != visitedEdges.end()) {
            visitedEdgeCount++;
          }
        }

        if (visitedEdges.size() && visitedEdgeCount == 0)
          continue;

        if (B.test_facet(faces[fi].begin(), faces[fi].end()) == false) {
          std::reverse(faces[fi].begin(), faces[fi].end());
        }

        if (B.test_facet(faces[fi].begin(), faces[fi].end()) == false) {
          std::cout << "Cannot add facet " << fi << std::endl;
        }
        else {
          HalfedgeDS::Halfedge_handle hit = B.add_facet(faces[fi].begin(), faces[fi].end());
          hit->facet()->id() = fi;

          for (int j = 0; j < (int)faces[fi].size(); j++) {
            std::pair edgeID{ faces[fi][j], faces[fi][(j + 1) % faces[fi].size()] };
            if (edgeID.first > edgeID.second)
              std::swap(edgeID.first, edgeID.second);

            visitedEdges.emplace(edgeID);
          }

          isFaceVisited[fi] = 1;
          count++;
        }
      }

      if (count == 0)
        break;
    }

    B.end_surface();
  }
};
}  // namespace MedialAxisRepresentation

void MedialAxisRepresentation::orientMesh(const std::vector<K::Point_3> &vertices, std::vector<std::vector<int>> &faces, Poly &finalMesh)
{
  finalMesh.clear();

  PolyhedronBuilderFace localBuilder(vertices, faces);
  finalMesh.delegate(localBuilder);

  if constexpr (0) {
    std::cout << finalMesh.size_of_facets() << ',' << faces.size() << std::endl;

    std::ofstream outfile("zzzz.obj");
    for (int vi = 0; vi < (int)vertices.size(); vi++) {
      ES::V3d p = toVec3(vertices[vi]);
      outfile << "v " << p[0] << ' ' << p[1] << ' ' << p[2] << std::endl;
    }

    for (int fi = 0; fi < (int)faces.size(); fi++) {
      outfile << "f";
      for (int f : faces[fi]) {
        outfile << " " << (f + 1);
      }
      outfile << std::endl;
    }
    outfile.close();
  }
}

bool MedialAxisRepresentation::triangulateMesh(Poly &finalMesh)
{
  std::vector<K::Point_3> polylines;
  std::vector<int> vertexIDs;
  polylines.reserve(1000);
  vertexIDs.reserve(1000);

  std::vector<CGAL::Triple<int, int, int>> triangles;
  triangles.reserve(1000);

  std::list<Poly::Halfedge_handle> currentFacets;

  // std::ofstream("tri0.off") << finalMesh;

  std::unordered_set<std::pair<int, int>, ES::IntPairHash, ES::IntPairEqual> visitedEdges;
  int baseVtxID = (int)finalMesh.size_of_vertices();
  for (auto fi = finalMesh.facets_begin(); fi != finalMesh.facets_end(); ++fi) {
    if (fi->is_triangle())
      continue;

    polylines.clear();
    vertexIDs.clear();
    triangles.clear();

    Poly::Halfedge_handle h0 = fi->halfedge();
    Poly::Halfedge_handle h1 = h0;

    do {
      polylines.emplace_back(h1->vertex()->point());
      vertexIDs.emplace_back(h1->vertex()->id());
      h1 = h1->next();
    } while (h1 != h0);

    CGAL::Polygon_mesh_processing::triangulate_hole_polyline(polylines, std::back_inserter(triangles));

    if (triangles.size() == 0ull)
      continue;

    currentFacets.clear();
    currentFacets.emplace_back(h0);

    for (const auto &tri : triangles) {
      int vids[3] = {
        tri.first,
        tri.second,
        tri.third
      };

      for (int j = 0; j < 3; j++) {
        if ((vids[j] == ((vids[(j + 1) % 3]) + 1) % (int)polylines.size()) ||
          ((vids[j] + 1) % (int)polylines.size() == (vids[(j + 1) % 3]))) {
          continue;
        }

        std::pair<int, int> edgeID{ vertexIDs[vids[j]], vertexIDs[vids[(j + 1) % 3]] };
        if (edgeID.first > edgeID.second) {
          std::swap(edgeID.first, edgeID.second);
        }

        if (edgeID.first == 1151 && edgeID.second == 321) {
          std::cout << "1" << std::endl;
        }

        if (edgeID.second == 1151 && edgeID.first == 321) {
          std::cout << "2" << std::endl;
        }

        if (visitedEdges.find(edgeID) != visitedEdges.end()) {
          std::cout << "Edge (" << edgeID.first << ',' << edgeID.second << ") has been visited. " << std::endl;
        }
        else {
          std::list<Poly::Halfedge_handle>::iterator selected_face = currentFacets.end();
          Poly::Halfedge_handle h_pair[2];

          for (auto it = currentFacets.begin(); it != currentFacets.end(); ++it) {
            Poly::Halfedge_handle hh0 = *it;
            Poly::Halfedge_handle hh1 = hh0;
            bool found[2] = { false, false };

            do {
              if (hh1->vertex()->id() == edgeID.first) {
                found[0] = true;
                h_pair[0] = hh1;
              }

              if (hh1->vertex()->id() == edgeID.second) {
                found[1] = true;
                h_pair[1] = hh1;
              }
              hh1 = hh1->next();
            } while (hh1 != hh0);

            if (found[0] && found[1]) {
              selected_face = it;
              break;
            }
          }

          PGO_ALOG(selected_face != currentFacets.end());
          Poly::Halfedge_handle hret = finalMesh.split_facet(h_pair[0], h_pair[1]);

          currentFacets.erase(selected_face);
          currentFacets.emplace_back(hret);
          currentFacets.emplace_back(hret->opposite());

          visitedEdges.emplace(edgeID);
        }
      }
    }
  }

  return finalMesh.is_pure_triangle();
}

#include <nlohmann/json.hpp>

void MedialAxisRepresentation::saveExact(Poly &finalMesh, const char *filename)
{
  nlohmann::json meshOutJson;

  int idx = 0;
  std::map<Poly::Vertex_handle, int> vertexIndices;
  for (Poly::Vertex_iterator itt = finalMesh.vertices_begin(); itt != finalMesh.vertices_end(); ++itt) {
    K::Point_3 &pt = itt->point();

    std::vector<char> val[3];
    for (int j = 0; j < 3; j++) {
      val[j].resize((mpz_sizeinbase(mpq_numref(pt[j].get_mpq_t()), 10) + mpz_sizeinbase(mpq_denref(pt[j].get_mpq_t()), 10) + 3) * 2);
      memset(val[j].data(), 0, sizeof(char) * val[j].size());
      mpq_get_str(val[j].data(), 10, pt[j].get_mpq_t());
    }

    meshOutJson[fmt::format("v{:05d}", idx)][0] = val[0].data();
    meshOutJson[fmt::format("v{:05d}", idx)][1] = val[1].data();
    meshOutJson[fmt::format("v{:05d}", idx)][2] = val[2].data();

    vertexIndices.insert(make_pair(itt, idx));
    idx++;
  }

  idx = 0;
  for (Poly::Facet_iterator itt = finalMesh.facets_begin(); itt != finalMesh.facets_end(); ++itt) {
    std::array<int, 3> face;
    int inc = 0;

    Poly::Halfedge_around_facet_circulator ctr = itt->facet_begin();
    do {
      auto idxItt = vertexIndices.find(ctr->vertex());
      if (idxItt == vertexIndices.end()) {
        std::cerr << "Cannot find the vertex" << std::endl;
        abort();
      }
      if (inc == 3) {
        std::cerr << "Not a triangle mesh" << std::endl;
        abort();
      }

      face[inc++] = idxItt->second;
      ctr++;
    } while (ctr != itt->facet_begin());

    meshOutJson[fmt::format("f{:05d}", idx)] = face;
    idx++;
  }

  std::ofstream(filename) << meshOutJson.dump(2);
}

#if 0
bool MedialAxisRepresentation::collapseSmallTriangle(Poly &finalMesh, const K::FT &eps2)
{
  constexpr int debug = 0;
  static int iii = 0;

  for (auto it : finalMesh.facet_handles()) {
    std::cout << "f ";

    if (it->id() & 1) {
      std::cout << "done" << std::endl;
      continue;
    }
    else {
      std::cout << "not done" << std::endl;
    }

    Poly::Halfedge_handle h[3] = { it->halfedge() };
    h[1] = h[0]->next();
    h[2] = h[1]->next();

    bool isModified = false;

    for (int j = 0; j < 3; j++) {
      K::Line_3 line(h[(j + 1) % 3]->vertex()->point(), h[(j + 2) % 3]->vertex()->point());
      K::FT d2 = CGAL::squared_distance(h[j]->vertex()->point(), line);
      if (d2 > eps2) {
        continue;
      }

      std::cout << "Find one bad triangle: " << CGAL::to_double(d2) << std::endl;

      // check intersection point location, if it is on the segment or not.
      // we only want the ones with it on
      // (t (p1 - p0) + p0 - x)T (p1 - p0)
      // t (p1- p0)^T(p1 - p0
      K::Vector_3 dir = h[(j + 2) % 3]->vertex()->point() - h[(j + 1) % 3]->vertex()->point();
      K::Vector_3 diff = h[(j + 1) % 3]->vertex()->point() - h[j]->vertex()->point();

      K::FT a = dir.squared_length();
      K::FT b = CGAL::scalar_product(dir, diff);
      K::FT t = -b / a;

      if (t <= 0 || t >= 1) {
        std::cout << "t not in [0, 1]; t=" << CGAL::to_double(t) << std::endl;
        continue;
      }
      else {
        std::cout << "t in [0, 1]; t=" << CGAL::to_double(t) << std::endl;
      }

      if constexpr (1) {
        //        Poly pp;
        // Poly::Halfedge_handle h0 = pp.make_triangle(h[j]->vertex()->point(), h[(j + 1) % 3]->vertex()->point(), h[(j + 2) % 3]->vertex()->point());
        // Poly::Halfedge_handle hh = pp.add_vertex_and_facet_to_border(h0->opposite(), h0->prev()->opposite());
        Poly::Halfedge_handle h_oppo = h[(j + 2) % 3]->opposite()->next();
        //        hh->vertex()->point() = h_oppo->vertex()->point();
        //
        //        std::ofstream("aa1.off") << pp;
        static int iii = 0;
        TriMeshGeo aar;
        pgo::asVec3d p0 = toVec3(h[j]->vertex()->point());
        pgo::asVec3d p1 = toVec3(h[(j + 1) % 3]->vertex()->point());
        pgo::asVec3d p2 = toVec3(h[(j + 2) % 3]->vertex()->point());
        pgo::asVec3d p3 = toVec3(h_oppo->vertex()->point());

        aar.addPos(p0);
        aar.addPos(p1);
        aar.addPos(p2);
        aar.addPos(p3);

        aar.addTri(pgo::asVec3i(0, 1, 2));
        aar.addTri(pgo::asVec3i(2, 1, 3));

        aar.save(fmt::format("aa{}.obj", iii));

        iii = (iii + 1) % 10;
      }

      // new vtx pos
      K::Point_3 newPt = h[(j + 1) % 3]->vertex()->point() + dir * t;
      bool isMovingValid = true;

      // check if the current vertex is moved to the target location,
      // all related triangles has small area
      Poly ::Halfedge_handle hh = h[j]->next()->opposite();
      do {
        Poly::Halfedge_handle hv0 = hh;
        Poly::Halfedge_handle hv1 = hv0->next();
        Poly::Halfedge_handle hv2 = hv1->next();

        if (hv0->is_border()) {
        }
        else {
          K::Point_3 v_new[3] = { newPt, hv1->vertex()->point(), hv2->vertex()->point() };
          K::Point_3 v_old[3] = { hv0->vertex()->point(), hv1->vertex()->point(), hv2->vertex()->point() };
          K::Triangle_3 tri_new(v_new[0], v_new[1], v_new[2]);
          K::FT area2_new = tri_new.squared_area() * 4;

          K::Triangle_3 tri_old(v_old[0], v_old[1], v_old[2]);
          K::FT area2_old = tri_old.squared_area() * 4;

          // check old
          bool isOldSmall = false;
          for (int k = 0; k < 3; k++) {
            K::FT d2 = CGAL::squared_distance(v_old[k], v_old[(k + 1) % 3]);
            K::FT h2 = area2_old / d2;
            if (d2 == 0) {
              // std::cout << "old" << CGAL::to_double(d2) << ',' << CGAL::to_double(area2_old) << std::endl;
              isOldSmall = true;
              break;
            }

            if (h2 < eps2) {
              isOldSmall = true;
              break;
            }
          }

          // check new
          bool isNewSmall = false;
          for (int k = 0; k < 3; k++) {
            K::FT d2 = CGAL::squared_distance(v_new[k], v_new[(k + 1) % 3]);
            if (d2 == 0) {
              // std::cout << "new" << CGAL::to_double(d2) << ',' << CGAL::to_double(area2_new) << std::endl;
              // Poly pp;
              // Poly::Halfedge_handle h0 = pp.make_triangle(v_new[0], v_new[1], v_new[2]);
              // std::ofstream("aaaaa.off") << pp;

              // pp.clear();
              // pp.make_triangle(v_old[0], v_old[1], v_old[2]);
              // std::ofstream("aaaaa1.off") << pp;

              isNewSmall = true;
              break;
            }
            K::FT h2 = area2_new / d2;

            if (h2 < eps2) {
              isNewSmall = true;
              break;
            }
          }
          // if the old triangle is not small but the new one is small
          if (isOldSmall == false && isNewSmall) {
            isMovingValid = false;
            break;
          }
        }

        hh = hh->next()->opposite();
      } while (hh != h[j]);

      if (isMovingValid == false) {
        continue;
      }

      // check opposite triangle, if the area is small too, then we skip this temporarily
      Poly::Halfedge_handle h_oppo = h[(j + 2) % 3]->opposite()->next();
      // if the opposite is a border
      if (h_oppo->is_border()) {
        // create a point in the middle of the segment
        Poly::Halfedge_handle hnew = finalMesh.split_edge(h[(j + 2) % 3]);
        hnew->vertex()->point() = newPt;

        // create two faces
        Poly::Halfedge_handle h1 = finalMesh.split_facet(hnew, h[j]);

        // mark dirty
        Poly ::Halfedge_handle hh = h[j];
        do {
          if (hh->is_border() == false)
            hh->facet()->id() = 0;

          hh = hh->next()->opposite();
        } while (hh != h[j]);

        /*
                if (collapseEdge(finalMesh, h1)) {
                  // std::cout << "After splitting: is triangle mesh: " << finalMesh.is_pure_triangle() << std::endl;
                  std::cout << "Triangle collapsed. " << std::endl;

                  //  TriMeshGeo mm;
                  //  pgo::CGALInterface::polyhedron2TriangleMesh(finalMesh, mm);

                  isModified = true;
                  break;
                }
                else {
                  std::cout << "Failed to collapse." << std::endl;
                }
                */
      }
      // otherwise
      else {
        // oppo triangle area all:
        d2 = CGAL::squared_distance(h_oppo->vertex()->point(), line);
        K::FT oppoArea2 = d2 * dir.squared_length();

        // check if all new opposite triangles are big enough
        K::Vector_3 segs[] = {
          newPt - h[(j + 1) % 3]->vertex()->point(),
          h[(j + 2) % 3]->vertex()->point() - newPt,
          h_oppo->vertex()->point() - newPt,
          h_oppo->vertex()->point() - h[(j + 1) % 3]->vertex()->point(),
          h_oppo->vertex()->point() - h[(j + 2) % 3]->vertex()->point()
        };

        bool isOppoValid = true;
        for (int i = 0; i < 5; i++) {
          K::FT h2 = oppoArea2 / segs->squared_length();
          if (h2 < eps2) {
            std::cout << "Opposite triangle is too small: h2=" << CGAL::to_double(h2) << std::endl;
            isOppoValid = false;
            break;
          }
          else {
            std::cout << "Opposite triangle h2=" << CGAL::to_double(h2) << std::endl;
          }
        }

        if (isOppoValid == false) {
          continue;
        }

        if constexpr (0) {
          K::FT edgeLength[3] = {
            CGAL::squared_distance(h[0]->vertex()->point(), h[0]->opposite()->vertex()->point()),
            CGAL::squared_distance(h[1]->vertex()->point(), h[1]->opposite()->vertex()->point()),
            CGAL::squared_distance(h[2]->vertex()->point(), h[2]->opposite()->vertex()->point()),
          };

          if (edgeLength[0] < edgeLength[1] && edgeLength[0] < edgeLength[2]) {
            CGAL::Euler::collapse_edge(CGAL::edge(h[0], finalMesh), finalMesh);
          }
          else if (edgeLength[1] < edgeLength[0] && edgeLength[1] < edgeLength[2]) {
            CGAL::Euler::collapse_edge(CGAL::edge(h[1], finalMesh), finalMesh);
          }
          else {
            CGAL::Euler::collapse_edge(CGAL::edge(h[2], finalMesh), finalMesh);
          }

          std::cout << "After splitting: is triangle mesh: " << finalMesh.is_pure_triangle() << std::endl;

          isModified = true;
          break;
        }
        else {
          // create a point in the middle of the segment
          Poly::Halfedge_handle hnew = finalMesh.split_edge(h[(j + 2) % 3]);
          Poly::Halfedge_handle hnew_oppo = hnew->next()->opposite();

          // std::cout << hnew->vertex()->point() << std::endl;
          // std::cout << hnew->next()->vertex()->point() << std::endl;
          // std::cout << hnew->prev()->vertex()->point() << std::endl;
          // std::cout << "==\n";

          // K::Point_3 pnew = hnew->vertex()->point() + dir * t;
          hnew->vertex()->point() = newPt;

          // std::cout << hnew->vertex()->point() << std::endl;
          // std::cout << hnew->next()->vertex()->point() << std::endl;
          // std::cout << hnew->prev()->vertex()->point() << std::endl;
          Poly::Halfedge_handle h1 = finalMesh.split_facet(hnew, h[j]);

          h[j]->vertex()->point() = newPt;

          // mark dirty
          Poly ::Halfedge_handle hh = h[j];
          do {
            hh->facet()->id() = 0;
            hh = hh->next()->opposite();
          } while (hh != h[j]);

          // Poly pp1;
          // pp1.make_triangle(h1->vertex()->point(), h1->next()->vertex()->point(), h1->next()->next()->vertex()->point());
          // std::ofstream("aa.off") << pp1;

          Poly::Halfedge_handle g = hnew_oppo->next()->next();
          Poly::Halfedge_handle h2 = finalMesh.split_facet(hnew_oppo, g);

          if constexpr (debug) {
            Poly pp;
            pp.make_triangle(h[0]->vertex()->point(), h[1]->vertex()->point(), h[2]->vertex()->point());
            pp.make_triangle(h[0]->vertex()->point(), h[1]->vertex()->point(), h[2]->vertex()->point());
            std::ofstream("aa1.off") << pp;
          }

          // mark dirty
          h2->facet()->id() = 0;
          h2->opposite()->facet()->id() = 0;

          /*
          if (collapseEdge(finalMesh, h1)) {
            // std::cout << "After splitting: is triangle mesh: " << finalMesh.is_pure_triangle() << std::endl;
            std::cout << "Triangle collapsed. " << std::endl;

            TriMeshGeo mm;
            pgo::CGALInterface::polyhedron2TriangleMesh(finalMesh, mm);
            iii++;
            if (iii > 10) {
              mm.save("rr.obj");
              iii = 0;
            }

            isModified = true;
            break;
          }
          else {
            std::cout << "Failed to collapse." << std::endl;
          }
          */
        }
      }  //
    }    // j = 0 to 3

    if (isModified)
      return true;
    else {
      it->id() = 1;
    }
  }

  return false;
}
#endif