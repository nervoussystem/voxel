#pragma once

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Delaunay_triangulation_3.h>
#include "CGAL/Triangulation_vertex_base_with_info_3.h"
#include "CGAL/Triangulation_cell_base_with_circumcenter_3.h"
#include <CGAL/HalfedgeDS_vector.h>
#include "CGAL/Polyhedron_3.h"
#include "CGAL/Simple_cartesian.h"
#include "CGAL/Bbox_3.h"

#include <CGAL/Surface_mesh_default_triangulation_3.h>
#include <CGAL/Surface_mesh_default_criteria_3.h>
#include <CGAL/Complex_2_in_triangulation_3.h>
#include <CGAL/IO/Complex_2_in_triangulation_3_file_writer.h>
#include <CGAL/make_surface_mesh.h>
#include <CGAL/Implicit_surface_3.h>

#include "VDB.h"
#include "ofMesh.h"

typedef CGAL::Exact_predicates_inexact_constructions_kernel K;

typedef CGAL::Triangulation_vertex_base_with_info_3<int, K, CGAL::Triangulation_vertex_base_3<K> > Vb2;
typedef CGAL::Surface_mesh_vertex_base_3<K, Vb2> Smvb;
typedef CGAL::Surface_mesh_cell_base_3<K> Smcb;
typedef CGAL::Triangulation_data_structure_3<Smvb, Smcb > Tds2;
typedef CGAL::Delaunay_triangulation_3<K, Tds2> Tr;

// c2t3
typedef CGAL::Complex_2_in_triangulation_3<Tr> C2t3;

typedef Tr::Geom_traits GT;

typedef Tr::Finite_facets_iterator Finite_facets_iterator;
typedef Tr::Finite_vertices_iterator Finite_vertices_iterator;
typedef Tr::Facet Facet;
typedef Tr::Edge Edge;
typedef Tr::Vertex_handle Vertex_handle;

struct VDBImplicit {
	VDB vdb;

	const GT::FT operator()(GT::Point_3 pt) const {
		return vdb.samplePt(ofVec3f(pt.x(), pt.y(), pt.z()));
	}
};

typedef CGAL::Implicit_surface_3<GT, VDBImplicit> Surface_3;


static void complexToMesh(C2t3 &c2t3, ofMesh & mesh) {
	mesh.clear();

	using CGAL::Surface_mesher::number_of_facets_on_surface;


	const Tr& tr = c2t3.triangulation();

	// Finite vertices coordinates.
	std::map<Vertex_handle, int> V;
	int inum = 0;
	int jnum = 0;
	for (Finite_vertices_iterator vit = tr.finite_vertices_begin();
	vit != tr.finite_vertices_end();
		++vit)
	{
		V[vit] = inum++;
		mesh.addVertex(ofVec3f(vit->point().x(), vit->point().y(), vit->point().z()));
	}

	//facets
	jnum = 0;
	Finite_facets_iterator fit = tr.finite_facets_begin();
	std::set<Facet> oriented_set;
	std::stack<Facet> stack;

	int number_of_facets = c2t3.number_of_facets();

	while (oriented_set.size() != number_of_facets)
	{
		while (fit->first->is_facet_on_surface(fit->second) == false ||
			oriented_set.find(*fit) != oriented_set.end() ||

			oriented_set.find(c2t3.opposite_facet(*fit)) !=
			oriented_set.end())
		{
			++fit;
		}
		oriented_set.insert(*fit);
		stack.push(*fit);
		while (!stack.empty())
		{
			Facet f = stack.top();
			stack.pop();
			for (int ih = 0; ih < 3; ++ih) {
				const int i1 = tr.vertex_triple_index(f.second, tr.cw(ih));
				const int i2 = tr.vertex_triple_index(f.second, tr.ccw(ih));

				const C2t3::Face_status face_status
					= c2t3.face_status(Edge(f.first, i1, i2));
				if (face_status == C2t3::REGULAR) {
					Facet fn = c2t3.neighbor(f, ih);
					if (oriented_set.find(fn) == oriented_set.end()) {
						if (oriented_set.find(c2t3.opposite_facet(fn)) == oriented_set.end())
						{
							oriented_set.insert(fn);
							stack.push(fn);
						}
					}
				}
			} // end "for each neighbor of f"
		} // end "stack non empty"
	} // end "oriented_set not full"

	cout << "get faces" << endl;

	for (std::set<Facet>::const_iterator fit =
		oriented_set.begin();
		fit != oriented_set.end();
		++fit)
	{
		const Tr::Cell_handle cell = fit->first;
		const int& index = fit->second;
		const int index1 = V[cell->vertex(tr.vertex_triple_index(index, 0))];
		const int index2 = V[cell->vertex(tr.vertex_triple_index(index, 1))];
		const int index3 = V[cell->vertex(tr.vertex_triple_index(index, 2))];
		//off.write_facet(index1, index2, index3);
		mesh.addIndex(index1);
		mesh.addIndex(index2);
		mesh.addIndex(index3);
	}

	cout << "done" << endl;
}


static void buildMesh(VDB & voxel, ofMesh & mesh, float maxTriangle = .45, float maxError = .1) {
	Tr tr;            // 3D-Delaunay triangulation

					  //use exact method with divided space;

	C2t3 c2t3(tr);   // 2D-complex in 3D-Delaunay triangulation

					 //get inside pt
	ofVec3f insidePt;
	typedef openvdb::FloatGrid GridType;
	typedef GridType::TreeType TreeType;
	// Iterate over references to const LeafNodes.
	for (TreeType::LeafCIter iter = voxel.grid->tree().cbeginLeaf(); iter; ++iter) {
		const TreeType::LeafNodeType& leaf = *iter;
		if (leaf.getFirstValue() < 0) {
			openvdb::CoordBBox box = leaf.getNodeBoundingBox();
			openvdb::Vec3d pt = voxel.grid->indexToWorld(box.getCenter());
			insidePt.x = pt.x();
			insidePt.y = pt.y();
			insidePt.z = pt.z();
		}
	}

	

	GT::Point_3 bounding_sphere_center(insidePt.x, insidePt.y, insidePt.z);
	auto bbox = voxel.bbox();
	GT::FT bounding_sphere_squared_radius = (bbox.second - bbox.first).lengthSquared();
	GT::Sphere_3 bounding_sphere(bounding_sphere_center,
		bounding_sphere_squared_radius);

	VDBImplicit srfFunc;
	srfFunc.vdb = voxel;
	Surface_3 surface(srfFunc, bounding_sphere, 1e-4);

	// defining meshing criteria
	//CGAL::Surface_mesh_default_criteria_3<Tr> criteria(20.,.5,.5); //mesh detail for lamps
	// CGAL::Surface_mesh_default_criteria_3<Tr> criteria(20.,.25,.25); //mesh detail for jewelry
	CGAL::Surface_mesh_default_criteria_3<Tr> criteria(20., maxTriangle, maxError);

	// meshing surface, with the "non-manifold" algorithm
	//CGAL::make_surface_mesh(c2t3, surface, criteria, CGAL::Non_manifold_tag());
	CGAL::make_surface_mesh(c2t3, surface, criteria, CGAL::Manifold_tag());

	cout << "convert to mesh" << endl;
	complexToMesh(c2t3, mesh);
}