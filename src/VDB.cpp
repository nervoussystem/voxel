#include "VDB.h"
#include <openvdb/tools/RayIntersector.h>
#include <openvdb/tools/VolumeToMesh.h>
#include <openvdb/tools/MeshToVolume.h>
#include <openvdb/tools/LevelSetFilter.h>
#include <openvdb/tools/Composite.h>
#include <openvdb/tools/LevelSetUtil.h>
#include <openvdb/tools/GridTransformer.h>

#include "ofMatrix4x4.h"
#include "ofImage.h"
#include "ofGraphics.h"

using namespace openvdb;
using namespace openvdb::tools;

struct MeshWrapper {
	ofMesh * mesh;
	size_t polygonCount() const {
		return mesh->getNumIndices() / 3;
	}
	size_t pointCount() {
		return mesh->getNumVertices();
	}
	size_t vertexCount(size_t n) const { // Vertex count for polygon n
		return 3;
	}

	// Return position pos in local grid index space for polygon n and vertex v
	void getIndexSpacePoint(size_t n, size_t v, openvdb::Vec3d& pos) const {
		ofVec3f pt = mesh->getVertex(mesh->getIndex(n * 3 + v));
		pos[0] = pt.x;
		pos[1] = pt.y;
		pos[2] = pt.z;
	}

	MeshWrapper(ofMesh * m) : mesh(m) {}
};

void VDB::loadMesh(ofMesh & toLoad, float resolution, int band) {
	MeshWrapper mesher(&toLoad);
	for (auto & v : mesher.mesh->getVertices()) {
		v /= resolution;
	}
	math::Transform trans;
	trans.preScale(resolution);

	grid = meshToVolume<openvdb::FloatGrid>(mesher, trans, band, band);

	for (auto & v : mesher.mesh->getVertices()) {
		v *= resolution;
	}

	isUpdated = false;
}

void VDB::offset(float amt) {
	LevelSetFilter<FloatGrid> filter(*grid);
	filter.offset(amt);

	isUpdated = false;
}

void VDB::load(string filename) {
}

void VDB::doUnion(VDB & vdb) {
	const math::Transform
		&sourceXform = vdb.grid->transform(),
		&targetXform = grid->transform();
	FloatGrid::Ptr cGrid = createLevelSet<FloatGrid>(grid->voxelSize()[0]);
	cGrid->transform() = grid->transform();
	// Compute a source grid to target grid transform.
	openvdb::Mat4R xform =
		sourceXform.baseMap()->getAffineMap()->getMat4() *
		targetXform.baseMap()->getAffineMap()->getMat4().inverse();
	// Create the transformer.
	openvdb::tools::GridTransformer transformer(xform);
	
	// Resample using trilinear interpolation.
	transformer.transformGrid<openvdb::tools::BoxSampler, openvdb::FloatGrid>(
		*vdb.grid, *cGrid);

	csgUnion(*grid, *cGrid);
	isUpdated = false;
}

void VDB::doDifference(VDB & vdb) {
	//copy
	FloatGrid::Ptr cGrid = FloatGrid::create();
	cGrid = vdb.grid->deepCopy();
	csgDifference(*grid, *cGrid);
	isUpdated = false;
}

void VDB::doUnion(FloatGrid::Ptr vdb) {
	csgUnion(*grid, *(vdb));
	isUpdated = false;
}

void VDB::doIntersect(VDB & vdb) {
	//copy
	FloatGrid::Ptr cGrid = FloatGrid::create();
	cGrid = vdb.grid->deepCopy(); 
	csgIntersection(*grid, *(vdb.grid));
	isUpdated = false;
}

void VDB::blur() {
	LevelSetFilter<FloatGrid> filter(*grid);
	//filter.gaussian();
	filter.laplacian();
	isUpdated = false;
}

void VDB::clear() {
	grid->clear();
}

ofMesh VDB::toMesh() {
	if (!isUpdated) {
		updateMesh();
	}
	return mesh;
}

void VDB::floodFill() {
	signedFloodFill(grid->tree());
}

void VDB::updateMesh() {
	
	//openvdb::tools::VolumeToMesh mesher(grid->getGridClass() == openvdb::GRID_LEVEL_SET ? 0.0 : 0.01);
	openvdb::tools::VolumeToMesh mesher(0.0);
	mesher(*grid);

	mesh.clear();

	openvdb::Coord ijk;

	for (Index64 n = 0, i = 0, N = mesher.pointListSize(); n < N; ++n) {
		const openvdb::Vec3s& p = mesher.pointList()[n];
		mesh.addVertex(ofVec3f(p[0], p[1], p[2]));
		mesh.addNormal(ofVec3f(0,0,1));
	}

	// Copy primitives
	openvdb::tools::PolygonPoolList& polygonPoolList = mesher.polygonPoolList();
	Index64 numQuads = 0;
	for (Index64 n = 0, N = mesher.polygonPoolListSize(); n < N; ++n) {
		numQuads += polygonPoolList[n].numQuads();
	}

	ofVec3f norm, e1, e2;

	for (Index64 n = 0, N = mesher.polygonPoolListSize(); n < N; ++n) {
		const openvdb::tools::PolygonPool& polygons = polygonPoolList[n];
		for (Index64 i = 0, I = polygons.numQuads(); i < I; ++i) {
			const openvdb::Vec4I& quad = polygons.quad(i);
			mesh.addIndex(quad[0]);
			mesh.addIndex(quad[2]);
			mesh.addIndex(quad[1]);
			mesh.addIndex(quad[0]);
			mesh.addIndex(quad[3]);
			mesh.addIndex(quad[2]);

			e1 = mesh.getVertex(quad[1]);
			e1 -= mesh.getVertex(quad[0]);
			e2 = mesh.getVertex(quad[2]);
			e2 -= mesh.getVertex(quad[1]);
			norm = e1.cross(e2);

			const float length = norm.length();
			if (length > 1.0e-6) norm /= length;
			//norm *= -1;
			for (int v = 0; v < 4; ++v) {
				mesh.setNormal(quad[v], norm);
			}
		}
	}
	tempTransform.makeIdentityMatrix();
	isUpdated = true;
}

void VDB::transform(ofMatrix4x4 & mat) {
	tempTransform.postMult(mat);
	math::Mat4d vMat(mat.getPtr());
	grid->transform().postMult(vMat);
}

void VDB::translate(ofVec3f dir) {
	tempTransform.translate(dir);
	grid->transform().postTranslate(Vec3d(dir.x, dir.y, dir.z));
}

void VDB::rotate(const ofVec3f & axis, float angle) {
	Mat4d mat;
	mat.setToRotation(Vec3R(axis.x, axis.y, axis.z), angle);
	tempTransform.rotateRad(angle, axis.x, axis.y, axis.z);
	grid->transform().postMult(mat);
}

void VDB::draw() {
	if (!isUpdated) updateMesh();
	ofPushMatrix();
	ofMultMatrix(tempTransform);
	mesh.draw();
	ofPopMatrix();
}

void VDB::save(string filename) {

}

VDB::VDB() {
	grid = FloatGrid::create(1.2);
	grid->setGridClass(GRID_LEVEL_SET);
	mesh.enableNormals();
	isUpdated = false;
}

VDB::VDB(ofMesh & m, float resolution) {
	grid = FloatGrid::create(resolution*3);
	mesh.enableNormals();
	isUpdated = false;
	loadMesh(m, resolution);
}

VDB::VDB(const VDB & vdb) {
	grid = vdb.grid->deepCopy();
	mesh.enableNormals();
	isUpdated = false;
}

void VDB::toEmber(string filename) {
	//transform grid
	FloatGrid::Ptr xGrid = FloatGrid::create();

	const openvdb::math::Transform
		&sourceXform = grid->transform();

	math::Transform targetXform;
	//TODO: deal with layer height
	targetXform.preScale(.05);

	Mat4R xform =
		sourceXform.baseMap()->getAffineMap()->getMat4() *
		targetXform.baseMap()->getAffineMap()->getMat4().inverse();

	GridTransformer transformer(xform);
	transformer.transformGrid<BoxSampler, FloatGrid>(*grid, *xGrid);

	//make bool
	BoolGrid::ConstPtr bGrid = sdfInteriorMask(*xGrid);

	openvdb::io::File file(filename + ".vdb");
	// Add the grid pointer to a container.
	openvdb::GridCPtrVec grids;
	grids.push_back(bGrid);
	// Write out the contents of the container.
	file.write(grids);
	file.close();

	CoordBBox bbox = bGrid->evalActiveVoxelBoundingBox();
	
	cout << bbox.getStart() << " " << bbox.getEnd() << endl;

	if (bbox.dim().x() > 1280 || bbox.dim().y() > 800) {
		cout << "TOO BIG FOR EMBER " << bbox.dim() << endl;
		return;
	}

	//TODO: CENTER
	int zStart = bbox.min().z(),
		zEnd = bbox.max().z(),
		xStart = bbox.min().x(),
		xEnd = bbox.max().x(),
		yStart = bbox.min().y(),
		yEnd = bbox.max().y();

	int xOff = -xStart;
	int yOff = -yStart;
	ofImage img;

	const int imgW = xEnd - xStart+1;
	const int imgH = yEnd - yStart+1;
	unsigned char * pval = new unsigned char[imgW * imgH];
	//img.setImageType(OF_IMAGE_GRAYSCALE);
	//img.allocate(1280, 800, OF_IMAGE_GRAYSCALE);

	BoolGrid::ConstAccessor acc = bGrid->getAccessor();
	openvdb::Coord ijk;
	int slice = 1;
	const ofColor white = ofColor::white;
	int &x = ijk[0], &y = ijk[1], &z = ijk[2];
	for (z = zStart; z <= zEnd; ++z) {
		for (int i = 0; i < imgW * imgH; ++i) pval[i] = 0;
		//pval = img.getPixels().getData();
		bool found = false;
		for (x = xStart; x <= xEnd; ++x) {
			for (y = yStart; y <= yEnd; ++y) {
				if (acc.getValue(ijk)) {
					found = true;
					pval[(y + yOff)*imgW + x + xOff] =  255;
				}
			}
		}
		if (found || slice != 1) {
			img.setFromPixels(pval, imgW, imgH, OF_IMAGE_GRAYSCALE);
			img.save(filename + "/slice_" + to_string(slice) + ".png");
			cout << slice << endl;
			slice++;
		}
	}

	delete[] pval;

}

pair<ofVec3f, ofVec3f> VDB::bbox() {
	math::CoordBBox bbox = grid->evalActiveVoxelBoundingBox();
	Coord minC = bbox.getStart();
	Coord maxC = bbox.getEnd();
	Vec3d minPt = grid->indexToWorld(minC);
	Vec3d maxPt = grid->indexToWorld(maxC);
	return make_pair(ofVec3f(minPt.x(), minPt.y(), minPt.z()), ofVec3f(maxPt.x(), maxPt.y(), maxPt.z()));
}

bool VDB::intersectRay(const float x, const float y, const float z, const float dx, const float dy, const float dz, float & ox, float &oy, float &oz) {
	LevelSetRayIntersector<FloatGrid> intersector(*grid);
	Vec3R pt;
	
	bool val = intersector.intersectsWS(math::Ray<Real>(Vec3R(x, y, z), Vec3R(dx, dy, dz)), pt);
	ox = pt[0];
	oy = pt[1];
	oz = pt[2];
	return val;
}

bool VDB::intersectRay(const ofVec3f & pt, const ofVec3f & dir, ofVec3f & out) {
	return intersectRay(pt.x, pt.y, pt.z, dir.x, dir.y, dir.z, out.x, out.y, out.z);
}

bool VDB::intersectRay(const float x, const float y, const float z, const float dx, const float dy, const float dz, float & ox, float &oy, float &oz, float & t) {
	LevelSetRayIntersector<FloatGrid> intersector(*grid);
	Vec3R pt;
	Real t0;
	bool val = intersector.intersectsWS(math::Ray<Real>(Vec3R(x, y, z), Vec3R(dx, dy, dz)), pt, t0);
	t = t0;
	ox = pt[0];
	oy = pt[1];
	oz = pt[2];
	return val;
}

bool VDB::intersectRay(const ofVec3f & pt, const ofVec3f & dir, ofVec3f & out, float &t) {
	return intersectRay(pt.x, pt.y, pt.z, dir.x, dir.y, dir.z, out.x, out.y, out.z, t);
}
