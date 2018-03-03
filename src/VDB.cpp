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
	filter.normalize();
	filter.prune();
	isUpdated = false;
}

void VDB::offset(float amt, VDB & mask) {
	LevelSetFilter<FloatGrid> filter(*grid);
	filter.offset(amt, &*(mask.grid));
	filter.normalize();
	filter.prune();
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
	csgDifference(*grid, *cGrid);
	isUpdated = false;
}

void VDB::doUnion(FloatGrid::Ptr vdb) {
	csgUnion(*grid, *(vdb));
	isUpdated = false;
}

void VDB::doIntersect(VDB & vdb) {
	//copy
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
	csgIntersection(*grid, *cGrid);
	isUpdated = false;
}

void VDB::blur() {
	LevelSetFilter<FloatGrid> filter(*grid);
	//filter.gaussian();
	filter.laplacian();
	isUpdated = false;
}

void VDB::smooth() {
	LevelSetFilter<FloatGrid> filter(*grid);
	//filter.gaussian();
	filter.meanCurvature();
	isUpdated = false;
}

void VDB::taubin() {
	LevelSetFilter<FloatGrid> filter(*grid);
	//filter.gaussian();
	//filter.taubin();
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
	openvdb::tools::VolumeToMesh mesher(0);
	mesher(*grid);

	mesh.clear();

	openvdb::Coord ijk;

	for (Index64 n = 0, i = 0, N = mesher.pointListSize(); n < N; ++n) {
		const openvdb::Vec3s& p = mesher.pointList()[n];
		mesh.addVertex(ofVec3f(p[0], p[1], p[2]));
		mesh.addNormal(ofVec3f(0,0,0));
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

			e1 = mesh.getVertex(quad[2]);
			e1 -= mesh.getVertex(quad[0]);
			e2 = mesh.getVertex(quad[3]);
			e2 -= mesh.getVertex(quad[1]);
			norm = e1.cross(e2);

			//const float length = norm.length();
			//if (length > 1.0e-6) norm /= length;
			//norm *= -1;
			for (int v = 0; v < 4; ++v) {
				mesh.setNormal(quad[v], mesh.getNormal(quad[v]) + norm);
				//mesh.setNormal(quad[v], norm);
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
	grid = FloatGrid::create(1.2f);
	grid->setGridClass(GRID_LEVEL_SET);
	mesh.enableNormals();
	isovalue = 0.0;
	isUpdated = false;
}

VDB::VDB(ofMesh & m, float resolution) {
	grid = FloatGrid::create(resolution*3);
	mesh.enableNormals();
	isUpdated = false;
	isovalue = 0.0;
	loadMesh(m, resolution);
}

VDB::VDB(const VDB & vdb) {
	grid = vdb.grid->deepCopy();
	mesh.enableNormals();
	isovalue = vdb.isovalue;
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
	LevelSetRayIntersector<FloatGrid> intersector(*grid, isovalue);
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
	//if (grid->getGridClass() == GRID_LEVEL_SET) {
		LevelSetRayIntersector<FloatGrid> intersector(*grid, 0);
		Vec3R pt;
		Real t0;
		bool val = intersector.intersectsWS(math::Ray<Real>(Vec3R(x, y, z), Vec3R(dx, dy, dz)), pt, t0);
		t = t0;
		ox = pt[0];
		oy = pt[1];
		oz = pt[2];
		return val;
		/*
		}

	else {
		VolumeRayIntersector<FloatGrid> intersector(*grid, isovalue);
		Vec3R pt;
		intersector.setWorldRay(math::Ray<Real>(Vec3R(x, y, z), Vec3R(dx, dy, dz)));
		Real t0 = 0, t1 = 0;
		while (intersector.march(t0,t1)) {
			t = t0;
			pt = intersector.getWorldPos(t0);
			ox = pt[0];
			oy = pt[1];
			oz = pt[2];
			return true;
		}
		return false;
	}*/
}

bool VDB::intersectRay(const ofVec3f & pt, const ofVec3f & dir, ofVec3f & out, float &t) {
	return intersectRay(pt.x, pt.y, pt.z, dir.x, dir.y, dir.z, out.x, out.y, out.z, t);
}

struct MatAdd {
	float f;
	MatAdd(float _f) : f(_f) {}
	inline void operator()(const FloatGrid::ValueOnIter& iter) const {
		iter.setValue(*iter + f);
	}
};


void VDB::setThreshold(float thresh) {

	foreach(grid->beginValueOn(), MatAdd(thresh-isovalue));
	isovalue = thresh;
	isUpdated = false;
}


void VDB::loadVol(ifstream & buf, int w, int h, int d, float resolution) {
	grid = FloatGrid::create(100);
	FloatGrid::Accessor acc = grid->getAccessor();
	Coord ijk;
	grid->setGridClass(GRID_LEVEL_SET);
	
	float threshold = .05;
	float f;
	float minV = 9e9;
	float maxV = -9e9;
	int &x = ijk[0], &y = ijk[1], &z = ijk[2];
	for (x = 0; x < w; ++x) {
		for (y = 0; y < h; ++y) {
			for (z = 0; z < d; ++z) {
				buf.read((char *)&f, sizeof(f));
				if (f < threshold) {
					f = 0;
				}
				else {
					f = -f;
				}
				minV = min(minV, f);
				maxV = max(maxV, f);
				acc.setValue(ijk, f);
			}
		}
	}
	cout << "Loaded volume. min value: " << minV << " max value: " << maxV << endl;
	//grid-> = max(abs(maxV), abs(minV));
	grid->pruneGrid();
	setThreshold(-.5);
	math::Transform trans;
	trans.preScale(resolution);
	grid->transform() = trans;
	isUpdated = false;
}

float VDB::samplePt(const ofVec3f & pt) const {
	GridSampler<FloatGrid, BoxSampler> sampler(*grid);
	return sampler.wsSample(Vec3R(pt.x, pt.y, pt.z));
}

void VDB::rasterSphere(ofVec3f center, float rad, float val) {
	FloatGrid::Accessor acc = grid->getAccessor();
	//const SdfT dx = SdfT(mParent.mDx), w = SdfT(mParent.mHalfWidth);
	const float max = rad;// maximum distance in voxel units
	const Coord a(math::Floor(center.x - max), math::Floor(center.y - max), math::Floor(center.z - max));
	const Coord b(math::Ceil(center.x + max), math::Ceil(center.y + max), math::Ceil(center.z + max));
	const float max2 = math::Pow2(max);//square of maximum distance in voxel units
	//const float min2 = math::Pow2(math::Max(0, rad));//square of minimum distance
	float v;
	for (Coord c = a; c.x() <= b.x(); ++c.x()) {
		float x2 = math::Pow2(c.x() - center.x);
		for (c.y() = a.y(); c.y() <= b.y(); ++c.y()) {
			float x2y2 = x2 + math::Pow2(c.y() - center.y);
			for (c.z() = a.z(); c.z() <= b.z(); ++c.z()) {
				float x2y2z2 = x2y2 + math::Pow2(c.z() - center.z);//square distance from c to P
				if (x2y2z2 >= max2)
					continue;//outside narrow band of the particle or inside existing level set
				//if (x2y2z2 <= min2) {//inside narrow band of the particle.
				//	acc.setValueOff(c, inside);
				//	continue;
				//}
				// convert signed distance from voxel units to world units
				//const ValueT d=dx*(math::Sqrt(x2y2z2) - R);
				acc.setValue(c, val);
			}//end loop over z
		}//end loop over y
	}//end loop over x
}