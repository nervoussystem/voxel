#pragma once

#include <openvdb/openvdb.h>

#include <openvdb/tools/GridTransformer.h>
#include "ofVboMesh.h"
#include "ofMatrix4x4.h"

class VDB {
public:
	ofMatrix4x4 tempTransform;
	typedef shared_ptr<VDB> Ptr;
	openvdb::FloatGrid::Ptr grid;
	ofVboMesh mesh;
	float isovalue;
	bool isUpdated;

	void loadMesh(ofMesh & toLoad, float resolution = 1.0, int band = 3);
	void loadVol(ifstream & buf, int w, int h , int d, float resolution = 1.0);
	void clear();
	void offset(float amt);
	void offset(float amt, VDB & mask);
	void load(string filename);
	void doUnion(VDB & vdb);
	void doDifference(VDB & vdb);
	void doUnion(openvdb::FloatGrid::Ptr vdb);
	void doIntersect(VDB & vdb);
	void setThreshold(float thresh);

	void floodFill();

	bool intersectRay(const float x, const float y, const float z, const float dx, const float dy, const float dz, float & ox, float &oy, float &oz);
	bool intersectRay(const ofVec3f & pt, const ofVec3f & dir, ofVec3f & out);
	bool intersectRay(const float x, const float y, const float z, const float dx, const float dy, const float dz, float & ox, float &oy, float &oz, float &t);
	bool intersectRay(const ofVec3f & pt, const ofVec3f & dir, ofVec3f & out, float &t);
	void blur();
	void smooth();
	void taubin();

	void updateMesh();

	void transform(ofMatrix4x4 & mat);
	void translate(ofVec3f dir);
	void rotate(const ofVec3f & axis, float angle);
	pair<ofVec3f, ofVec3f> bbox();
	void draw();
	void save(string filename);
	void toEmber(string filename);
	float samplePt(const ofVec3f & pt) const;

	ofMesh toMesh();

	VDB();
	VDB(ofMesh & m, float resolution = 1.0);
	VDB(const VDB & vdb);

	void rasterSphere(ofVec3f center, float rad, float val);
};