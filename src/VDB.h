#pragma once

#include <openvdb/openvdb.h>
#include "ofVboMesh.h"


class VDB {
public:
	typedef shared_ptr<VDB> Ptr;
	openvdb::FloatGrid::Ptr grid;
	ofVboMesh mesh;
	bool isUpdated;

	void loadMesh(ofMesh & toLoad, float resolution = 1.0, int band = 3);
	void clear();
	void offset(float amt);
	void load(string filename);
	void doUnion(VDB & vdb);
	void doDifference(VDB & vdb);
	void doUnion(openvdb::FloatGrid::Ptr vdb);
	void doIntersect(VDB & vdb);

	bool intersectRay(const float x, const float y, const float z, const float dx, const float dy, const float dz, float & ox, float &oy, float &oz);
	bool intersectRay(const ofVec3f & pt, const ofVec3f & dir, ofVec3f & out);
	void blur();
	void updateMesh();

	void draw();
	void save(string filename);
	void toEmber(string filename);

	ofMesh toMesh();

	VDB();
	VDB(const VDB & vdb);
};