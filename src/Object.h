#pragma once
#include "ofVboMesh.h"

class MeshObject{
public:
	ofVboMesh mesh;
	vector<float> thicknesses;

	void draw() {
		mesh.draw();
	}
};