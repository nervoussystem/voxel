#pragma once
#include "ofVbo.h"
class ObjMesh {
public:
	vector<GLuint> triangleIndices, lineIndices, sphereIndices;
	vector<float> thickness;
	vector<ofVec3f> positions;
	bool hasThickness = false;

	template<class ContainerT>
	void tokenize(const std::string& str, ContainerT& tokens,
		const std::string& delimiters = " ", bool trimEmpty = false);
	void clear();
	void load(string path);
};