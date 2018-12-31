#pragma once
#include "ObjMesh.h"

template<class ContainerT>
void ObjMesh::tokenize(const std::string& str, ContainerT& tokens,const std::string& delimiters, bool trimEmpty)
{
	std::string::size_type pos, lastPos = 0, length = str.length();

	using value_type = typename ContainerT::value_type;
	using size_type = typename ContainerT::size_type;

	while (lastPos < length + 1)
	{
		pos = str.find_first_of(delimiters, lastPos);
		if (pos == std::string::npos)
		{
			pos = length;
		}

		if (pos != lastPos || !trimEmpty)
			tokens.push_back(value_type(str.data() + lastPos,
			(size_type)pos - lastPos));

		lastPos = pos + 1;
	}
}


void ObjMesh::clear() {
	triangleIndices.clear();
	lineIndices.clear();
	sphereIndices.clear();
	thickness.clear();
	positions.clear();
}

void ObjMesh::load(string path) {
	clear();
	ifstream in(path);
	hasThickness = false;
	string line;
	list<string> tokens;
	while (getline(in, line)) {
		tokens.clear();
		tokenize(line, tokens);
		if (tokens.size() > 1) {
			if (tokens.front() == "v") {
				if (tokens.size() > 3) {
					auto it = tokens.begin();
					it++;
					ofVec3f pos;
					pos.x = stof(*it);
					it++;
					pos.y = stof(*it);
					it++;
					pos.z = stof(*it);
					it++;
					if (it != tokens.end()) {
						float t = stof(*it);
						thickness.push_back(t);
						hasThickness = true;
					}
					else {
						thickness.push_back(-1);
					}
					positions.push_back(pos);
				}
			}
			else if (tokens.front() == "l") {
				auto it = tokens.begin();
				it++;
				lineIndices.push_back(stoi(*it) - 1);
				it++;
				lineIndices.push_back(stoi(*it) - 1);
			}
			else if (tokens.front() == "f") {
				auto it = tokens.begin();
				it++;
				int firstIndex = stoi(*it) - 1;
				it++;
				int prevIndex = stoi(*it) - 1;
				it++;
				while (it != tokens.end()) {
					int cIndex = stoi(*it) - 1;
					triangleIndices.push_back(firstIndex);
					triangleIndices.push_back(prevIndex);
					triangleIndices.push_back(cIndex);
					prevIndex = cIndex;
					it++;
				}
			}
			else if (tokens.front() == "p") {
				auto it = tokens.begin();
				it++;
				sphereIndices.push_back(stoi(*it) - 1);
			}
		}
	}

}