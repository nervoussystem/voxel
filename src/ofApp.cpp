#include "ofApp.h"
#include "LevelSetCapsule.h"

//--------------------------------------------------------------
void ofApp::setup(){
	grid.grid = openvdb::tools::createLevelSetCapsule<openvdb::FloatGrid>(3, openvdb::Vec3f(0, 0, 0), openvdb::Vec3f(20, -20, 20), 0.5);
	grid.isUpdated = false;
	resolution = .3;
}

int step = 200;
string folder = "C:\\Users\\nervous system\\of_v0.8.0_vs_release\\apps\\myApps\\3DVoronoi\\bin\\data\\opt2\\lines\\SLAB_sz9.5_b10_";
//--------------------------------------------------------------
void ofApp::update(){
	if (step < 150) {
		int frame = (int) (0.5*.86*pow(step,1.5)+1);
		loadLines(folder + to_string(frame) + ".csv");
		grid.updateMesh();
		grid.mesh.save(folder + to_string(step) + ".ply");
		step++;
	}
}

//--------------------------------------------------------------
void ofApp::draw(){
	ofBackground(255);
	cam.begin();
	ofEnableDepthTest();
	ofEnableLighting();
	ofLight light0;
	light0.enable();
	ofSetColor(100);
	grid.draw();
	ofVec2f mousePt(ofGetMouseX(), ofGetMouseY());
	ofVec3f sPt = cam.screenToWorld(ofVec3f(mousePt.x,mousePt.y, -.9)) ;
	ofVec3f sPt2 = cam.screenToWorld(ofVec3f(mousePt.x, mousePt.y,.9)) ;
	ofVec3f srfPt;
	if (grid.intersectRay(sPt, (sPt2 - sPt), srfPt)) {
		cout << sPt << " : " << sPt2 << " : " << srfPt << endl;
		ofPushMatrix();
		ofTranslate(srfPt);
		ofSphere(2);
		ofPopMatrix();
	}
	grid.isUpdated = false;
	cam.end();
}

void ofApp::loadLines(string filename) {
	ifstream in(filename);
	string cLine, token;
	float f;
	openvdb::Vec3f v1, v2;
	float cthick;
	grid.clear();
	int count = 0;
	while (getline(in, cLine)) {
		stringstream ss(cLine);
		int num = 0;
		if (cLine.length() > 12) {
			while (getline(ss, token, ',')) {
				f = ::atof(token.c_str());
				switch (num) {
				case 0:
					v1[0] = f;
					break;
				case 1:
					v1[1] = f;
					break;
				case 2:
					v1[2] = f;
					break;
				case 3:
					v2[0] = f;
					break;
				case 4:
					v2[1] = f;
					break;
				case 5:
					v2[2] = f;
					break;
				case 6:
					cthick = f;
					break;
				}
				num++;
			}
			//ofVec3f dir = v2-v1;
			//float len = dir.length();
			//dir /= len;
			//float offset = min(cthick*0.25, len/3.0);
			//check total len
			//v1 += dir*offset;
			//v2 -= dir*offset;

			//openvdb::FloatGrid::Ptr mgrid = openvdb::tools::createLevelSetCapsule<openvdb::FloatGrid>(cthick, v1,v2, 0.5);
			//grid.doUnion(mgrid);
			openvdb::tools::LevelSetCapsule<openvdb::FloatGrid> factory(cthick, v1, v2);
			factory.mGrid = grid.grid;
			factory.rasterCapsule(0.5, 3);
			cout << count++ << endl;
		}
	}
	grid.isUpdated = false;

	in.close();
}

void ofApp::dragEvent(ofDragInfo info) {
	VDB subBox;
	ofMesh box = ofMesh::box(60, 320, 45);
	for (auto & v : box.getVertices()) {
		v -= ofVec3f(60-30, 7-160, 9-22.5);
	}
	subBox.loadMesh(box, 0.5);
	if (info.files.size() > 0) {
		for (int i = 0; i < info.files.size();i++) {
			if (info.files[i].substr(info.files[i].size() - 3) == "ply") {
				ofMesh mesh;
				cout << info.files[i] << endl;
				mesh.load(info.files[i]);
				if (mesh.getNumIndices() > 3) {
					cout << "vertices " << mesh.getNumVertices() << endl;
					if (ofGetKeyPressed(OF_KEY_CONTROL)) {
						VDB temp;
						temp.loadMesh(mesh, resolution);
						grid.doDifference(temp);
					}
					else {
						grid.loadMesh(mesh, resolution);
						cout << "loaded mesh to voxels" << endl;
					}
				}
			}
			else if (info.files[i].substr(info.files[i].size() - 3) == "csv") {
				loadLines(info.files[i]);
				//grid.doDifference(subBox);
				grid.updateMesh();
				//grid.mesh.save(info.files[i].substr(0,info.files[i].size() - 3) + "ply");
			}
		}
	}
}

//--------------------------------------------------------------
void ofApp::keyPressed(int key){
	if (key == 's') {
		grid.mesh.save(ofToDataPath("test.ply"));
	}
	else if (key == 'e') {
		grid.toEmber("cats");
	}
	else if (key == 'b') {
		grid.blur();
	}
	else if (key == 'o') {
		grid.offset(-.05);
	}
}

//--------------------------------------------------------------
void ofApp::keyReleased(int key){
}

//--------------------------------------------------------------
void ofApp::mouseMoved(int x, int y ){

}

//--------------------------------------------------------------
void ofApp::mouseDragged(int x, int y, int button){

}

//--------------------------------------------------------------
void ofApp::mousePressed(int x, int y, int button){

}

//--------------------------------------------------------------
void ofApp::mouseReleased(int x, int y, int button){

}

//--------------------------------------------------------------
void ofApp::mouseEntered(int x, int y){

}

//--------------------------------------------------------------
void ofApp::mouseExited(int x, int y){

}

//--------------------------------------------------------------
void ofApp::windowResized(int w, int h){

}

//--------------------------------------------------------------
void ofApp::gotMessage(ofMessage msg){

}
