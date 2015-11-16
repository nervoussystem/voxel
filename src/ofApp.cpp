#include "ofApp.h"
#include "LevelSetCapsule.h"

//--------------------------------------------------------------
void ofApp::setup(){
	openvdb::initialize();
	gumball.setCamera(cam);
	ofAddListener(gumball.gumballEvent, this, &ofApp::gumballEvent);
	resolution = .3;
	setupGui();
}

int step = 200;
string folder = "C:\\Users\\nervous system\\of_v0.8.0_vs_release\\apps\\myApps\\3DVoronoi\\bin\\data\\opt2\\lines\\SLAB_sz9.5_b10_";
//--------------------------------------------------------------
void ofApp::update() {
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
	for (auto g : grids) {
		if (isSelected(g)) {
			ofSetColor(200, 0, 0);
		}
		else {
			ofSetColor(120);
		}
		g->draw();
	}
	ofDisableDepthTest();
	ofDisableLighting();
	gumball.draw();

	cam.end();

	gui.draw();
}

bool ofApp::isSelected(VDB::Ptr g) {
	return find(selected.begin(), selected.end(), g) != selected.end();
}

void ofApp::doDelete() {
	for (auto sel : selected) {
		auto it = find(grids.begin(), grids.end(), sel);
		if (it != grids.end()) grids.erase(it);
	}
	selected.clear();
}

void ofApp::loadLines(string filename) {
	ifstream in(filename);
	string cLine, token;
	float f;
	openvdb::Vec3f v1, v2,v3;
	float cthick;
	VDB::Ptr newGrid(new VDB());
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
			v3 = v2 - v1;
			float len = v3.lengthSqr();
			if (len > 1e-8 && len < 100) {
				openvdb::tools::LevelSetCapsule<openvdb::FloatGrid> factory(cthick, v1, v2);
				factory.mGrid = newGrid->grid;
				factory.rasterCapsule(resolution, 3);
			}
			cout << count++ << endl;
		}
	}
	cout << "background " <<  newGrid->grid->background() << endl;
	in.close(); 
	newGrid->floodFill();
	newGrid->isUpdated = false;
	newGrid->grid->transform().preScale(resolution);
	if(!newGrid->grid->empty())
		grids.push_back(newGrid);	
}

void ofApp::setupGui() {
	resolutionSlider.addListener(this, &ofApp::resolutionChanged);
	unionButton.addListener(this, &ofApp::doUnion);
	intersectButton.addListener(this, &ofApp::doIntersection);
	differenceButton.addListener(this, &ofApp::doDifference);
	offsetButton.addListener(this, &ofApp::doOffset);

	gui.setup();
	gui.add(resolutionSlider.setup("resolution", resolution, 0.01, 2));
	gui.add(unionButton.setup("union"));
	gui.add(intersectButton.setup("intersect"));
	gui.add(differenceButton.setup("difference"));
	gui.add(offsetButton.setup("offset"));
	gui.add(offsetSlider.setup("offset amt", 0, -5, 5));
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
					VDB::Ptr newGrid(new VDB(mesh, resolution));
					cout << "background " << newGrid->grid->background() << endl;
					grids.push_back(newGrid);
				}
			}
			else if (info.files[i].substr(info.files[i].size() - 3) == "csv") {
				loadLines(info.files[i]);
				//grid.doDifference(subBox);
				//grid.updateMesh();
				//grid.mesh.save(info.files[i].substr(0,info.files[i].size() - 3) + "ply");
			}
			else if (info.files[i].substr(info.files[i].size() - 3) == "vdb") {
				openvdb::io::File file(info.files[i]);
				cout << file.filename() << endl;
				file.open();
				openvdb::GridPtrVecPtr ogrids =  file.getGrids();
				for (auto g : *ogrids) {
					
					VDB::Ptr newGrid(new VDB());
					newGrid->grid = openvdb::gridPtrCast<openvdb::FloatGrid>(g);
					if (!newGrid->grid->empty())
						grids.push_back(newGrid);
				}
			}
		}
	}
}

//--------------------------------------------------------------
void ofApp::keyPressed(int key){
	if (key == 's') {
		grids.front()->mesh.save("solid.ply");
	} if (key == 'v') {
		openvdb::io::File file(ofToDataPath("mygrids.vdb"));
		// Add the grid pointer to a container.
		openvdb::GridPtrVec sgrids;
		for(auto &g : grids) 
			sgrids.push_back(g->grid);
		// Write out the contents of the container.
		file.write(sgrids);
		file.close();
	} else if (key == 'e') {
		grid.toEmber("cats");
	}
	else if (key == 'b') {
		for (auto g : selected) {
			g->blur();
		}
	}
	else if (key == 'o') {
		grid.offset(-.05);
	}
	else if (key == OF_KEY_DEL) {
		doDelete();
	}
}

//--------------------------------------------------------------
void ofApp::keyReleased(int key){
}

//--------------------------------------------------------------
void ofApp::mouseMoved(int x, int y ){
	isMouseClick = false;
}

//--------------------------------------------------------------
void ofApp::mouseDragged(int x, int y, int button){
	isMouseClick = false;
}

//--------------------------------------------------------------
void ofApp::mousePressed(int x, int y, int button){
	isMouseClick = true;
}

//--------------------------------------------------------------
void ofApp::mouseReleased(int x, int y, int button){
	if (isMouseClick) {
		if (button == OF_MOUSE_BUTTON_1) {
			ofVec2f mousePt(x, y);
			ofVec3f sPt = cam.screenToWorld(ofVec3f(mousePt.x, mousePt.y, -.9));
			ofVec3f sPt2 = cam.screenToWorld(ofVec3f(mousePt.x, mousePt.y, .9));
			ofVec3f srfPt;
			float t, minT = 9e9;
			VDB::Ptr sel = nullptr;
			for (auto g : grids) {
				if (g->intersectRay(sPt, (sPt2 - sPt), srfPt, t)) {
					if (t < minT) {
						minT = t;
						sel = g;
					}
				}
			}
			if (ofGetKeyPressed(OF_KEY_SHIFT)) {
				if (sel != nullptr) {
					if (!isSelected(sel)) {
						selected.push_back(sel);
					}
				}
			}
			else if (ofGetKeyPressed(OF_KEY_CONTROL)) {
				if (sel != nullptr) {
					auto it = find(selected.begin(), selected.end(), sel);
					if (it != selected.end()) {
						selected.erase(it);
					}
				}
			}
			else {
				selected.clear();
				if (sel != nullptr) selected.push_back(sel);
			}
		}
	}
	if (selected.size() > 0) {
		gumball.enable();
		//get center
		pair<ofVec3f, ofVec3f> bbox = selected.front()->bbox();
		gumball.setPosition((bbox.first + bbox.second)*0.5);
	}
	else {
		gumball.disable();
		gumball.setOrientation(ofVec3f(0, 0, 0));
	}
}

void ofApp::gumballEvent(GumballInfo & args) {
	for (auto g : selected) {
		if (args.type == GumballEventType::GUMBALL_TRANSLATE) {
			g->translate(args.dir*args.val);
		}
		else if (args.type == GumballEventType::GUMBALL_ROTATE) {
			g->translate(-gumball.getPosition());
			g->rotate(args.dir, args.val);
			g->translate(gumball.getPosition());
		}
	}
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

void ofApp::resolutionChanged(float & val) {
	resolution = val;
}


void ofApp::doUnion() {
	if (selected.size() > 1) {
		auto it = selected.begin();
		VDB::Ptr newGrid(new VDB(**it));
		it++;
		while (it != selected.end()) {
			newGrid->doUnion(**it);
			it++;
		}
		grids.push_back(newGrid);
	}
}


void ofApp::doIntersection() {
	if (selected.size() > 1) {
		auto it = selected.begin();
		VDB::Ptr newGrid(new VDB(**it));
		it++;
		while (it != selected.end()) {
			newGrid->doIntersect(**it);
			it++;
		}
		grids.push_back(newGrid);
	}
}

void ofApp::doDifference() {
	if (selected.size() > 1) {
		auto it = selected.begin();
		VDB::Ptr newGrid(new VDB(**it));
		it++;
		while (it != selected.end()) {
			newGrid->doDifference(**it);
			it++;
		}
		grids.push_back(newGrid);
	}
}

void ofApp::doOffset() {
	float offsetAmt = offsetSlider.getParameter().cast<float>().get();
	for(auto g : selected) {
		g->offset(offsetAmt);
	}
}

