#include "ofApp.h"
#include "LevelSetCapsule.h"
#include "LevelSetWedge.h"
#include "LevelSetWedgeV.h"
#include "ObjMesh.h"

float volume;
float defaultThickness = 0.3f;
float targetVolume = 31000;// 13928.9;
float volumeEps = 100;

bool targetingVolume = false;
bool autosave = true;
float surfaceThickness = 1.0;
bool doThickening = false;
float radiusScaling = 1.0;
float minRad;

ofVboMesh loadMesh;
ofCamera prevCam;
bool importing = true;
bool autoResolution = false;
string filename;
ObjMesh loadingMesh;

float computeVolume(ofMesh & mesh) {
	float volume = 0;
	for (int i = 0; i < mesh.getNumIndices();) {
		int i1 = mesh.getIndex(i++);
		int i2 = mesh.getIndex(i++);
		int i3 = mesh.getIndex(i++);

		ofVec3f p1 = mesh.getVertex(i1);
		ofVec3f p2 = mesh.getVertex(i2);
		ofVec3f p3 = mesh.getVertex(i3);

		volume += (p1.y*p2.z - p1.z*p2.y)*p3.x + (p1.z*p2.x - p1.x*p2.z)*p3.y + (p1.x*p2.y - p1.y*p2.x)*p3.z;
	}
	return volume / 6;
}

ofVec3f min(ofVec3f & a, ofVec3f & b) {
	return ofVec3f(min(a.x, b.x), min(a.y, b.y), min(a.z, b.z));
}
ofVec3f max(ofVec3f & a, ofVec3f & b) {
	return ofVec3f(max(a.x, b.x), max(a.y, b.y), max(a.z, b.z));
}
//--------------------------------------------------------------
void ofApp::setup(){
	ofSetWindowTitle("Nervous System voxelator");
	openvdb::initialize();
	gumball.setCamera(cam);
	ofAddListener(gumball.gumballEvent, this, &ofApp::gumballEvent);
	resolution = .3;
	maxTriangle = .81;
	maxError = .02;
	maskMode = false;
	maskRadius = 10;
	mask.grid->setGridClass(openvdb::GridClass::GRID_FOG_VOLUME);
	isHover = false;
	setupGui();
	VDB::Ptr gr = thickenSrf();
	if (!gr->grid->empty()) {
		grids.push_back(gr);
	}
	else {
		cout << "empty" << endl;
	}
}

int step = 200;
string folder = "C:\\Users\\nervous system\\of_v0.8.0_vs_release\\apps\\myApps\\3DVoronoi\\bin\\data\\opt2\\lines\\SLAB_sz9.5_b10_";
//--------------------------------------------------------------
void ofApp::update() {
	if (step < 150) {
		int frame = (int)(0.5*.86*pow(step, 1.5) + 1);
		loadLines(folder + to_string(frame) + ".csv");
		grid.updateMesh();
		grid.mesh.save(folder + to_string(step) + ".ply");
		step++;
	}
}

//--------------------------------------------------------------
void ofApp::draw() {
	ofBackground(255);
	
	guiFunc();

	cam.begin();
	ofEnableDepthTest();
	//ofEnableLighting();
	ofEnableAlphaBlending();
	ofLight light0;
	light0.enable();
	ofSetColor(100);
	bool doVolume = false;
	if (grids.size() > 0) doVolume = !grids.front()->isUpdated;
	for (auto g : grids) {
		if (isSelected(g)) {
			ofSetColor(200, 0, 0);
		}
		else {
			ofSetColor(120);
		}
		g->draw();
	}
	if (doVolume) volume = computeVolume(grids.front()->mesh);

	if (importing) {
		loadMesh.draw();
	}
	ofDisableDepthTest();
	ofDisableLighting();
	gumball.draw();

	if (maskMode) {
		if (isHover) {
			//ofSetColor(200, 200, 200,100);
			//ofPushMatrix();
			//ofTranslate(intersectionPt);
			//ofDrawSphere(maskRadius);
			//ofPopMatrix();
		}
	}

	cam.end();

	gui.end();
}

void ofApp::guiFunc() {
	gui.begin();
	bool doImport = false;
	importing = false;
	ImFont * font = ImGui::GetFont();
	font->Scale = 3;
	//font->FontSize *= 3;
	ImGui::PushFont(font);
	if (ImGui::BeginMainMenuBar())
	{
		if (ImGui::BeginMenu("File"))
		{
			if (ImGui::MenuItem("Open", "CTRL+O")) {
				auto result = ofSystemLoadDialog();
				if (result.bSuccess) {
					filename = result.filePath;
					string filetype = ofToLower(filename.substr(filename.size() - 3));
					if (filetype == "obj") {
						loadingMesh.load(filename);
						doImport = true;
						prevCam = cam;
						ofVec3f bboxMin, bboxMax;
						auto & verts = loadingMesh.positions;
						bboxMin = bboxMax = verts[0];
						for (int iVert = 1; iVert < verts.size(); ++iVert) {
							bboxMin = min(bboxMin, verts[iVert]);
							bboxMax = max(bboxMax, verts[iVert]);
						}
						zoom(bboxMin, bboxMax);

					} else if (filetype == "ply" || filetype == "stl") {
						loadMesh.load(filename);
						if (loadMesh.getNumVertices()>0) {
							doImport = true;
							prevCam = cam;
							ofVec3f bboxMin, bboxMax;
							auto & verts = loadMesh.getVertices();
							bboxMin = bboxMax = verts[0];
							for (int iVert = 1; iVert < verts.size(); ++iVert) {
								bboxMin = min(bboxMin, verts[iVert]);
								bboxMax = max(bboxMax, verts[iVert]);
							}
							zoom(bboxMin, bboxMax);
						}
					}
				}
			}
			if (ImGui::MenuItem("Save", "CTRL+S")) {
				ofSystemSaveDialog("volume.vdb", "Save VDB file");
			}
			if (ImGui::MenuItem("Export", "CTRL+E")) {
				auto result = ofSystemSaveDialog("mesh.obj", "Export mesh");
				if (result.bSuccess) {
					saveMeshFast(result.filePath);
				}
			}
			ImGui::EndMenu();
		}
		if (ImGui::BeginMenu("Edit"))
		{
			if (ImGui::MenuItem("Undo", "CTRL+Z")) {}
			if (ImGui::MenuItem("Redo", "CTRL+Y", false, false)) {}  // Disabled item
			ImGui::Separator();
			if (ImGui::MenuItem("Cut", "CTRL+X")) {}
			if (ImGui::MenuItem("Copy", "CTRL+C")) {}
			if (ImGui::MenuItem("Paste", "CTRL+V")) {}
			ImGui::EndMenu();
		}
		ImGui::EndMainMenuBar();
	}

	ImGui::Begin("tools");
	ImGui::Text("boolean");
	if (ImGui::Button("union")) {
		doUnion();
	}
	if (ImGui::Button("difference")) {
		doDifference();
	}
	if (ImGui::Button("intersection")) {
		doIntersection();
	}
	ImGui::InputFloat("distance", &offsetAmt);
	if (ImGui::Button("offset")) {
		doOffset();
	}

	ImGui::Text("view");
	if(ImGui::Button("Top")) {
		cam.enableOrtho();
		cam.setPosition(cam.getTarget().getPosition() + ofVec3f(0, 0, cam.getDistance()));
		cam.lookAt(cam.getTarget(), ofVec3f(0, 1, 0));
		//cam.

	}
	if (ImGui::Button("Left")) {
	}
	if (ImGui::Button("Front")) {
	}
	if (ImGui::Button("Zoom Selected")) {
		if (selected.size() > 0) {
			ofVec3f bboxMin, bboxMax;
			auto it = selected.begin();
			auto box = (*it)->bbox();
			bboxMin = box.first;
			bboxMax = box.second;
			it++;
			while(it != selected.end()){
				box = (*it)->bbox();
				bboxMin = min(bboxMin, box.first);
				bboxMax = max(bboxMax, box.second);
				it++;

			}
			zoom(bboxMin, bboxMax);
		}
	}
	if (ImGui::Button("Zoom Extents")) {
		if (grids.size() > 0) {
			ofVec3f bboxMin, bboxMax;
			auto it = grids.begin();
			auto box = (*it)->bbox();
			bboxMin = box.first;
			bboxMax = box.second;
			it++;
			while (it != grids.end()) {
				box = (*it)->bbox();
				bboxMin = min(bboxMin, box.first);
				bboxMax = max(bboxMax, box.second);
				it++;

			}
			zoom(bboxMin, bboxMax);
		}
	}

	ImGui::End();
	
	ImGui::Begin("info");
		ImGui::Text("volume");
		ImGui::TextColored(ofColor(255, 0, 0), to_string(volume).c_str());
	ImGui::End();
	if (doImport) {
		ImGui::OpenPopup("Import Mesh");
	}
	if (ImGui::BeginPopupModal("Import Mesh")) {
		importing = true;
		
		ImGui::InputFloat("resolution", &resolution);

		ImGui::Checkbox("automatic", &autoResolution);
		if (ImGui::Button("Load Volume")) {
			VDB::Ptr newGrid(new VDB(loadMesh, resolution));
			grids.push_back(newGrid);

			ImGui::CloseCurrentPopup();
			cam.setTransformMatrix(prevCam.getGlobalTransformMatrix());
		}
		if (ImGui::Button("Close")) {
			ImGui::CloseCurrentPopup();
			cam.setTransformMatrix(prevCam.getGlobalTransformMatrix());
		}
		ImGui::EndPopup();
	}
}

void ofApp::zoom(ofVec3f bboxMin, ofVec3f bboxMax) {
	ofVec3f dir = cam.getLookAtDir();
	ofVec3f up = cam.getUpDir();
	ofVec3f center = (bboxMin + bboxMax)*0.5;
	float dist = bboxMin.distance(bboxMax);
	dir.normalize();
	cam.setTarget(center);
	cam.setPosition(center-dir*dist);
	cam.lookAt(center, up);
	cam.setNearClip(dist*.1);
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
	/*
	minRad = 9e9;
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
			minRad = min(cthick*radiusScaling, minRad);
			//cout << count++ << endl;
		}
	}

	in.clear();
	in.seekg(0);
	*/
	//resolution = min(resolution,minRad *1.2f / 3.0f);
	cout << minRad << endl;
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
			//cthick = defaultThickness*0.5;
			if (len > 1e-8) {
				openvdb::tools::LevelSetCapsule<openvdb::FloatGrid> factory(cthick*radiusScaling, v1, v2);
				factory.mGrid = newGrid->grid;
				factory.rasterCapsule(resolution, 3);
			}
			//cout << count++ << endl;
		}
	}
	cout << "background " <<  newGrid->grid->background() << endl;
	in.close(); 
	newGrid->floodFill();
	newGrid->isUpdated = false;
	newGrid->grid->pruneGrid();
	newGrid->grid->transform().preScale(resolution);
	if(!newGrid->grid->empty())
		grids.push_back(newGrid);	
}

VDB::Ptr ofApp::thickenSrf(ofMesh & mesh, float thickness) {
	openvdb::Vec3f v1, v2, v3;
	VDB::Ptr newGrid(new VDB());
	for (int i = 0; i < mesh.getNumIndices();) {
		int i1 = mesh.getIndex(i++);
		int i2 = mesh.getIndex(i++);
		int i3 = mesh.getIndex(i++);

		ofVec3f p1 = mesh.getVertex(i1);
		ofVec3f p2 = mesh.getVertex(i2);
		ofVec3f p3 = mesh.getVertex(i3);

		v1 = openvdb::Vec3f(p1.x, p1.y, p1.z);
		v2 = openvdb::Vec3f(p2.x, p2.y, p2.z);
		v3 = openvdb::Vec3f(p3.x, p3.y, p3.z);
		openvdb::tools::LevelSetWedge<openvdb::FloatGrid> factory(v1, v2,v3, thickness);
		factory.mGrid = newGrid->grid;
		factory.rasterWedge(resolution, 3);

	}
	newGrid->floodFill();
	newGrid->isUpdated = false;
	newGrid->grid->pruneGrid();
	newGrid->grid->transform().preScale(resolution);
	return newGrid;
}

VDB::Ptr ofApp::thickenSrf() {
	openvdb::Vec3f v1, v2, v3;
	VDB::Ptr newGrid(new VDB());
	v1 = openvdb::Vec3f(0,0,0);
	v2 = openvdb::Vec3f(0,10,0);
	v3 = openvdb::Vec3f(10,10,0);
	openvdb::tools::LevelSetWedgeV<openvdb::FloatGrid> factory(v1, v2, v3, 1.,2.,3.);
	//openvdb::tools::LevelSetWedge<openvdb::FloatGrid> factory(v1, v2, v3, 2.);
	factory.mGrid = newGrid->grid;
	factory.rasterWedge(resolution, 3);

	newGrid->floodFill();
	newGrid->isUpdated = false;
	newGrid->grid->pruneGrid();
	newGrid->grid->transform().preScale(resolution);
	return newGrid;
}

void ofApp::loadVol(string filename) {
	ifstream in(filename, ios_base::binary);
	//get size
	in.seekg(0, in.end);
	size_t length = in.tellg();
	in.seekg(0, in.beg);
	//assume cube and 32-bit floats

	int dim = pow(length / 4, 1.0 / 3.0)+.5;
	cout << "volume dimension: " << length << " " << dim << endl;
	
	VDB::Ptr newGrid(new VDB());
	newGrid->loadVol(in, dim, dim, dim, .1725);
	if (!newGrid->grid->empty())
		grids.push_back(newGrid);
	else
		cout << "Error: empty grid" << endl;
}

void ofApp::setupGui() {
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
			string filetype = ofToLower(info.files[i].substr(info.files[i].size() - 3));
			if (filetype == "ply" || filetype == "stl" || filetype == "obj") {
				if (doThickening) {
					ofMesh mesh;
					cout << info.files[i] << endl;
					mesh.load(info.files[i]);
					if (mesh.getNumIndices() > 3) {
						VDB::Ptr gr = thickenSrf(mesh, surfaceThickness);
						if (!gr->grid->empty()) {
							grids.push_back(gr);
						}
					}
				}
				else {
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
			} else if (filetype == "csv") {
				//loadLines(info.files[i]);
				process(info.files[i]);
				//grid.doDifference(subBox);
				//grid.updateMesh();
				//grid.mesh.save(info.files[i].substr(0,info.files[i].size() - 3) + "ply");
			}
			else if (filetype == "vdb") {
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
			else if (filetype == "vol") {
				loadVol(info.files[i]);
			}
			else {
				process(info.files[i]);
			}

		}
	}
}

void ofApp::loadFile(string filename) {
	string filetype = ofToLower(filename.substr(filename.size() - 3));
	if (filetype == "ply" || filetype == "stl" || filetype == "obj") {
		if (doThickening) {
			ofMesh mesh;
			mesh.load(filename);
			if (mesh.getNumIndices() > 3) {
				VDB::Ptr gr = thickenSrf(mesh, surfaceThickness);
				if (!gr->grid->empty()) {
					grids.push_back(gr);
				}
			}
		}
		else {
			ofMesh mesh;
			mesh.load(filename);
			if (mesh.getNumIndices() > 3) {
				cout << "vertices " << mesh.getNumVertices() << endl;
				VDB::Ptr newGrid(new VDB(mesh, resolution));
				cout << "background " << newGrid->grid->background() << endl;
				grids.push_back(newGrid);
			}
		}
	}
	else if (filetype == "csv") {
		//loadLines(info.files[i]);
		process(filename);
		//grid.doDifference(subBox);
		//grid.updateMesh();
		//grid.mesh.save(info.files[i].substr(0,info.files[i].size() - 3) + "ply");
	}
	else if (filetype == "vdb") {
		openvdb::io::File file(filename);
		file.open();
		openvdb::GridPtrVecPtr ogrids = file.getGrids();
		for (auto g : *ogrids) {

			VDB::Ptr newGrid(new VDB());
			newGrid->grid = openvdb::gridPtrCast<openvdb::FloatGrid>(g);
			if (!newGrid->grid->empty())
				grids.push_back(newGrid);
		}
	}
	else if (filetype == "vol") {
		loadVol(filename);
	}
}
//--------------------------------------------------------------
void ofApp::keyPressed(int key){
	if (key == 's') {
		//grids.front()->mesh.save("solid.ply");
	} if (key == 'v') {
		
	} else if (key == 'e') {
		//grid.toEmber("cats");
	}
	else if (key == 'b') {
		
	}
	else if (key == 'o') {
		//grid.offset(-.05);
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
	if (maskMode) {
		VDB::Ptr g;
		ofVec3f selPt = selectScreen(x, y, g, &ofApp::nullSelect);
		if (g != nullptr) {
			intersectionPt = selPt;
			isHover = true;
		}
		else {
			isHover = false;
		}
	}
}

//--------------------------------------------------------------
void ofApp::mouseDragged(int x, int y, int button){
	if (!ImGui::GetIO().WantCaptureMouse) {
		isMouseClick = false;
		if (button == OF_MOUSE_BUTTON_1) {
			/*
			VDB::Ptr g;
			ofVec3f selPt = selectScreen(x, y, g, &ofApp::nullSelect);
			if (g != nullptr) {
				intersectionPt = selPt;
				isHover = true;
				mask.rasterSphere(selPt, maskRadius, 1);
				colorByMask();
			}
			*/
		}
	}
	
}

//--------------------------------------------------------------
void ofApp::mousePressed(int x, int y, int button){
	isMouseClick = true;
}

//--------------------------------------------------------------
void ofApp::mouseReleased(int x, int y, int button){
	//if(gui->hitTest(ofPoint(x,y))) return;
	if (isMouseClick && !ImGui::GetIO().WantCaptureMouse) {
		if (button == OF_MOUSE_BUTTON_1) {
			VDB::Ptr sel = nullptr;
			if (ofGetKeyPressed(OF_KEY_CONTROL)) {
				selectScreen(x, y, sel, &ofApp::getSelect);	
			}
			else {
				selectScreen(x, y, sel, &ofApp::getUnselect);
			}

			if (ofGetKeyPressed(OF_KEY_SHIFT)) {
				if (sel != nullptr) {
					selected.push_back(sel);
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

ofVec3f ofApp::selectScreen(float x, float y, VDB::Ptr & out, bool(ofApp::*selFunc)(VDB::Ptr)) {
	ofVec2f mousePt(x, y);
	ofVec3f sPt = cam.screenToWorld(ofVec3f(mousePt.x, mousePt.y, -.9));
	ofVec3f sPt2 = cam.screenToWorld(ofVec3f(mousePt.x, mousePt.y, .9));
	ofVec3f srfPt, selPt;
	float t, minT = 9e9;
	VDB::Ptr sel = nullptr;
	for (auto g : grids) {
		if ((this->*selFunc)(g)) {
			if (g->intersectRay(sPt, (sPt2 - sPt), srfPt, t)) {

				if (t < minT) {
					minT = t;
					sel = g;
					selPt = srfPt;
				}
			}
		}
	}
	out = sel;
	return selPt;
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

void ofApp::colorByMask() {
	ofColor c1(120);
	ofColor c2(255, 0, 0);
	openvdb::FloatGrid::ConstAccessor accessor = mask.grid->getConstAccessor();
	// Instantiate the GridSampler template on the accessor type and on
	// a box sampler for accelerated trilinear interpolation.
	openvdb::tools::GridSampler<openvdb::FloatGrid::ConstAccessor, openvdb::tools::BoxSampler>
		sampler(accessor, mask.grid->transform());

	//openvdb::tools::GridSampler<openvdb::FloatGrid, openvdb::tools::BoxSampler> sampler(*mask.grid);
	for (auto grid : grids) {
		ofVboMesh & mesh = grid->mesh;
		mesh.enableColors();
		vector<ofFloatColor> & colors = mesh.getColors();
		colors.resize(mesh.getNumVertices());
		for (int i = 0, l = mesh.getNumVertices(); i < l; ++i) {
			//float v = mask.samplePt(grid->mesh.getVertex(i));
			ofVec3f pt = mesh.getVertex(i);
			float v = sampler.wsSample(openvdb::Vec3R(pt.x, pt.y, pt.z));
			colors[i] = c1.getLerped(c2, v);
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
/*
void ofApp::buttonEvent(ofxDatGuiButtonEvent e) {
	cout << e.target->getName() << endl;
	if (ofGetKeyPressed(OF_KEY_CONTROL)) {
		if (e.target->is("union")) {
			for(auto sel : selected) operations.push_back(MeshOp("union", sel));
			gui->getFolder("operations")->addLabel("union");
		}
		else if (e.target->is("difference")) {
			for (auto sel : selected) operations.push_back(MeshOp("difference", sel));
			gui->getFolder("operations")->addLabel("difference");
		}
		else if (e.target->is("offset")) {
			operations.push_back(MeshOp("offset", 0));
			gui->getFolder("operations")->addLabel("offset");
		}
		else if (e.target->is("blur")) {
			operations.push_back(MeshOp("blur"));
			gui->getFolder("operations")->addLabel("blur");

		}
		else if (e.target->is("smooth")) {
			operations.push_back(MeshOp("smooth"));
			gui->getFolder("operations")->addLabel("smooth");

		}
	} else {
		if (e.target->is("union")) {
			doUnion();
		}
		else if (e.target->is("intersection")) {
			doIntersection();
		}
		else if (e.target->is("difference")) {
			doDifference();
		}
		else if (e.target->is("offset")) {
			doOffset();
		}
		else if (e.target->is("blur")) {
			doLaplacianBlur();
		}
		else if (e.target->is("smooth")) {
			doSmooth();
		}
		else if (e.target->is("taubin")) {
			doTaubin();
		}
		else if (e.target->is("save")) {
			saveVDB(gui->getTextInput("filename")->getText());
		}
		else if (e.target->is("process")) {
			process(gui->getTextInput("filename")->getText());
		}
		else if (e.target->is("export mesh")) {
			saveMesh(gui->getTextInput("filename")->getText());
		}
		else if (e.target->is("export mesh fast")) {
			saveMeshFast(gui->getTextInput("filename")->getText());
		}
		else if (e.target->is("clear ops")) {
			operations.clear();
			gui->getFolder("operations")->children.clear();
		}
		else if (e.target->is("thickening")) {
			doThickening = ((ofxDatGuiToggle *)e.target)->getEnabled();
		}
	}
}
*/
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
	//float offsetAmt = offsetSlider.getParameter().cast<float>().get();
	for(auto g : selected) {
		g->offset(-offsetAmt, mask);
	}
}

void ofApp::doLaplacianBlur() {
	for (auto g : selected) {
		g->blur();
	}
}

void ofApp::doSmooth() {
	for (auto g : selected) {
		for (int i = 0; i < 10;++i)
			g->smooth();
	}
}

void ofApp::doTaubin() {
	for (auto g : selected) {
		g->taubin();
	}
}

void ofApp::saveVDB(string filename) {
	if (filename.size() < 5) return;
	openvdb::io::File file(ofToDataPath(filename));
	// Add the grid pointer to a container.
	openvdb::GridPtrVec sgrids;
	for (auto &g : grids)
		sgrids.push_back(g->grid);
	// Write out the contents of the container.
	file.write(sgrids);
	file.close();
}

void ofApp::saveMesh(string filename) {
	if (filename.size() < 5) return;
	ofMesh m;
	cout << filename << endl;
	for (auto g : selected) {
		//g->updateMesh();
		buildMesh(*g, m, maxTriangle, maxError);
		//ofIndexType baseIndex = m.getNumVertices();
		//m.addVertices(m.getVertices());
		//for (auto i : g->mesh.getIndices()) {
		//	m.addIndex(i + baseIndex);
		//}
	}
	m.save(filename);
}


void ofApp::saveMeshFast(string filename) {
	if (filename.size() < 5) return;
	ofMesh m;
	cout << filename << endl;
	for (auto g : selected) {
		g->updateMesh();
		ofIndexType baseIndex = m.getNumVertices();
		m.addVertices(g->mesh.getVertices());
		for (auto i : g->mesh.getIndices()) {
			m.addIndex(i + baseIndex);
		}
	}
	m.save(filename);
}

void ofApp::process(string filename) {
	string filetype = ofToLower(filename.substr(filename.size() - 3));
	if (filetype == "csv" || filetype == "txt") {
		defaultThickness = .3;
		//while (defaultThickness <= .71) {
			loadLines(filename);
			for (int i = 0; i < 3; ++i) grids.back()->smooth();
			/*
			grids.back()->updateMesh();
			volume = computeVolume(grids.back()->mesh);
			while (targetingVolume && abs(volume - targetVolume) > volumeEps) {
				radiusScaling *= (sqrt(targetVolume / volume)-1.0)*0.75+1.0;
				grids.erase(--grids.end());
				loadLines(filename);
				for (int i = 0; i < 8; ++i) grids.back()->smooth();
				grids.back()->updateMesh();
				volume = computeVolume(grids.back()->mesh);
			}
			*/
			//ofMesh mesh = process(grids.back());
			//mesh.save(filename.substr(0, filename.size() - 3) + "obj");
			//grids.erase(grids.end()--);
			if (autosave) {
				ofMesh mesh;
				//for (int i = 0; i < 10; ++i) grids.back()->smooth();
				buildMesh(*grids.back(), mesh, maxTriangle, maxError);
				stringstream ss;
				volume = computeVolume(mesh);
				ss << filename.substr(0, filename.size() - 4) << "_T" << (int)(minRad * 1000) << "_V" << ((int)(volume * 10)) / 10.0 << ".obj";
				if (volume < 0) {
					//flip
					for (int i = 0; i < mesh.getNumIndices(); i += 3) {
						unsigned int temp = mesh.getIndex(i + 1);
						mesh.setIndex(i + 1, mesh.getIndex(i + 2));
						mesh.setIndex(i + 2, temp);
					}
				}
				mesh.save(ss.str());
				if (grids.size() > 0) {
					auto gEnd = grids.end();
					gEnd--;
					grids.erase(gEnd);
				}
				//defaultThickness += 0.05;
			}
		//}
	}
}

void ofApp::doOp(VDB::Ptr grid, MeshOp & op) {
	if (op.name == "smooth") {
		for (int i = 0; i < 10;++i) grid->smooth();
	}
	else if (op.name == "blur") {
		grid->blur();
	}
	else if (op.name == "union") {
		grid->doUnion(*op.ptr);
	}
	else if (op.name == "difference") {
		grid->doDifference(*op.ptr);

	}
	else if (op.name == "offset") {
		grid->offset(op.val);
	}
}
ofMesh ofApp::process(VDB::Ptr grid) {
	//operations
	//do operations
	for (auto & op : operations) {
		doOp(grid, op);
	}
	ofMesh out;
	buildMesh(*grid, out, maxTriangle, maxError);
	return out;
}