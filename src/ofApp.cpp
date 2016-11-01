#include "ofApp.h"
#include "LevelSetCapsule.h"

float volume;
float defaultThickness = 0.3f;

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

//--------------------------------------------------------------
void ofApp::setup(){
	openvdb::initialize();
	gumball.setCamera(cam);
	ofAddListener(gumball.gumballEvent, this, &ofApp::gumballEvent);
	resolution = .3;
	maxTriangle = .4;
	maxError = .2;
	maskMode = true;
	maskRadius = 10;
	mask.grid->setGridClass(openvdb::v3_1_0::GridClass::GRID_FOG_VOLUME);
	isHover = false;
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

	ofDrawBitmapString(volume, ofVec2f(500, 20));
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
			//cthick = defaultThickness*0.5;
			if (len > 1e-8 && len < 100) {
				openvdb::tools::LevelSetCapsule<openvdb::FloatGrid> factory(cthick, v1, v2);
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
	gui = new ofxDatGui(ofxDatGuiAnchor::TOP_LEFT);
	gui->addHeader(":: VOXELS ::");
	ofxDatGuiSlider * resolutionSlider = gui->addSlider("resolution", 0, 2);
	gui->addBreak();
	gui->addLabel("BOOLEAN");
	gui->addButton("union");
	gui->addButton("intersection");
	gui->addButton("difference");
	gui->addBreak();
	gui->addLabel("modification");
	ofxDatGuiSlider * offsetSlider = gui->addSlider("offset amt", -10, 10);
	gui->addButton("offset");
	gui->addButton("blur");
	gui->addButton("smooth");
	gui->addButton("taubin");
	gui->addBreak();
	gui->addLabel("saving");
	ofxDatGuiSlider * triSlider = gui->addSlider("max triangle", 0, 2);
	ofxDatGuiSlider * errorSlider = gui->addSlider("max error", 0, 2);
	gui->addTextInput("filename");
	gui->addButton("process");
	gui->addButton("save");
	gui->addButton("export mesh");
	gui->addButton("clear ops");
	gui->addFolder("operations");

	resolutionSlider->bind(resolution, 0, 2);
	triSlider->bind(maxTriangle, 0, 2);
	errorSlider->bind(maxError, 0, 2);
	offsetSlider->bind(offsetAmt,-15,15);
	
	gui->onButtonEvent(this, &ofApp::buttonEvent);
	/*
	resolutionSlider.addListener(this, &ofApp::resolutionChanged);
	unionButton.addListener(this, &ofApp::doUnion);
	intersectButton.addListener(this, &ofApp::doIntersection);
	differenceButton.addListener(this, &ofApp::doDifference);
	offsetButton.addListener(this, &ofApp::doOffset);
	laplacianBlurButton.addListener(this, &ofApp::doLaplacianBlur);
	saveButton.addListener(this, &ofApp::saveVDB);
	exportButton.addListener(this, &ofApp::saveMesh);

	gui.setup();
	gui.add(resolutionSlider.setup("resolution", resolution, 0.01, 2));
	gui.add(unionButton.setup("union"));
	gui.add(intersectButton.setup("intersect"));
	gui.add(differenceButton.setup("difference"));
	gui.add(offsetButton.setup("offset"));
	gui.add(offsetSlider.setup("offset amt", 0, -5, 5));
	gui.add(laplacianBlurButton.setup("laplacian blur"));
	
	gui.add(saveButton.setup("save"));
	gui.add(exportButton.setup("export mesh"));
	*/
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
				ofMesh mesh;
				cout << info.files[i] << endl;
				mesh.load(info.files[i]);
				if (mesh.getNumIndices() > 3) {
					cout << "vertices " << mesh.getNumVertices() << endl;
					VDB::Ptr newGrid(new VDB(mesh, resolution));
					cout << "background " << newGrid->grid->background() << endl;
					grids.push_back(newGrid);
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
	isMouseClick = false;
	if (button == OF_MOUSE_BUTTON_1) {
		VDB::Ptr g;
		ofVec3f selPt = selectScreen(x, y, g, &ofApp::nullSelect);
		if (g != nullptr) {
			intersectionPt = selPt;
			isHover = true;
			mask.rasterSphere(selPt, maskRadius, 1);
			colorByMask();
		}
	}
	
}

//--------------------------------------------------------------
void ofApp::mousePressed(int x, int y, int button){
	isMouseClick = true;
}

//--------------------------------------------------------------
void ofApp::mouseReleased(int x, int y, int button){
	if(gui->hitTest(ofPoint(x,y))) return;
	if (isMouseClick) {
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

void ofApp::buttonEvent(ofxDatGuiButtonEvent e) {
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
		else if (e.target->is("clear ops")) {
			operations.clear();
			gui->getFolder("operations")->children.clear();
		}
	}
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
	//float offsetAmt = offsetSlider.getParameter().cast<float>().get();
	for(auto g : selected) {
		g->offset(offsetAmt, mask);
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
		//buildMesh(*g, m, maxTriangle, maxError);
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
			//ofMesh mesh = process(grids.back());
			//mesh.save(filename.substr(0, filename.size() - 3) + "obj");
			//grids.erase(grids.end()--);
			ofMesh mesh;
			//for (int i = 0; i < 10; ++i) grids.back()->smooth();
			buildMesh(*grids.back(), mesh, .4, 0.01);
			stringstream ss;
			ss << filename.substr(0, filename.size() - 4) << "_T" << (int)(defaultThickness * 1000) << ".obj";
			mesh.save(ss.str());
			//grids.erase(grids.end()--);
			defaultThickness += 0.05;
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