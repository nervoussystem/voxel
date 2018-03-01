#pragma once

#include "ofMain.h"
#include "ofxGui.h"
#include "ofxDatGui.h"
#include "VDB.h"
#include "Meshing.h"
#include "Gumball.h"

struct MeshOp {
	string name;
	VDB::Ptr ptr;
	float val;
	MeshOp() {

	}
	MeshOp(string _name, VDB::Ptr _ptr = nullptr) {
		name = _name;
		ptr = _ptr;
	}

	MeshOp(string _name, float _val) {
		name = _name;
		ptr = nullptr;
		val = _val;
	}
};

class ofApp : public ofBaseApp{

	public:
		ofEasyCam cam;
		Gumball gumball;
		VDB grid;
		VDB mask;
		bool maskMode;
		float maskRadius;
		bool isHover;
		ofVec3f intersectionPt;
		float maxTriangle;
		float maxError;

		list<VDB::Ptr> grids;
		list<VDB::Ptr> selected;
		stack<list<VDB *> > state;

		vector<MeshOp> operations;
		float resolution;
		float offsetAmt;
		void setup();
		void update();
		void draw();
		void loadLines(string filename);
		void loadVol(string filename);
		VDB::Ptr thickenSrf(ofMesh & mesh, float thickness);
		bool isSelected(VDB::Ptr g);

		void colorByMask();

		void doDelete();

		ofxDatGui * gui;
		//ofxPanel gui;
		ofxButton unionButton;
		ofxButton differenceButton;
		ofxButton intersectButton;
		ofxButton offsetButton;
		ofxFloatSlider offsetSlider;
		ofxFloatSlider resolutionSlider;
		ofxButton laplacianBlurButton;
		ofxButton saveButton;
		ofxButton exportButton;
		
		ofVec3f selectScreen(float x, float y, VDB::Ptr & out, bool(ofApp::*selFunc)(VDB::Ptr));
		void setupGui();
		bool isMouseClick;
		
		void keyPressed(int key);
		void keyReleased(int key);
		void mouseMoved(int x, int y );
		void mouseDragged(int x, int y, int button);
		void mousePressed(int x, int y, int button);
		void mouseReleased(int x, int y, int button);
		void mouseEntered(int x, int y);
		void mouseExited(int x, int y);
		void gumballEvent(GumballInfo & args);
		void dragEvent(ofDragInfo info);
		void windowResized(int w, int h);
		void gotMessage(ofMessage msg);

		void buttonEvent(ofxDatGuiButtonEvent e);
		void resolutionChanged(float & val);
		void doUnion();
		void doIntersection();
		void doDifference();
		void doOffset();
		void saveVDB(string filename);
		void saveMesh(string filename);
		void saveMeshFast(string filename);
		void doLaplacianBlur();
		void doSmooth();
		void doTaubin();

		void process(string filename);
		ofMesh process(VDB::Ptr grid);
		void doOp(VDB::Ptr grid, MeshOp & op);
		bool nullSelect(VDB::Ptr g) { return true; }
		bool getUnselect(VDB::Ptr g) { return !isSelected(g); }
		bool getSelect(VDB::Ptr g) { return isSelected(g); }

};
