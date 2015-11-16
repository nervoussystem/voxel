#pragma once

#include "ofMain.h"
#include "VDB.h"
#include "Gumball.h"
#include "ofxGui.h"

class ofApp : public ofBaseApp{

	public:
		ofEasyCam cam;
		Gumball gumball;
		VDB grid;
		VDB mask;

		list<VDB::Ptr> grids;
		list<VDB::Ptr> selected;
		stack<list<VDB *> > state;

		float resolution;
		void setup();
		void update();
		void draw();
		void loadLines(string filename);
		bool isSelected(VDB::Ptr g);

		void doDelete();

		ofxPanel gui;
		ofxButton unionButton;
		ofxButton differenceButton;
		ofxButton intersectButton;
		ofxButton offsetButton;
		ofxFloatSlider offsetSlider;
		ofxFloatSlider resolutionSlider;
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

		void resolutionChanged(float & val);
		void doUnion();
		void doIntersection();
		void doDifference();
		void doOffset();
};
