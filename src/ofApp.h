#pragma once

#include "ofMain.h"
#include "VDB.h"

class ofApp : public ofBaseApp{

	public:
		ofEasyCam cam;

		VDB grid;

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
		
		void keyPressed(int key);
		void keyReleased(int key);
		void mouseMoved(int x, int y );
		void mouseDragged(int x, int y, int button);
		void mousePressed(int x, int y, int button);
		void mouseReleased(int x, int y, int button);
		void mouseEntered(int x, int y);
		void mouseExited(int x, int y);
		void dragEvent(ofDragInfo info);
		void windowResized(int w, int h);
		void gotMessage(ofMessage msg);
		
};
