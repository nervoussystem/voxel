#pragma once

#include "ofEvents.h"
#include "ofNode.h"
#include "ofGraphics.h"
#include "ofCamera.h"

enum GumballEventType {
	GUMBALL_DRAGGING,
	GUMBALL_DOWN,
	GUMBALL_UP,
};

struct GumballInfo {
	ofMatrix4x4 transform;
	GumballEventType type;

	GumballInfo(ofMatrix4x4 & trans, GumballEventType t) {
		trans = transform;
		type = t;
	}

	GumballInfo() {
		type = GUMBALL_DRAGGING;
	}
};

class Gumball : public ofNode {

public:
	void draw();
	Gumball();

	void disable();
	void enable();

	void setCamera(ofCamera & cam);

	ofEvent<GumballInfo> gumballEvent;
private:
	float size;
	float scaling;
	float mouseSensitivity;
	void update(ofEventArgs & args);
	ofCamera * camera;
	bool enabled;

	bool bXClicked, bYClicked,bZClicked;
	bool bXHover, bYHover, bZHover;
	//arcs
	bool bXYClicked, bYZClicked,bZXClicked;
	bool bXYHover, bYZHover, bZXHover;

	bool bClicked;
	ofVec2f pMouse;
	GumballInfo gumballInfo;
};