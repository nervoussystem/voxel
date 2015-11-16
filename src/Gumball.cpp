#include "Gumball.h"
#include "ofPolyline.h"

Gumball::Gumball() {
	size = 35;
	mouseSensitivity = 4;
	enabled = true;
	bXClicked = false;
	bYClicked = false;
	bZClicked = false;
	bClicked = false;
	bXHover = bYHover = bZHover = false;
	bXYHover = bYZHover = bZXHover = false;
	bXYClicked = bYZClicked = bZXClicked = false;
	scaling = 1;
	ofAddListener(ofEvents().update , this, &Gumball::update);
}

ofPolyline arcShape;
void drawArc(float startAngle, float endAngle, float radius, ofVec3f center, ofVec3f axis1 = ofVec3f(1,0,0),ofVec3f axis2 = ofVec3f(0,1,0)) {
	int divs = abs(endAngle-startAngle)/(PI*.03);
	arcShape.clear();
	for(int i=0;i<divs;++i) {
		float angle = ofLerp(startAngle,endAngle, i*1.0/(divs-1.0));
		float cosA = cos(angle);
		float sinA = sin(angle);

		arcShape.addVertex(center+radius*cosA*axis1+radius*sinA*axis2);

	}
	arcShape.draw();
}

void Gumball::draw() {
	if(enabled) {
		float offset = size*.2;

		//x
		ofSetColor(255,0,0);
		ofSetLineWidth(1);
		if(bXHover || bXClicked) {
			ofSetLineWidth(3);
		}
		ofLine(getPosition()+scaling*offset*getXAxis(),getPosition()+scaling*(size)*getXAxis());

		ofSetLineWidth(1);
		if(bYZHover || bYZClicked) {
			ofSetLineWidth(3);
		}
		drawArc(PI*.05,PI*.45,size*scaling, getPosition(), getYAxis(), getZAxis());

		//y
		ofSetColor(0,255,0);
		ofSetLineWidth(1);
		if(bYHover || bYClicked) {
			ofSetLineWidth(3);
		}
		ofLine(getPosition()+scaling*offset*getYAxis(),getPosition()+scaling*(size)*getYAxis());

		ofSetLineWidth(1);
		if(bZXHover || bZXClicked) {
			ofSetLineWidth(3);
		}
		drawArc(PI*.05,PI*.45,size*scaling, getPosition(), getZAxis(), getXAxis());

		//z
		ofSetLineWidth(1);
		if(bZHover || bZClicked) {
			ofSetLineWidth(3);
		}
		ofSetColor(0,0,255);
		ofLine(getPosition()+scaling*offset*getZAxis(),getPosition()+scaling*(size)*getZAxis());

		ofSetLineWidth(1);
		if(bXYHover || bXYClicked) {
			ofSetLineWidth(3);
		}
		drawArc(PI*.05,PI*.45,size*scaling, getPosition(), getXAxis(), getYAxis());

	}
}

void Gumball::enable() {
	enabled = true;
	ofAddListener(ofEvents().update , this, &Gumball::update);
}

void Gumball::disable() {
	enabled = false;
	ofRemoveListener(ofEvents().update, this, &Gumball::update);
}

void Gumball::setCamera(ofCamera & cam) {
	camera = &cam;
}

//return normalized coordinate with 0=p1 and 1=p2
float projectPtToSegment(const ofVec2f & p1, const ofVec2f & p2, const ofVec2f & p3) {
	ofVec2f v = p2-p1;
	float dot = v.dot(p3-p1);
	dot /= v.lengthSquared();
	return dot;
}

float sqDistance(const ofVec2f & p1, const ofVec2f & p2, const ofVec2f & p3) {
	//project
	ofVec2f v = p2-p1;
	float dot = v.dot(p3-p1);
	dot /= v.lengthSquared();
	dot = ofClamp(dot,0,1);
	v = v*dot+p1;

	return v.distanceSquared(p3);
}

ofVec3f linePlaneIntersect(const ofVec3f &p1, const ofVec3f & p2, const ofVec3f &center, const ofVec3f & norm) {
	float d = (center-p1).dot(norm)/(p2-p1).dot(norm);
	return p1+(p2-p1)*d;
}

void Gumball::update(ofEventArgs &args) {
	if(enabled) {
		ofVec3f screenPos = camera->worldToScreen(getPosition());
		ofVec3f modPos = camera->screenToWorld(screenPos+ofVec3f(1,0,0));
		ofVec2f mouse(ofGetMouseX(),ofGetMouseY());
		scaling = modPos.distance(getPosition());
		if(!ofGetMousePressed(OF_MOUSE_BUTTON_LEFT) || !bClicked) {
			float offset = size*.2;

			//do hover
			bXHover = bYHover = bZHover = false;
			bXYHover = bYZHover = bZXHover = false;
			if(screenPos.distanceSquared(mouse) < size*size) {
				//get distance to axis
				float minDist = mouseSensitivity*mouseSensitivity+1;
				ofVec2f px1 = camera->worldToScreen(getPosition()+scaling*offset*getXAxis());
				ofVec2f px2 = camera->worldToScreen(getPosition()+scaling*(size)*getXAxis());
				float d = sqDistance(px1,px2,mouse);
				if(d < minDist) {
					minDist = d;
					bXHover = true;
				}

				ofVec2f py1 = camera->worldToScreen(getPosition()+scaling*offset*getYAxis());
				ofVec2f py2 = camera->worldToScreen(getPosition()+scaling*(size)*getYAxis());
				d = sqDistance(py1,py2,mouse);
				if(d < minDist) {
					bXHover = false;
					bYHover = true;
					minDist = d;
				}

				ofVec2f pz1 = camera->worldToScreen(getPosition()+scaling*offset*getZAxis());
				ofVec2f pz2 = camera->worldToScreen(getPosition()+scaling*(size)*getZAxis());
				d = sqDistance(pz1,pz2,mouse);
				if(d < minDist) {
					bXHover = bYHover = false;
					bZHover = true;
					minDist = d;
				}

				//if no hit on axis go to arcs
				if(minDist == mouseSensitivity*mouseSensitivity+1) {
					/*ofVec2f xPos = camera->worldToScreen(getPosition()+getXAxis());
					float xDir = projectPtToSegment(screenPos, xPos, mouse);

					ofVec2f yPos = camera->worldToScreen(getPosition()+getYAxis());
					float yDir = projectPtToSegment(screenPos, yPos, mouse);

					ofVec2f zPos = camera->worldToScreen(getPosition()+getZAxis());
					float zDir = projectPtToSegment(screenPos, zPos, mouse);
					*/

					//get plane intersection
					ofVec3f p1 = camera->screenToWorld(ofVec3f(mouse.x,mouse.y,0));
					ofVec3f p2 = camera->screenToWorld(ofVec3f(mouse.x,mouse.y,5));
					ofVec3f intersect = linePlaneIntersect(p1,p2,getPosition(), getZAxis());
					float xDir = getXAxis().dot(intersect-getPosition());
					float yDir = getYAxis().dot(intersect-getPosition());
					//inside quadrant
					if(xDir > 0 && yDir > 0) {
						float sqRad = xDir*xDir+yDir*yDir;
						if(sqRad > pow(size*scaling-mouseSensitivity*scaling,2) && sqRad < pow(size*scaling+mouseSensitivity*scaling,2)) {
							bXYHover = true;
						}
					}
					if(!bXYHover) {
						intersect = linePlaneIntersect(p1,p2,getPosition(), getXAxis());
						float zDir = getZAxis().dot(intersect-getPosition());
						yDir = getYAxis().dot(intersect-getPosition());
						//inside quadrant
						if(zDir > 0 && yDir > 0) {
							float sqRad = zDir*zDir+yDir*yDir;
							if(sqRad > pow(size*scaling-mouseSensitivity*scaling,2) && sqRad < pow(size*scaling+mouseSensitivity*scaling,2)) {
								bYZHover = true;
							}
						}
						if(!bYZHover) {
							intersect = linePlaneIntersect(p1,p2,getPosition(), getYAxis());
							zDir = getZAxis().dot(intersect-getPosition());
							xDir = getXAxis().dot(intersect-getPosition());
							//inside quadrant
							if(zDir > 0 && xDir > 0) {
								float sqRad = zDir*zDir+xDir*xDir;
								if(sqRad > pow(size*scaling-mouseSensitivity*scaling,2) && sqRad < pow(size*scaling+mouseSensitivity*scaling,2)) {
									bZXHover = true;
								}
							}
						}
					}

				}

			}
		}

		if(ofGetMousePressed(OF_MOUSE_BUTTON_LEFT)) {
			if(!bClicked) {
				bXClicked = bXHover;
				bYClicked = bYHover;
				bZClicked = bZHover;
				bXYClicked = bXYHover;
				bYZClicked = bYZHover;
				bZXClicked = bZXHover;
				bClicked = true;
			} else {
				ofMatrix4x4 transMat;

				if(bXClicked) {
					ofVec2f pt = camera->worldToScreen(getPosition()+getXAxis());
					float proj = projectPtToSegment(screenPos, pt, mouse);
					float prevProj = projectPtToSegment(screenPos, pt, pMouse);
					transMat.setTranslation(getXAxis()*(proj-prevProj));
					move(getXAxis()*(proj-prevProj));
					gumballInfo.transform = transMat;
					gumballInfo.type = GumballEventType::GUMBALL_TRANSLATE;
					gumballInfo.dir = getXAxis();
					gumballInfo.val = proj - prevProj;
					ofNotifyEvent(gumballEvent, gumballInfo ,this);
				} else if(bYClicked) {
					ofVec2f pt = camera->worldToScreen(getPosition()+getYAxis());
					float proj = projectPtToSegment(screenPos, pt, mouse);
					float prevProj = projectPtToSegment(screenPos, pt, pMouse);
					transMat.setTranslation(getYAxis()*(proj-prevProj));
					move(getYAxis()*(proj-prevProj));
					gumballInfo.transform = transMat;
					gumballInfo.type = GumballEventType::GUMBALL_TRANSLATE;
					gumballInfo.dir = getYAxis();
					gumballInfo.val = proj - prevProj;
					ofNotifyEvent(gumballEvent, gumballInfo ,this);
				} else if(bZClicked) {
					ofVec2f pt = camera->worldToScreen(getPosition()+getZAxis());
					float proj = projectPtToSegment(screenPos, pt, mouse);
					float prevProj = projectPtToSegment(screenPos, pt, pMouse);
					transMat.setTranslation(getZAxis()*(proj-prevProj));
					move(getZAxis()*(proj-prevProj));
					gumballInfo.transform = transMat;
					gumballInfo.type = GumballEventType::GUMBALL_TRANSLATE;
					gumballInfo.dir = getZAxis();
					gumballInfo.val = proj - prevProj;
					ofNotifyEvent(gumballEvent, gumballInfo ,this);
				} else if(bXYClicked) {
					if(mouse != pMouse) {
						ofVec3f p1 = camera->screenToWorld(ofVec3f(mouse.x,mouse.y,0));
						ofVec3f p2 = camera->screenToWorld(ofVec3f(mouse.x,mouse.y,5));
						ofVec3f intersect = linePlaneIntersect(p1,p2,getPosition(), getZAxis());
						float xDir = getXAxis().dot(intersect-getPosition());
						float yDir = getYAxis().dot(intersect-getPosition());

						p1 = camera->screenToWorld(ofVec3f(pMouse.x,pMouse.y,0));
						p2 = camera->screenToWorld(ofVec3f(pMouse.x,pMouse.y,5));
						intersect = linePlaneIntersect(p1,p2,getPosition(), getZAxis());
						float pxDir = getXAxis().dot(intersect-getPosition());
						float pyDir = getYAxis().dot(intersect-getPosition());

						float radians = atan2(yDir,xDir)-atan2(pyDir,pxDir);
						ofQuaternion rot(radians*180/PI,getZAxis());
						rotate(rot);
						transMat.translate(-getPosition());
						transMat.rotate(rot);
						transMat.translate(getPosition());
						gumballInfo.transform = transMat;
						gumballInfo.type = GumballEventType::GUMBALL_ROTATE;
						gumballInfo.dir = getZAxis();
						gumballInfo.val = radians;
						ofNotifyEvent(gumballEvent, gumballInfo ,this);
					}
				} else if(bYZClicked) {
					if(mouse != pMouse) {
						ofVec3f p1 = camera->screenToWorld(ofVec3f(mouse.x,mouse.y,0));
						ofVec3f p2 = camera->screenToWorld(ofVec3f(mouse.x,mouse.y,5));
						ofVec3f intersect = linePlaneIntersect(p1,p2,getPosition(), getXAxis());
						float zDir = getZAxis().dot(intersect-getPosition());
						float yDir = getYAxis().dot(intersect-getPosition());

						p1 = camera->screenToWorld(ofVec3f(pMouse.x,pMouse.y,0));
						p2 = camera->screenToWorld(ofVec3f(pMouse.x,pMouse.y,5));
						intersect = linePlaneIntersect(p1,p2,getPosition(), getXAxis());
						float pzDir = getZAxis().dot(intersect-getPosition());
						float pyDir = getYAxis().dot(intersect-getPosition());

						float radians = atan2(zDir,yDir)-atan2(pzDir,pyDir);
						ofQuaternion rot(radians*180/PI,getXAxis());
						rotate(rot);
						transMat.translate(-getPosition());
						transMat.rotate(rot);
						transMat.translate(getPosition());
						gumballInfo.transform = transMat;
						gumballInfo.type = GumballEventType::GUMBALL_ROTATE;
						gumballInfo.dir = getXAxis();
						gumballInfo.val = radians;
						ofNotifyEvent(gumballEvent, gumballInfo ,this);
					}
				} else if(bZXClicked) {
					if(mouse != pMouse) {
						ofVec3f p1 = camera->screenToWorld(ofVec3f(mouse.x,mouse.y,0));
						ofVec3f p2 = camera->screenToWorld(ofVec3f(mouse.x,mouse.y,5));
						ofVec3f intersect = linePlaneIntersect(p1,p2,getPosition(), getYAxis());
						float xDir = getXAxis().dot(intersect-getPosition());
						float zDir = getZAxis().dot(intersect-getPosition());

						p1 = camera->screenToWorld(ofVec3f(pMouse.x,pMouse.y,0));
						p2 = camera->screenToWorld(ofVec3f(pMouse.x,pMouse.y,5));
						intersect = linePlaneIntersect(p1,p2,getPosition(), getYAxis());
						float pxDir = getXAxis().dot(intersect-getPosition());
						float pzDir = getZAxis().dot(intersect-getPosition());

						float radians = atan2(xDir,zDir)-atan2(pxDir,pzDir);
						ofQuaternion rot(radians*180/PI,getYAxis());
						rotate(rot);
						transMat.translate(-getPosition());
						transMat.rotate(rot);
						transMat.translate(getPosition());

						gumballInfo.transform = transMat;
						gumballInfo.type = GumballEventType::GUMBALL_ROTATE;
						gumballInfo.dir = getYAxis();
						gumballInfo.val = radians;

						ofNotifyEvent(gumballEvent, gumballInfo ,this);
					}
				}
			}
		} else {
			bClicked = bXClicked = bYClicked = bZClicked = bXYClicked = bYZClicked = bZXClicked = false;
		}
		pMouse = mouse;
	}
}