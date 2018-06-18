#include "ofMain.h"
#include "ofApp.h"
//#include "../resource.h"

//========================================================================
int main( ){
//int WINAPI WinMain(HINSTANCE hInstance, HINSTANCE hPrevInstance, LPSTR lpCmdLine, int nShowCmd) {
	ofSetupOpenGL(1024*2, 768*2, OF_WINDOW);			// <-------- setup the GL context
	//HWND hwnd = ofGetWin32Window();
	//HICON hMyIcon = LoadIcon(hInstance, MAKEINTRESOURCE(IDI_ICON1));// = LoadIcon(hInstance, MAKEINTRESOURCE(IDI_ICON1));
	//SendMessage(hwnd, WM_SETICON, ICON_BIG, (LPARAM)hMyIcon);
	//SetIcon(hMyIcon, FALSE);
	// this kicks off the running of my app
	// can be OF_WINDOW or OF_FULLSCREEN
	// pass in width and height too:
	ofRunApp(new ofApp());

}
