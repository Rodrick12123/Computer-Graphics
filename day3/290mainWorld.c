


/* On macOS, compile with...
    clang 290mainWorld.c 040pixel.o -lglfw -framework OpenGL -framework Cocoa -framework IOKit
*/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <GLFW/glfw3.h>

#include "040pixel.h"

#include "250vector.c"
#include "280matrix.c"
#include "150texture.c"
#include "260shading.c"
#include "260depth.c"
#include "270triangle.c"
#include "280mesh.c"
#include "190mesh2D.c"
#include "250mesh3D.c"
#include "280camera.c"

#define ATTRX 0
#define ATTRY 1
#define ATTRZ 2
#define ATTRS 3
#define ATTRT 4
#define ATTRN 5
#define ATTRO 6
#define ATTRP 7
#define VARYX 0
#define VARYY 1
#define VARYZ 2
#define VARYW 3
#define VARYS 4
#define VARYT 5
#define UNIFR 0
#define UNIFG 1
#define UNIFB 2
#define UNIFMODELING 3
#define UNIFPROJ 19
#define TEXR 0
#define TEXG 1
#define TEXB 2

void shadeVertex(
        int unifDim, const double unif[], int attrDim, const double attr[], 
        int varyDim, double vary[]) {
	double attrHomog[4] = {attr[0], attr[1], attr[2], 1.0};
	double world[4];
	
	mat441Multiply((double(*)[4])(&unif[UNIFMODELING]), attrHomog, world);
	mat441Multiply((double(*)[4])(&unif[UNIFPROJ]), world, vary);

	vary[VARYS] = attr[ATTRS];
	vary[VARYT] = attr[ATTRT];
}



void shadeFragment(
        int unifDim, const double unif[], int texNum, const texTexture *tex[], 
        int varyDim, const double vary[], double rgbd[4]) {
	double sample[tex[1]->texelDim];
	texSample(tex[1], vary[VARYS], vary[VARYT], sample);
	rgbd[0] = 1.0;//sample[TEXR] * unif[UNIFR];
	rgbd[1] = 1.0; //sample[TEXG] * unif[UNIFG];
	rgbd[2] = 0.0;//sample[TEXB] * unif[UNIFB];
	rgbd[3] = vary[VARYZ];
}

void shadeFragment2(
        int unifDim, const double unif[], int texNum, const texTexture *tex[], 
        int varyDim, const double vary[], double rgbd[4]) {
	double sample[tex[0]->texelDim];
	texSample(tex[0], vary[VARYS], vary[VARYT], sample);
	rgbd[0] = sample[TEXR];
	rgbd[1] = sample[TEXG]; 
	rgbd[2] = sample[TEXB]; 
	rgbd[3] = vary[VARYZ];
}

shaShading sha;
shaShading sha2;
texTexture texture;
texTexture texture2;
const texTexture *textures[2] = {&texture, &texture2};
const texTexture **tex = textures;
meshMesh mesh;
meshMesh mesh2;
depthBuffer buf;
double unif[3 + 16 +16] = {
    1.0, 1.0, 1.0, 
	1.0, 0.0, 0.0, 0.0, 
	0.0, 1.0, 0.0, 0.0, 
	0.0, 0.0, 1.0, 0.0, 
	0.0, 0.0, 0.0, 1.0,
	
	1.0, 0.0, 0.0, 0.0, 
	0.0, 1.0, 0.0, 0.0, 
	0.0, 0.0, 1.0, 0.0, 
	0.0, 0.0, 0.0, 1.0};

double rotationAngle = 45.0;
double translationVector[3] = {0.0, 0.0, -10};

camCamera cam;
double viewport[4][4], proj[4][4];

void render(void) {
	pixClearRGB(0.0, 0.0, 0.0);
	depthClearDepths(&buf, 100000000000);
	meshRender(&mesh, &buf, viewport, &sha, unif, tex);
	meshRender(&mesh2, &buf, viewport, &sha2, unif, tex);
}

void handleKeyUp(
        int key, int shiftIsDown, int controlIsDown, int altOptionIsDown, 
        int superCommandIsDown) {
	if (key == GLFW_KEY_ENTER) {
		if (texture.filtering == texLINEAR)
			texSetFiltering(&texture, texNEAREST);
		else
			texSetFiltering(&texture, texLINEAR);
		render();
	}
	if (key == GLFW_KEY_P) {
		if (cam.projectionType == camORTHOGRAPHIC){
			camSetProjectionType(&cam, camPERSPECTIVE);
			camSetFrustum(&cam, M_PI/6.0, 10.0, 10.0, 512.0, 512.0);
			camGetPerspective(&cam, proj);
			vecCopy(16, (double *)proj, &unif[UNIFPROJ]);
		}
		else {
			camSetProjectionType(&cam, camORTHOGRAPHIC);
			camSetFrustum(&cam, M_PI/6.0, 10.0, 10.0, 512.0, 512.0);
			camGetOrthographic(&cam, proj);
			vecCopy(16, (double *)proj, &unif[UNIFPROJ]);
		}
		render();
	}	
}

void handleTimeStep(double oldTime, double newTime) {
	if (floor(newTime) - floor(oldTime) >= 1.0)
		printf("handleTimeStep: %f frames/sec\n", 1.0 / (newTime - oldTime));
/* 	unif[UNIFR] = sin(newTime);
	unif[UNIFG] = cos(oldTime);
	rotationAngle += (newTime - oldTime);
	double isom[4][4];
	double rotation[3][3];
	double axis[3] = {1.0 / sqrt(3.0), 1.0 / sqrt(3.0), 1.0 / sqrt(3.0)};
	mat33AngleAxisRotation(rotationAngle, axis, rotation);
	mat44Isometry(rotation, translationVector, isom);
	vecCopy(16, (double *)isom, &unif[UNIFMODELING]);
	render(); */
}

int main(void) {
	if (pixInitialize(512, 512, "Three Dimensions") != 0)
		return 1;
	if (texInitializeFile(&texture, "download.jpeg") != 0) {
	    pixFinalize();
		return 2;
	}
		if (texInitializeFile(&texture2, "download2.jpeg") != 0) {
	    texFinalize(&texture);
	    pixFinalize();
		return 3;
	}
	if (mesh3DInitializeBox(&mesh, -0.5, 0.5, -1.0, 1.0, -1.0, 1.0) != 0) {
	    texFinalize(&texture2);
		texFinalize(&texture);
	    pixFinalize();
		return 3;
	}
	if (mesh3DInitializeSphere(&mesh2, 1.5, 50, 100)) {
		meshFinalize(&mesh);
	    texFinalize(&texture2);
 		texFinalize(&texture);
	    pixFinalize();
		return 4;
	}
	if (depthInitialize(&buf, 512, 512) != 0){
		meshFinalize(&mesh2);
		meshFinalize(&mesh);
		texFinalize(&texture2);
	    texFinalize(&texture);
	    pixFinalize();
		return 5;
	}

    texSetFiltering(&texture, texNEAREST);
    texSetLeftRight(&texture, texREPEAT);
    texSetTopBottom(&texture, texREPEAT);
	texSetFiltering(&texture2, texNEAREST);
    texSetLeftRight(&texture2, texREPEAT);
    texSetTopBottom(&texture2, texREPEAT);

    sha.unifDim = 3 + 16 + 16;
    sha.attrDim = 3 + 2 + 3;
    sha.varyDim = 4 + 2;
    sha.shadeVertex = shadeVertex;
    sha.shadeFragment = shadeFragment;
    sha.texNum = 1;

	sha2.unifDim = 3 + 16 + 16;
    sha2.attrDim = 3 + 2 + 3;
    sha2.varyDim = 4 + 2;
    sha2.shadeVertex = shadeVertex;
    sha2.shadeFragment = shadeFragment2;
    sha2.texNum = 1;

	//setting up camera
	camSetProjectionType(&cam, camORTHOGRAPHIC);
	camSetFrustum(&cam, M_PI/6.0, 10.0, 10.0, 512.0, 512.0);
	camGetOrthographic(&cam, proj);

	//viewport
	mat44Viewport(512, 512, viewport);

	double isom[4][4];
	double rotation[3][3];
	double axis[3] = {1.0 / sqrt(3.0), 1.0 / sqrt(3.0), 1.0 / sqrt(3.0)};
	mat33AngleAxisRotation(rotationAngle, axis, rotation);
	mat44Isometry(rotation, translationVector, isom);
	vecCopy(16, (double *)isom, &unif[UNIFMODELING]);
	vecCopy(16, (double *)proj, &unif[UNIFPROJ]);

    render();
    pixSetKeyUpHandler(handleKeyUp);
    pixSetTimeStepHandler(handleTimeStep);
    pixRun();
    meshFinalize(&mesh2);
	meshFinalize(&mesh);
	texFinalize(&texture2);
    texFinalize(&texture);
	depthFinalize(&buf);
    pixFinalize();
    return 0;
}


