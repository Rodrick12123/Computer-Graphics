
//Jeremiah  Vic

/* On macOS, compile with...
    clang 670mainBody.c 040pixel.o -lglfw -framework OpenGL -framework Cocoa -framework IOKit
On Ubuntu, compile with...
    cc 640mainSpheres.c 040pixel.o -lglfw -lGL -lm -ldl
*/
#include <stdio.h>
#include <math.h>
#include <GLFW/glfw3.h>
#include "040pixel.h"

#include "650vector.c"
#include "280matrix.c"
#include "300isometry.c"
#include "300camera.c"
#include "150texture.c"
#include "660ray.c"
#include "670body.c"
#include "670sphere.c"


#define SCREENWIDTH 512
#define SCREENHEIGHT 512



/*** SPHERES ******************************************************************/





/*** ARTWORK ******************************************************************/



camCamera camera;
double cameraTarget[3] = {0.0, 0.0, 0.0};
double cameraRho = 10.0, cameraPhi = M_PI / 3.0, cameraTheta = M_PI / 3.0;

/* Four spheres. */
#define BODYNUM 4
// isoIsometry isoms[BODYNUM];
// double radii[BODYNUM] = {1.0, 0.5, 0.5, 0.5};

bodyBody bodyArray[BODYNUM]; 
double cAmbient[3] = {1.0/4.0, 1.0/4.0, 1.0/4.0};


texTexture texture;
texTexture texture2;
texTexture texture3;
texTexture texture4;
const texTexture *textures[4] = {&texture, &texture2, &texture3, &texture4};
const texTexture **tex = textures;

rayMaterial material;

/* Based on the uniforms, textures, rayIntersection, and texture coordinates, 
outputs a material. */
void getMaterial(
        int unifDim, const double unif[], int texNum, const texTexture *tex[], 
        const rayIntersection *inter, const double texCoords[2], 
        rayMaterial *material){
            material->hasDiffuse = 0;
            material->hasSpecular = 0;
            material->hasTransmission = 0;
            material->hasMirror = 0;
            material->hasAmbient = 1;
            double sampleTex[3];
            texSample(tex[0], texCoords[0], texCoords[1], sampleTex);
            vecCopy(3, sampleTex, material->cDiffuse);
        }

int initializeArtwork(void) {
    camSetProjectionType(&camera, camPERSPECTIVE);
    camSetFrustum(
        &camera, M_PI / 6.0, cameraRho, 10.0, SCREENWIDTH, SCREENHEIGHT);
    camLookAt(&camera, cameraTarget, cameraRho, cameraPhi, cameraTheta);
    double rot[3][3] = {{1.0, 0.0, 0.0}, {0.0, 1.0, 0.0}, {0.0, 0.0, 1.0}};
    for (int k = 0; k < BODYNUM; k += 1)
        isoSetRotation(&(bodyArray[k].isometry), rot);

    // double transl[3] = {0.0, 0.0, 0.0};
    // isoSetTranslation(&(isoms[0]), transl);
    // vec3Set(1.0, 0.0, 0.0, transl);
    // isoSetTranslation(&(isoms[1]), transl);
    // vec3Set(0.0, 1.0, 0.0, transl);
    // isoSetTranslation(&(isoms[2]), transl);
    // vec3Set(0.0, 0.0, 1.0, transl);
    // isoSetTranslation(&(isoms[3]), transl);
    if ( texInitializeFile(&texture, "borealis.jpeg") != 0) {
        pixFinalize();
        return 1;
    }
    if ( texInitializeFile(&texture2, "awesome.png") != 0) {
        pixFinalize();
        return 1;
    }
    if ( texInitializeFile(&texture3, "download.jpeg") != 0) {
        pixFinalize();
        return 1;
    }
    if ( texInitializeFile(&texture4, "download2.jpeg") != 0) {
        pixFinalize();
        return 1;
    }
    
    bodyInitialize(&bodyArray[0], 1, 0, 1, sphGetIntersection, 
    sphGetTexCoordsAndNormal, getMaterial);
    
    bodySetTexture(&bodyArray[0], 0, &texture);
    // bodyArray[0]->textures = tex[0];
    //set radious from radii
    double data[1] = {1.0};
    bodySetGeometryUniforms(&bodyArray[0], 0, data, 1);

    bodyInitialize(&bodyArray[1], 1, 0, 1, sphGetIntersection, 
    sphGetTexCoordsAndNormal, getMaterial);

    bodySetTexture(&bodyArray[1], 0, &texture2);
    // bodyArray[1]->textures = tex[1];
    //set radious from radii
    double data2[1] = {0.5};
    bodySetGeometryUniforms(&bodyArray[1], 0, data2, 1);

    bodyInitialize(&bodyArray[2], 1, 0, 1, sphGetIntersection, 
    sphGetTexCoordsAndNormal, getMaterial);

    bodySetTexture(&bodyArray[2], 0, &texture3);
    // bodyArray[1]->textures = tex[2];
    //set radious from radii
    double data3[1] = {0.5};
    bodySetGeometryUniforms(&bodyArray[2], 0, data3, 1);

    bodyInitialize(&bodyArray[3], 1, 0, 1, sphGetIntersection, 
    sphGetTexCoordsAndNormal, getMaterial);

    bodySetTexture(&bodyArray[3], 0, &texture4);
    // bodyArray[3]->textures = tex[3];
    //set radious from radii
    double data4[1] = {0.5};
    bodySetGeometryUniforms(&bodyArray[3], 0, data4, 1);


    double transl[3] = {0.0, 0.0, 0.0};
    isoSetTranslation(&bodyArray[0].isometry, transl);
    vec3Set(1.0, 0.0, 0.0, transl);
    isoSetTranslation(&bodyArray[1].isometry, transl);
    vec3Set(0.0, 1.0, 0.0, transl);
    isoSetTranslation(&bodyArray[2].isometry, transl);
    vec3Set(0.0, 0.0, 1.0, transl);
    isoSetTranslation(&bodyArray[3].isometry, transl);

    return 0;
}

void finalizeArtwork(void) {
    bodyFinalize(&bodyArray[0]);
    bodyFinalize(&bodyArray[1]);
    bodyFinalize(&bodyArray[2]);
    bodyFinalize(&bodyArray[3]);
    return;
}

/*** RENDERING ****************************************************************/


/*** RENDERING ****************************************************************/

/* Given a ray x(t) = p + t d. Finds the color where that ray hits the scene (or 
the background) and loads the color into the rgb parameter. */
void getSceneColor(
        int bodyNum, const bodyBody bodies[], const double p[3], 
        const double d[3], double rgb[3]){
    double bound = rayINFINITY;
    rayIntersection bestInter;
    int bestI = -1;
    for (int i = 0; i < BODYNUM; i++)
    {
        rayIntersection inter;
        bodyGetIntersection(&bodies[i], p, d, bound, &inter);
        if (inter.t != rayNONE)
        {
            bound = inter.t;
            bestI = i;
            bestInter = inter;
        }
    }
    if (bestI == -1)
    {
        vec3Set(1.0,1.0,1.0,rgb);
    }else
    {
        double texCoor[2];
        double normal[3];
        double sampleTex[tex[0]->texelDim];
        double unif[0];
        bodyGetTexCoordsAndNormal(&bodies[bestI], p, d, &bestInter, texCoor, normal);
        bodyGetMaterial(&bodies[bestI], &bestInter, texCoor, &material);

        
        vecModulate(3, cAmbient, material.cDiffuse,  rgb);
    }
    

    
}

void render(void) {
    /* Build a 4x4 matrix that (along with homogeneous division) takes screen 
    coordinates (x0, x1, 0, 1) to the corresponding world coordinates. */
    double invView[4][4];
    double invProj[4][4];
    double Xscreen[4];
    double invProjXinvView[4][4];
    double camMatrix[4][4];
    double camXprojview[4][4];
    if(camera.projectionType == camPERSPECTIVE){
        camGetInversePerspective(&camera, invProj);
    }
    else if (camera.projectionType == camORTHOGRAPHIC){
        camGetInverseOrthographic(&camera, invProj);
    }
    mat44InverseViewport(SCREENWIDTH, SCREENHEIGHT, invView);
    mat444Multiply(invProj, invView,invProjXinvView);
    isoGetHomogeneous(&camera.isometry, camMatrix);
    mat444Multiply(camMatrix, invProjXinvView,camXprojview);
    /* Declare p and maybe compute d. */
    double p[4], d[3];


    
    
    double camRot[3] = {0.0,0.0,-1.0};
    if (camera.projectionType == camORTHOGRAPHIC)
    {
        isoRotateDirection(&camera.isometry, camRot, d);
    }
    
    

    /* Each screen point is chosen to be on the near plane. */
    double screen[4] = {0.0, 0.0, 0.0, 1.0};
    for (int i = 0; i < SCREENWIDTH; i += 1) {
        screen[0] = i;
        for (int j = 0; j < SCREENHEIGHT; j += 1) {
            screen[1] = j;
            /* Compute p and maybe also d. */
            mat441Multiply(camXprojview,screen,p);
            vecScale(4, 1.0/p[3], p, p);
            double pCam[3];
            vec3Set(camera.isometry.translation[0], camera.isometry.translation[1], camera.isometry.translation[2], pCam);
            if (camera.projectionType == camPERSPECTIVE)
            {
                //pcamera is translation
                vecSubtract(3, p, pCam, d);
            }
            /* Set the pixel to the color of that ray. */
            double rgb[3];
            getSceneColor(BODYNUM, bodyArray, p, d, rgb);
            pixSetRGB(i, j, rgb[0], rgb[1], rgb[2]);
            
            
        }
    }
}



/*** USER INTERFACE ***********************************************************/

void handleKey(
        int key, int shiftIsDown, int controlIsDown, int altOptionIsDown, 
        int superCommandIsDown) {
    if (key == GLFW_KEY_I)
        cameraPhi -= 0.1;
    else if (key == GLFW_KEY_K)
        cameraPhi += 0.1;
    else if (key == GLFW_KEY_J)
        cameraTheta -= 0.1;
    else if (key == GLFW_KEY_L)
        cameraTheta += 0.1;
    else if (key == GLFW_KEY_U)
        cameraRho *= 1.1;
    else if (key == GLFW_KEY_O)
        cameraRho *= 0.9;
    else if (key == GLFW_KEY_P) {
        if (camera.projectionType == camORTHOGRAPHIC)
            camSetProjectionType(&camera, camPERSPECTIVE);
        else
            camSetProjectionType(&camera, camORTHOGRAPHIC);
    }
    camSetFrustum(
        &camera, M_PI / 6.0, cameraRho, 10.0, SCREENWIDTH, SCREENHEIGHT);
    camLookAt(&camera, cameraTarget, cameraRho, cameraPhi, cameraTheta);
}

void handleTimeStep(double oldTime, double newTime) {
    if (floor(newTime) - floor(oldTime) >= 1.0)
        printf(
            "info: handleTimeStep: %f frames/s\n", 1.0 / (newTime - oldTime));
    double rotAxis[3] = {1.0 / sqrt(3.0), 1.0 / sqrt(3.0), 1.0 / sqrt(3.0)};
    double rotMatrix[3][3];
    mat33AngleAxisRotation(newTime, rotAxis, rotMatrix);
    for (int k = 0; k < BODYNUM; k += 1)
        isoSetRotation(&(bodyArray[k].isometry), rotMatrix);
    render();
}

int main(void) {
    if (pixInitialize(SCREENWIDTH, SCREENHEIGHT, "670mainBody") != 0)
        return 1;
    if (initializeArtwork() != 0) {
        pixFinalize();
        return 2;
    }
    pixSetKeyDownHandler(handleKey);
    pixSetKeyRepeatHandler(handleKey);
    pixSetTimeStepHandler(handleTimeStep);
    pixRun();
    texFinalize(&texture);
    texFinalize(&texture2);
    texFinalize(&texture3);
    texFinalize(&texture4);
    finalizeArtwork();
    pixFinalize();
    return 0;
}


