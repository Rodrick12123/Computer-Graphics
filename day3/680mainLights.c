
//Jeremiah  Vic

/* On macOS, compile with...
    clang 680mainLights.c 040pixel.o -lglfw -framework OpenGL -framework Cocoa -framework IOKit
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
#include "680light.c"
#define SCREENWIDTH 512
#define SCREENHEIGHT 512



/*** SPHERES ******************************************************************/





/*** ARTWORK ******************************************************************/



camCamera camera;
double cameraTarget[3] = {0.0, 0.0, 0.0};
double cameraRho = 10.0, cameraPhi = M_PI / 3.0, cameraTheta = M_PI / 3.0;

/* Four spheres. */
#define BODYNUM 4
#define LIGHTNUM 0

bodyBody bodyArray[BODYNUM]; 
double cAmbient[3] = {1.0/4.0, 1.0/4.0, 1.0/4.0};


texTexture texture;
texTexture texture2;
texTexture texture3;
texTexture texture4;
const texTexture *textures[4] = {&texture, &texture2, &texture3, &texture4};
const texTexture **tex = textures;



lightLight lights[LIGHTNUM]; 
/* Based on the uniforms, textures, rayIntersection, and texture coordinates, 
outputs a material. */
void getMaterial(
        int unifDim, const double unif[], int texNum, const texTexture *tex[], 
        const rayIntersection *inter, const double texCoords[2], 
        rayMaterial *material){
            material->hasDiffuse = 1;
            material->hasSpecular = 1;
            material->hasTransmission = 0;
            material->hasMirror = 0;
            material->hasAmbient = 1;
            double sampleTex[tex[0]->texelDim];
            texSample(tex[0], texCoords[0], texCoords[1], sampleTex);
            vecCopy(3, sampleTex, material->cDiffuse);
            //setting specular reflection
            vecCopy(3, unif, material->cSpecular);
            vecCopy(1, &unif[3], &material->shininess);
        }
void getDirectionalLighting(int unifDim, const double unif[], double distance, lightLighting lighting){
    distance = rayINFINITY;
    //calculating uLight
    double d[3] = {0.0,0.0,1.0};
    isoRotateDirection(&lights[0].isometry, d, lighting.uLight);
    //calculating cLight
    vecCopy(3, unif, lighting.cLight);
}
int initializeArtwork(void) {
    camSetProjectionType(&camera, camPERSPECTIVE);
    camSetFrustum(
        &camera, M_PI / 6.0, cameraRho, 10.0, SCREENWIDTH, SCREENHEIGHT);
    camLookAt(&camera, cameraTarget, cameraRho, cameraPhi, cameraTheta);
    double rot[3][3] = {{1.0, 0.0, 0.0}, {0.0, 1.0, 0.0}, {0.0, 0.0, 1.0}};
    for (int k = 0; k < BODYNUM; k += 1)
        isoSetRotation(&(bodyArray[k].isometry), rot);

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
    //material uniforms
    double matUnif[3] = {1.0,1.0,1.0,64.0};

    bodyInitialize(&bodyArray[0], sphUNIFDIM, 4, 1, sphGetIntersection, 
    sphGetTexCoordsAndNormal, getMaterial);
    
    bodySetTexture(&bodyArray[0], 0, &texture);
    //set radious from radii
    double data[1] = {1.0};
    bodySetGeometryUniforms(&bodyArray[0], 0, data, 1);
    bodySetMaterialUniforms(&bodyArray[0], 0, matUnif, 4);

    bodyInitialize(&bodyArray[1], sphUNIFDIM, 4, 1, sphGetIntersection, 
    sphGetTexCoordsAndNormal, getMaterial);

    bodySetTexture(&bodyArray[1], 0, &texture2);
    //set radious from radii
    double data2[1] = {0.5};
    bodySetGeometryUniforms(&bodyArray[1], 0, data2, 1);
    bodySetMaterialUniforms(&bodyArray[1], 0, matUnif, 4);

    bodyInitialize(&bodyArray[2], sphUNIFDIM, 4, 1, sphGetIntersection, 
    sphGetTexCoordsAndNormal, getMaterial);

    bodySetTexture(&bodyArray[2], 0, &texture3);
    //set radious from radii
    double data3[1] = {0.5};
    bodySetGeometryUniforms(&bodyArray[2], 0, data3, 1);
    bodySetMaterialUniforms(&bodyArray[2], 0, matUnif, 4);

    bodyInitialize(&bodyArray[3], sphUNIFDIM, 4, 1, sphGetIntersection, 
    sphGetTexCoordsAndNormal, getMaterial);

    bodySetTexture(&bodyArray[3], 0, &texture4);
    //set radious from radii
    double data4[1] = {0.5};
    bodySetGeometryUniforms(&bodyArray[3], 0, data4, 1);
    bodySetMaterialUniforms(&bodyArray[3], 0, matUnif, 4);


    double transl[3] = {0.0, 0.0, 0.0};
    isoSetTranslation(&bodyArray[0].isometry, transl);
    vec3Set(1.0, 0.0, 0.0, transl);
    isoSetTranslation(&bodyArray[1].isometry, transl);
    vec3Set(0.0, 1.0, 0.0, transl);
    isoSetTranslation(&bodyArray[2].isometry, transl);
    vec3Set(0.0, 0.0, 1.0, transl);
    isoSetTranslation(&bodyArray[3].isometry, transl);

    //initalizing lights
    int LightunifDim = 3;
    lightInitialize(&lights[0], LightunifDim, lights[0].getLighting);
    double lightUnif[3] = {};
    lightSetUniforms(&lights[0], 0, lightUnif, 3);
    getDirectionalLighting();
    isoSetRotation(&lights[0].isometry, lights[0].isometry.rotation);

    return 0;
}

void finalizeArtwork(void) {
    bodyFinalize(&bodyArray[0]);
    bodyFinalize(&bodyArray[1]);
    bodyFinalize(&bodyArray[2]);
    bodyFinalize(&bodyArray[3]);
    texFinalize(&texture);
    texFinalize(&texture2);
    texFinalize(&texture3);
    texFinalize(&texture4);
    lightFinalize(&lights[0]);
    return;
}

/*** RENDERING ****************************************************************/


/*** RENDERING ****************************************************************/

/* Given a ray x(t) = p + t d. Finds the color where that ray hits the scene (or 
the background) and loads the color into the rgb parameter. */
void getSceneColor(
        int bodyNum, const bodyBody bodies[], const double cAmbient[3], 
        int lightNum, const lightLight lights[], const double p[3], 
        const double d[3], double rgb[3]){
    rayMaterial material;
    double bound = rayINFINITY;
    rayIntersection bestInter;
    lightLighting lighting;
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
        //if body has ambient start rgb with ambient term otherwise set body to black
        if (material.hasAmbient == 1){
             vecModulate(3, cAmbient, material.cDiffuse,  rgb);
        }
        else{
            vec3Set(0.0,0.0,0.0,rgb);
        }
        if (material.hasDiffuse == 1 || material.hasSpecular == 1)
        {
            //compute x
            double tTimesd[3];
            double x[3];
            vecScale(3, bestInter.t, d, tTimesd);
            vecAdd(3, p, tTimesd, x);

            //loop through all the lights
            for (int i = 0; i < LIGHTNUM; i++)
            {
                //ask each light for its lighting
                lightGetLighting(&lights[i], x, &lighting);

                //Compute diffuse intensity for specular or diffuse
                double iDiff = vecDot(3,normal,lighting.uLight);

                //if it has specular compute specular and add it onto rgb
                if (material.hasSpecular == 1 && iDiff > 0)
                {
                    double unit[1];
                    double ucamera[3];
                    //vecScale d by -1 to get -d
                    double negativeD[3];
                    vecScale(3, -1, d, negativeD);
                    vecUnit(3, negativeD, ucamera);

                    double modnormlight[3];
                    double urefl[3];
                    double subUlight[3];

                    vecModulate(3, normal, lighting.uLight, modnormlight);
                    vecSubtract(3, normal, lighting.uLight, subUlight);
                    vecScale(3, 2, modnormlight, modnormlight);
                    vecModulate(3, modnormlight, subUlight, urefl);

                    double ispec = vecDot(3, urefl, ucamera);
                    ispec = pow(ispec, material.shininess);

                    //if ispec is less than 0 set it to 0 else if i diff is less <= 0 set i spec to 0
                    if (ispec < 0){
                        ispec = 0;
                    }
                    else if(iDiff <= 0){
                        ispec = 0;
                    }
                    
                    double ispecTimescLight[3];
                    vecScale(3, ispec, lighting.cLight, ispecTimescLight);
                    vec3Set(1.0,1.0,1.0, material.cSpecular);
                    vecModulate(3, ispecTimescLight, material.cSpecular, rgb);
                }
                //if it has diffuse compute and it onto rgb
                if (material.hasDiffuse == 1)
                {
                    double iDiffTimescLight[3];
                    vecScale(3, iDiff, lighting.cLight, iDiffTimescLight);
                    vecModulate(3, iDiffTimescLight, material.cDiffuse,  rgb);
                }
                
                
            }
        }
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
            getSceneColor(BODYNUM, bodyArray, cAmbient, 0, lights, p, d, rgb);
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
    
    finalizeArtwork();
    pixFinalize();
    return 0;
}


