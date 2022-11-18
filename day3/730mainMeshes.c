
//Jeremiah  Vic

/* On macOS, compile with...
    clang 730mainMeshes.c 040pixel.o -lglfw -framework OpenGL -framework Cocoa -framework IOKit
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
#include "730ray.c"
#include "730body.c"
#include "730sphere.c"
#include "680light.c"
#include "730plane.c"
#include "730mesh.c"
#include "250mesh3D.c"
#include "740resh.c"
#define SCREENWIDTH 512
#define SCREENHEIGHT 512



/*** SPHERES ******************************************************************/





/*** ARTWORK ******************************************************************/



camCamera camera;
double cameraTarget[3] = {0.0, 0.0, 0.0};
double cameraRho = 10.0, cameraPhi = M_PI / 3.0, cameraTheta = M_PI / 3.0;

/* Four spheres. */
#define BODYNUM 5
#define LIGHTNUM 2

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
void getPhongMaterial(
        int unifDim, const double unif[], const void *data, int texNum, const texTexture *tex[], 
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

void getMirrorMaterial(
        int unifDim, const double unif[], const void *data, int texNum, const texTexture *tex[], 
        const rayIntersection *inter, const double texCoords[2], 
        rayMaterial *material){

            material->hasDiffuse = 0;
            material->hasSpecular = 0;
            material->hasAmbient = 0;
            material->hasTransmission = 0;

            material->hasMirror = 1;
            double cMirror[3] = {1,1,1};
            vecCopy(3, cMirror, material->cMirror);
        }

void getDirectionalLighting(
        int unifDim, const double unif[], const isoIsometry *isometry, const double x[3],  lightLighting *lighting){
        lighting->distance = rayINFINITY;
        //calculating uLight
        double d[3] = {0.0,0.0,1.0};
        isoRotateDirection(isometry, d, lighting->uLight);
        //calculating cLight
        vecCopy(3, unif, lighting->cLight);
    }

void getPositionalLighting(
        int unifDim, const double unif[], const isoIsometry *isometry, const double x[3],  lightLighting *lighting){
        double subx[3];
        double plight[3];
        vecCopy(3, isometry->translation, plight);
        vecSubtract(3, plight, x, subx);
        //calc d
        lighting->distance = vecLength(3,subx); //figure this out
        //calculating uLight
        double ulight[3];
        vecUnit(3,subx,ulight);
        vecCopy(3, ulight, lighting->uLight);
        //calculating cLight
        vecCopy(3, unif, lighting->cLight);
    }

/* Casts the ray x(t) = p + t d into the scene. Returns 0 if it hits no body or 
1 if it hits any body. Used to determine whether a fragment is in shadow. */
int getSceneShadow(
        int bodyNum, const bodyBody bodies[], const double p[3], 
        const double d[3]){

    double bound = rayINFINITY;
    for (int i = 0; i < bodyNum; i++)
    {
        rayIntersection inter;
        bodyGetIntersection(&bodies[i], p, d, bound, &inter);
        if (inter.t != rayNONE)
        {
            
            return 1;
        }
    }
    return 0;
    
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
    if ( texInitializeFile(&texture2, "download2.jpeg") != 0) {
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
    double matUnif[4] = {1.0,1.0,1.0,64.0};

    bodyInitialize(&bodyArray[0], sphUNIFDIM, 4, 1, sphGetIntersection, 
    sphGetTexCoordsAndNormal, getPhongMaterial);
    
    bodySetTexture(&bodyArray[0], 0, &texture);
    //set radious from radii
    double data[1] = {1.0};
    bodySetGeometryUniforms(&bodyArray[0], 0, data, 1);
    bodySetMaterialUniforms(&bodyArray[0], 0, matUnif, 4);
    //body 1 is mirror material now
    bodyInitialize(&bodyArray[1], sphUNIFDIM, 4, 1, sphGetIntersection, 
    sphGetTexCoordsAndNormal, getPhongMaterial);

    bodySetTexture(&bodyArray[1], 0, &texture2);
    //set radious from radii
    double data2[1] = {0.5};
    bodySetGeometryUniforms(&bodyArray[1], 0, data2, 1);
    bodySetMaterialUniforms(&bodyArray[1], 0, matUnif, 4);
    //body 2 is mirror material now
    bodyInitialize(&bodyArray[2], sphUNIFDIM, 4, 1, sphGetIntersection, 
    sphGetTexCoordsAndNormal, getMirrorMaterial);

    bodySetTexture(&bodyArray[2], 0, &texture3);
    //set radious from radii
    double data3[1] = {0.5};
    bodySetGeometryUniforms(&bodyArray[2], 0, data3, 1);
    bodySetMaterialUniforms(&bodyArray[2], 0, matUnif, 4);

    bodyInitialize(&bodyArray[3], sphUNIFDIM, 4, 1, sphGetIntersection, 
    sphGetTexCoordsAndNormal, getMirrorMaterial);
    bodySetTexture(&bodyArray[3], 0, &texture4);
    //set radious from radii
    double data4[1] = {0.5};
    bodySetGeometryUniforms(&bodyArray[3], 0, data4, 1);
    bodySetMaterialUniforms(&bodyArray[3], 0, matUnif, 4);

    //5th body
    bodyInitialize(&bodyArray[4], plaUNIFDIM, 4, 1, plaGetIntersection, 
    plaGetTexCoordsAndNormal, getPhongMaterial);
    bodySetTexture(&bodyArray[4], 0, &texture4);
    bodySetMaterialUniforms(&bodyArray[4], 0, matUnif, 4);
    


    double transl[3] = {0.0, 0.0, 0.0};
    isoSetTranslation(&bodyArray[0].isometry, transl);
    vec3Set(1.5, 0.0, 0.0, transl);
    isoSetTranslation(&bodyArray[1].isometry, transl);
    vec3Set(0.0, 1.5, 0.0, transl);
    isoSetTranslation(&bodyArray[2].isometry, transl);
    vec3Set(0.0, 0.0, 1.5, transl);
    isoSetTranslation(&bodyArray[3].isometry, transl);
    vec3Set(0.0,0,-1, transl);
    isoSetTranslation(&bodyArray[4].isometry, transl);

    //initalizing lights
    int LightunifDim = 3;
    lightInitialize(&lights[0], LightunifDim, getDirectionalLighting);
    double lightUnif[3] = {1.0,1.0,1.0};
    lightSetUniforms(&lights[0], 0, lightUnif, 3);
    //setting the lights isometry rotation
    double axis[3] = {1/sqrt(2),1/sqrt(2),0};
    double rott[3][3];
    mat33AngleAxisRotation(M_PI/4, axis, rott);
    isoSetRotation(&lights[0].isometry, rott);

    int LightunifDim2 = 3;
    lightInitialize(&lights[1], LightunifDim2, getPositionalLighting);
    double lightUnif2[3] = {1.0,1.0,1.0};
    lightSetUniforms(&lights[1], 0, lightUnif2, 3);
    //setting the lights isometry rotation

    double trans[3] = {1,1,2};
    isoSetTranslation(&lights[1].isometry, trans);
    
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
    lightFinalize(&lights[1]);
    return;
}

/*** RENDERING ****************************************************************/


/*** RENDERING ****************************************************************/

/* Given a ray x(t) = p + t d. Finds the color where that ray hits the scene (or 
the background) and loads the color into the rgb parameter. */
void getSceneColor(
        int recDepth, int bodyNum, const bodyBody bodies[], const double cAmbient[3], 
        int lightNum, const lightLight lights[], const double p[3], 
        const double d[3], double rgb[3]){
    rayMaterial material;
    double bound = rayINFINITY;
    rayIntersection bestInter;
    lightLighting lighting;
    int bestI = -1;
    recDepth = 3;
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
        vec3Set(0.0,0.0,0.0,rgb);
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

                //if it has specular compute specular and add it onto rgb
                if (getSceneShadow(BODYNUM, bodies, x, lighting.uLight) == 0)
                {
                    double iDiff = vecDot(3,normal,lighting.uLight);
                
                
                
                    if (iDiff <= 0)
                    {
                        iDiff = 0;
                    }
                
                    if (material.hasSpecular == 1 && iDiff > 0)
                    {
                        double ucamera[3];
                        //vecScale d by -1 to get -d
                        double negativeD[3];
                        vecScale(3, -1, d, negativeD);
                        vecUnit(3, negativeD, ucamera);

                        double dotnormlight;
                        double urefl[3];
                        double scaleN[3];

                        dotnormlight = 2 * vecDot(3, lighting.uLight, normal);
                        vecScale(3, dotnormlight, normal, scaleN);
                        vecSubtract(3, scaleN, lighting.uLight, urefl);

                        double ispec = vecDot(3, urefl, ucamera);
                        

                        //if ispec is less than 0 set it to 0 else if i diff is less <= 0 set i spec to 0
                        if (ispec < 0){
                            ispec = 0;
                        }
                        
                        ispec = pow(ispec, material.shininess);
                        double rgb2[3];
                        double ispecTimescLight[3];
                        vecScale(3, ispec, lighting.cLight, ispecTimescLight);
                        
                        vecModulate(3, ispecTimescLight, material.cSpecular, rgb2);
                        vecAdd(3, rgb, rgb2,  rgb);
                    }
                    //if it has diffuse compute and it onto rgb
                    if (material.hasDiffuse == 1)
                    {
                        double iDiffTimescLight[3];
                        double rgb2[3];
                        vecScale(3, iDiff, lighting.cLight, iDiffTimescLight);
                        vecModulate(3, iDiffTimescLight, material.cDiffuse,  rgb2);
                        vecAdd(3, rgb, rgb2,  rgb);
                    }
                            
                }

            }
        }
           //if recDepth is 0 turn off mirroring
                if (recDepth == 0)
                {
                    material.hasMirror = 0;
                }
                // checking to see if body has a mirror material
                if (material.hasMirror == 1 && recDepth > 0)
                {
                    // compute x
                    double tTimesd[3];
                    double x[3];
                    vecScale(3, bestInter.t, d, tTimesd);
                    vecAdd(3, p, tTimesd, x);
                    double cFromMirror[3];
                    //computing new d
                    double negativeD[3];
                    double refl[3];
                    double newD[3];
                    vecScale(3, -1, d, negativeD);
                    double negDTimesNormal = 2 * vecDot(3, negativeD, normal);
                    vecScale(3, negDTimesNormal, normal, refl);
                    vecSubtract(3, refl, negativeD, newD);
                    // call getSceneColor to compute cFromMirror
                    getSceneColor(recDepth - 1, BODYNUM, bodies, cAmbient, LIGHTNUM, lights, x, newD, cFromMirror);
                    double rgb2[3];
                    // mirror contribution
                    vecModulate(3, material.cMirror, cFromMirror, rgb2);
                    vecAdd(3, rgb, rgb2, rgb);
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
            //fix for infinite recursion
            int recDepth = 3;
            if(recDepth > 0){
                getSceneColor(recDepth - 1, BODYNUM, bodyArray, cAmbient, 0, lights, p, d, rgb);

            }
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
    for (int k = 0; k < BODYNUM - 1; k += 1)
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


