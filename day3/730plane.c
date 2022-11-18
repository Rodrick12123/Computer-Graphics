


/* A plane has no geometry uniforms. */
#define plaUNIFDIM 0

/* An implementation of getIntersection for bodies that are planes. */
void plaGetIntersection(
        int unifDim, const double unif[], const void *data, const isoIsometry *isom, 
        const double p[3], const double d[3], double bound, 
        rayIntersection* inter) {
                double locP[3];
                double locD[3];
                //transform p and d to local coords
                isoUntransformPoint(isom, p, locP);
                isoUnrotateDirection(isom, d, locD);
                if(locD[2] == 0){
                        inter->t = rayNONE;
                        return;
                } 
                double t = -locP[2]/locD[2];
                //if local d[2] equals 0 set t to rayNONE 
                
                //if t is between epsilon and bound return it 
                if (rayEPSILON > t || t > bound){
                        inter->t = rayNONE;
                        return;
                }
                inter->t = t;
                return;
    
}

/* An implementation of getTexCoordsAndNormal for bodies that are planes. */
void plaGetTexCoordsAndNormal(
        int unifDim, const double unif[], const void *data, const isoIsometry *isom, 
        const double p[3], const double d[3], const rayIntersection *inter, 
        double texCoords[2], double normal[3]) {

                double winningT = inter->t;
                double locP[3];
                double locD[3];
                //transform p and d to local coords
                isoUntransformPoint(isom, p, locP);
                isoUnrotateDirection(isom, d, locD);
                //transform normal to global
                double localGlobal[3] = {0,0,1};
                isoRotateDirection(isom, localGlobal, normal);

                //computing tex coords
                texCoords[0] = locP[0] + (winningT * locD[0]);
                texCoords[1] = locP[1] + (winningT * locD[1]);

    
}


