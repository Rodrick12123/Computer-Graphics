


/* A plane has no geometry uniforms. */
#define plaUNIFDIM 0

/* An implementation of getIntersection for bodies that are planes. */
void plaGetIntersection(
        int unifDim, const double unif[], const isoIsometry *isom, 
        const double p[3], const double d[3], double bound, 
        rayIntersection* inter) {
    
}

/* An implementation of getTexCoordsAndNormal for bodies that are planes. */
void plaGetTexCoordsAndNormal(
        int unifDim, const double unif[], const isoIsometry *isom, 
        const double p[3], const double d[3], const rayIntersection *inter, 
        double texCoords[2], double normal[3]) {
    
}


