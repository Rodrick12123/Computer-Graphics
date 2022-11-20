


/* A resh is a ray-tracing mesh. It has no geometry uniforms outside the 
attached meshMesh. */
#define reshUNIFDIM 0

/* Given vectors a, b - a, and c - a describing a triangle, with the first three 
entries being XYZ. Given point x, with the first three entries being XYZ, such 
that (up to numerical precision) x - a = p (b - a) + q (c - a). Computes p and 
q. Returns 1 on success or 0 on failure. */
int reshGetPQ(
        const double a[], const double bMinusA[], const double cMinusA[], 
        const double x[], double pq[2]) {
    /* For the 3x2 matrix A with columns b - a, c - a, compute A^T A. */
    double aTA[2][2];
    aTA[0][0] = vecDot(3, bMinusA, bMinusA);
    aTA[0][1] = vecDot(3, bMinusA, cMinusA);
    aTA[1][0] = aTA[0][1];
    aTA[1][1] = vecDot(3, cMinusA, cMinusA);
    /* Compute the 2x2 matrix (A^T A)^-1 if possible. */
    double aTAInv[2][2];
    if (mat22Invert(aTA, aTAInv) == 0.0)
        return 0;
    /* Compute the 2x3 matrix (A^T A)^-1 A^T. */
    double aTAInvAT[2][3];
    for (int i = 0; i < 2; i += 1)
        for (int j = 0; j < 3; j += 1)
            aTAInvAT[i][j] = 
                aTAInv[i][0] * bMinusA[j] + aTAInv[i][1] * cMinusA[j];
    /* Then pq = (A^T A)^-1 A^T (x - a). */
    double xMinusA[3];
    vecSubtract(3, x, a, xMinusA);
    pq[0] = vecDot(3, aTAInvAT[0], xMinusA);
    pq[1] = vecDot(3, aTAInvAT[1], xMinusA);
    return 1;
}

/* An implementation of getIntersection for bodies that are reshes. Assumes that 
the data parameter points to an underlying meshMesh with attribute structure 
XYZSTNOP. */
void reshGetIntersection(
        int unifDim, const double unif[], const void *data, 
        const isoIsometry *isom, const double p[3], const double d[3], 
        double bound, rayIntersection* inter) {
    meshMesh *mesh = (meshMesh *)data;
    double locP[3];
    double locD[3];
    rayIntersection bestInter;
    bestInter.t = rayNONE;
    // transform p and d to local coords
    isoUntransformPoint(isom, p, locP);
    isoUnrotateDirection(isom, d, locD);
    //intersect with each triangle keeping track of which triangle is winning
    for (int i = 0; i < mesh->triNum; i++)
    {
        int *triangle = meshGetTrianglePointer(mesh, i);
		double *a = meshGetVertexPointer(mesh, triangle[0]);
		double *b = meshGetVertexPointer(mesh, triangle[1]);
		double *c = meshGetVertexPointer(mesh, triangle[2]);

    

        //compute n
        double n[3];
        double bMINa[3];
        double cMINa[3];
        vecSubtract(3, b, a, bMINa);
        vecSubtract(3, c, a, cMINa);
        vec3Cross(bMINa, cMINa, n);
        
        //compute t
        double aMINp[3];
        vecSubtract(3, a, locP, aMINp);
        double nTIMESsubP = vecDot(3, n, aMINp);
        double nTimesd = vecDot(3, n, locD);

        if(nTimesd != 0){
            //compute ratio
            double t = nTIMESsubP / nTimesd;
            // compute x
            double x[3];
            double tTimesd[3];
            double xMinusc[3];

            vecScale(3, t, locD, tTimesd);
            vecAdd(3, locP, tTimesd, x);
            // get p and q
            double pq[2];
            reshGetPQ(a, bMINa, cMINa, x, pq);

            if (pq[0] >= 0 && pq[1] >= 0 && pq[0] + pq[1] <= 1)
            {
                bound = t;
                bestInter.t = t;
                bestInter.index = i;
            }
            //return rayIntersection
            inter->t = bestInter.t;
            inter->index = bestInter.index;

        }

    }


                
    
}

/* An implementation of getTexCoordsAndNormal for bodies that are reshes. 
Assumes that the data parameter points to an underlying meshMesh with attribute 
structure XYZSTNOP. */
void reshGetTexCoordsAndNormal(
        int unifDim, const double unif[], const void *data, 
        const isoIsometry *isom, const double p[3], const double d[3], 
        const rayIntersection *inter, double texCoords[2], double normal[3]) {
    meshMesh *mesh = (meshMesh *)data;
    /* REPLACE THESE THREE LINES WITH YOUR CODE. (MINE IS 17 LINES.) */
    // texCoords[0] = 0.5;
    // texCoords[1] = 0.5;

    
    double locP[3];
    double locD[3];
    //transform p and d to local coords
    isoUntransformPoint(isom, p, locP);
    isoUnrotateDirection(isom, d, locD);
    //calculate x(t)
    double x[3];
    double tTimesd[3];
    double xMinusc[3];

    vecScale(3, inter->t, locD, tTimesd);
    vecAdd(3, locP, tTimesd, x);
    //triangle info
    int *triangle = meshGetTrianglePointer(mesh, inter->index);
    
    double *vertexA = meshGetVertexPointer(mesh, triangle[0]);
    double *vertexB = meshGetVertexPointer(mesh, triangle[1]);
    double *vertexC = meshGetVertexPointer(mesh, triangle[2]);

    double bmina[8];
    vecSubtract(8, vertexB, vertexA, bmina);
    double cmina[8];
    vecSubtract(8, vertexC, vertexA, cmina);

    double pq[2];
    reshGetPQ(vertexA, bmina, cmina, x, pq);

    double X[8];
    double multp[8];
    vecScale(8, pq[0], bmina, multp);
    double multq[8];
    vecScale(8, pq[1], cmina, multq);
    double vertAdd[8];
    vecAdd(8, vertexA, multp, vertAdd);
    vecAdd(8, vertAdd, multq, X);

    double normalizedX[3];
    double nop[3] = {X[5], X[6], X[7]};
    vecUnit(3, nop, normalizedX);
    isoRotateDirection(isom, normalizedX, normal);
    texCoords[0] = X[3];
    texCoords[1] = X[4];
}


