


/*** 2 x 2 Matrices ***/

/* Pretty-prints the given matrix, with one line of text per row of matrix. */
void mat22Print(const double m[2][2]) {
    int i, j;
    for (i = 0; i < 2; i += 1) {
        for (j = 0; j < 2; j += 1)
            printf("%f    ", m[i][j]);
        printf("\n");
    }
}

/* Returns the determinant of the matrix m. If the determinant is 0.0, then the 
matrix is not invertible, and mInv is untouched. If the determinant is not 0.0, 
then the matrix is invertible, and its inverse is placed into mInv. The output 
CANNOT safely alias the input. */
double mat22Invert(const double m[2][2], double mInv[2][2]) {
    double d = m[0][0]*m[1][1] - m[0][1]*m[1][0];
    if (d != 0) {
        mInv[0][0]=m[1][1]/d;
        mInv[0][1]=-m[0][1]/d;
        mInv[1][0]=-m[1][0]/d;
        mInv[1][1]=m[0][0]/d;

    }
    return d;
}

/* Multiplies a 2x2 matrix m by a 2-column v, storing the result in mTimesV. 
The output CANNOT safely alias the input. */
void mat221Multiply(const double m[2][2], const double v[2], 
        double mTimesV[2]) {
    mTimesV[0] = m[0][0]*v[0] + m[0][1]*v[1];
    mTimesV[1] = m[1][0]*v[0] + m[1][1]*v[1];
    
}

/* Fills the matrix m from its two columns. The output CANNOT safely alias the 
input. */
void mat22Columns(const double col0[2], const double col1[2], double m[2][2]) {
    m[0][0] = col0[0];
    m[0][1] = col1[0];
    m[1][0] = col0[1];
    m[1][1] = col1[1];
}

/* The theta parameter is an angle in radians. Sets the matrix m to the 
rotation matrix corresponding to counterclockwise rotation of the plane through 
the angle theta. */
void mat22Rotation(double theta, double m[2][2]) {
    m[0][0] = cos(theta);
    m[0][1] = -sin(theta);
    m[1][0] = sin(theta);
    m[1][1] = cos(theta);  
}

/* Multiplies the 3x3 matrix m by the 3x3 matrix n. The output CANNOT safely 
alias the input. */
void mat333Multiply(
        const double m[3][3], const double n[3][3], double mTimesN[3][3]){
            mTimesN[0][0] = m[0][0]*n[0][0] + m[0][1]*n[1][0] + m[0][2]*n[2][0];
            mTimesN[0][1] = m[0][0]*n[0][1] + m[0][1]*n[1][1] + m[0][2]*n[2][1];
            mTimesN[0][2] = m[0][0]*n[0][2] + m[0][1]*n[1][2] + m[0][2]*n[2][2];
            mTimesN[1][0] = m[1][0]*n[0][0] + m[1][1]*n[1][0] + m[1][2]*n[2][0];
            mTimesN[1][1] = m[1][0]*n[0][1] + m[1][1]*n[1][1] + m[1][2]*n[2][1];
            mTimesN[1][2] = m[1][0]*n[0][2] + m[1][1]*n[1][2] + m[1][2]*n[2][2];
            mTimesN[2][0] = m[2][0]*n[0][0] + m[2][1]*n[1][0] + m[2][2]*n[2][0];
            mTimesN[2][1] = m[2][0]*n[0][1] + m[2][1]*n[1][1] + m[2][2]*n[2][1];
            mTimesN[2][2] = m[2][0]*n[0][2] + m[2][1]*n[1][2] + m[2][2]*n[2][2];
        }

/* Multiplies the 3x3 matrix m by the 3x1 matrix v. The output CANNOT safely 
alias the input. */
void mat331Multiply(
        const double m[3][3], const double v[3], double mTimesV[3]){
            mTimesV[0] = m[0][0]*v[0] + m[0][1]*v[1] + m[0][2]*v[2];
            // mTimesN[0][1] = m[0][0]*v[0] + m[0][1]*n[1][1] + m[0][2]*n[2][1];
            // mTimesN[0][2] = m[0][0]*v[0] + m[0][1]*n[1][2] + m[0][2]*n[2][2];
            mTimesV[1] = m[1][0]*v[0] + m[1][1]*v[1] + m[1][2]*v[2];
            // mTimesN[1][1] = m[1][0]*v[0] + m[1][1]*n[1][1] + m[1][2]*n[2][1];
            // mTimesN[1][2] = m[1][0]*v[0] + m[1][1]*n[1][2] + m[1][2]*n[2][2];
            mTimesV[2] = m[2][0]*v[0] + m[2][1]*v[1] + m[2][2]*v[2];
            // mTimesN[2][1] = m[2][0]*v[0] + m[2][1]*n[1][1] + m[2][2]*n[2][1];
            // mTimesN[2][2] = m[2][0]*v[0] + m[2][1]*n[1][2] + m[2][2]*n[2][2];
        }

/* Builds a 3x3 matrix representing 2D rotation and translation in homogeneous 
coordinates. More precisely, the transformation first rotates through the angle 
theta (in radians, counterclockwise), and then translates by the vector t. */
void mat33Isometry(double theta, const double t[2], double isom[3][3]){
        isom[0][0] = cos(theta);
        isom[0][1] = -sin(theta);
        isom[0][2] = t[0];
        isom[1][0] = sin(theta);
        isom[1][1] = cos(theta);
        isom[1][2] = t[1];
        isom[2][0] = 0;
        isom[2][1] = 0;
        isom[2][2] = 1;

}

/* Given a length-1 3D vector axis and an angle theta (in radians), builds the 
rotation matrix for the rotation about that axis through that angle. */
void mat33AngleAxisRotation(
        double theta, const double axis[3], double rot[3][3]){
        //Make U
        double U[3][3], U2[3][3];
        U[0][0] = 0;
        U[0][1] = -axis[2];
        U[0][2] = axis[1];
        U[1][0] = axis[2];
        U[1][1] = 0;
        U[1][2] = -axis[0];
        U[2][0] = -axis[1];
        U[2][1] = axis[0];
        U[2][2] = 0;

        //Make U^2
        mat333Multiply(U, U, U2);

        //Make I
        double I[3][3];
        I[0][0] = 1.0;
        I[0][1] = 0.0;
        I[0][2] = 0.0;
        I[1][0] = 0.0;
        I[1][1] = 1.0;
        I[1][2] = 0.0;
        I[2][0] = 0.0;
        I[2][1] = 0.0;
        I[2][2] = 1.0;

        rot[0][0] = I[0][0] + sin(theta)*U[0][0] + (1-cos(theta))*U2[0][0];
        rot[0][1] = I[0][1] + sin(theta)*U[0][1] + (1-cos(theta))*U2[0][1];
        rot[0][2] = I[0][2] + sin(theta)*U[0][2] + (1-cos(theta))*U2[0][2];
        rot[1][0] = I[1][0] + sin(theta)*U[1][0] + (1-cos(theta))*U2[1][0];
        rot[1][1] = I[1][1] + sin(theta)*U[1][1] + (1-cos(theta))*U2[1][1];
        rot[1][2] = I[1][2] + sin(theta)*U[1][2] + (1-cos(theta))*U2[1][2];
        rot[2][0] = I[2][0] + sin(theta)*U[2][0] + (1-cos(theta))*U2[2][0];
        rot[2][1] = I[2][1] + sin(theta)*U[2][1] + (1-cos(theta))*U2[2][1];
        rot[2][2] = I[2][2] + sin(theta)*U[2][2] + (1-cos(theta))*U2[2][2];

}  

void mat33Transpose(const double m[3][3], double mT[3][3]) {
    for (int i = 0; i < 3; i += 1)
        for (int j = 0; j < 3; j += 1)
            mT[i][j] = m[j][i];
}

/* Given two length-1 3D vectors u, v that are perpendicular to each other. 
Given two length-1 3D vectors a, b that are perpendicular to each other. Builds 
the rotation matrix that rotates u to a and v to b. */
void mat33BasisRotation(
        const double u[3], const double v[3], const double a[3], 
        const double b[3], double rot[3][3]){
        double w[3];
        vec3Cross(u, v, w);
        double rT[3][3];
        double r[3][3] = {{u[0], v[0], w[0]},
                        {u[1], v[1], w[1]},
                        {u[2], v[2], w[2]}};
        mat33Transpose(r, rT);

        double ab[3];
        vec3Cross(a, b, ab);
          
        double s[3][3] = {{a[0], b[0], ab[0]},
                        {a[1], b[1], ab[1]},
                        {a[2], b[2], ab[2]}};


        mat333Multiply(s, rT, rot);
        }

/* Computes the transpose M^T of the given matrix M. The output CANNOT safely 
alias the input. */
void mat44Transpose(const double m[4][4], double mT[4][4]) {
    for (int i = 0; i < 4; i += 1)
        for (int j = 0; j < 4; j += 1)
            mT[i][j] = m[j][i];
}

/* Multiplies m by n, placing the answer in mTimesN. The output CANNOT safely 
alias the input. */
void mat444Multiply(
        const double m[4][4], const double n[4][4], double mTimesN[4][4]){
    mTimesN[0][0] = m[0][0]*n[0][0] + m[0][1]*n[1][0] + m[0][2]*n[2][0] + m[0][3]*n[3][0];
    mTimesN[0][1] = m[0][0]*n[0][1] + m[0][1]*n[1][1] + m[0][2]*n[2][1] + m[0][3]*n[3][1];
    mTimesN[0][2] = m[0][0]*n[0][2] + m[0][1]*n[1][2] + m[0][2]*n[2][2] + m[0][3]*n[3][2];
    mTimesN[0][3] = m[0][0]*n[0][3] + m[0][1]*n[1][3] + m[0][2]*n[2][3] + m[0][3]*n[3][3];
    mTimesN[1][0] = m[1][0]*n[0][0] + m[1][1]*n[1][0] + m[1][2]*n[2][0] + m[1][3]*n[3][0];
    mTimesN[1][1] = m[1][0]*n[0][1] + m[1][1]*n[1][1] + m[1][2]*n[2][1] + m[1][3]*n[3][1];
    mTimesN[1][2] = m[1][0]*n[0][2] + m[1][1]*n[1][2] + m[1][2]*n[2][2] + m[1][3]*n[3][2];
    mTimesN[1][3] = m[1][0]*n[0][3] + m[1][1]*n[1][3] + m[1][2]*n[2][3] + m[1][3]*n[3][3];
    mTimesN[2][0] = m[2][0]*n[0][0] + m[2][1]*n[1][0] + m[2][2]*n[2][0] + m[2][3]*n[3][0];
    mTimesN[2][1] = m[2][0]*n[0][1] + m[2][1]*n[1][1] + m[2][2]*n[2][1] + m[2][3]*n[3][1];
    mTimesN[2][2] = m[2][0]*n[0][2] + m[2][1]*n[1][2] + m[2][2]*n[2][2] + m[2][3]*n[3][2];
    mTimesN[2][3] = m[2][0]*n[0][3] + m[2][1]*n[1][3] + m[2][2]*n[2][3] + m[2][3]*n[3][3];
    mTimesN[3][0] = m[3][0]*n[0][0] + m[3][1]*n[1][0] + m[3][2]*n[2][0] + m[3][3]*n[3][0];
    mTimesN[3][1] = m[3][0]*n[0][1] + m[3][1]*n[1][1] + m[3][2]*n[2][1] + m[3][3]*n[3][1];
    mTimesN[3][2] = m[3][0]*n[0][2] + m[3][1]*n[1][2] + m[3][2]*n[2][2] + m[3][3]*n[3][2];
    mTimesN[3][3] = m[3][0]*n[0][3] + m[3][1]*n[1][3] + m[3][2]*n[2][3] + m[3][3]*n[3][3];

}

/* Multiplies m by v, placing the answer in mTimesV. The output CANNOT safely 
alias the input. */
void mat441Multiply(
        const double m[4][4], const double v[4], double mTimesV[4]){
    mTimesV[0] = m[0][0]*v[0] + m[0][1]*v[1] + m[0][2]*v[2] + m[0][3]*v[3];
    mTimesV[1] = m[1][0]*v[0] + m[1][1]*v[1] + m[1][2]*v[2] + m[1][3]*v[3];
    mTimesV[2] = m[2][0]*v[0] + m[2][1]*v[1] + m[2][2]*v[2] + m[2][3]*v[3];
    mTimesV[3] = m[3][0]*v[0] + m[3][1]*v[1] + m[3][2]*v[2] + m[3][3]*v[3];
}

/* Given a rotation and a translation, forms the 4x4 homogeneous matrix 
representing the rotation followed in time by the translation. */
void mat44Isometry(
        const double rot[3][3], const double trans[3], double isom[4][4]){
            isom[0][0] = rot[0][0];
            isom[0][1] = rot[0][1];
            isom[0][2] = rot[0][2];
            isom[0][3] = trans[0];
            isom[1][0] = rot[1][0];
            isom[1][1] = rot[1][1];
            isom[1][2] = rot[1][2];
            isom[1][3] = trans[1];
            isom[2][0] = rot[2][0];
            isom[2][1] = rot[2][1];
            isom[2][2] = rot[2][2];
            isom[2][3] = trans[2];
            isom[3][0] = 0;
            isom[3][1] = 0;
            isom[3][2] = 0;
            isom[3][3] = 1;
        }

/* Sets its argument to the 4x4 zero matrix (which consists entirely of 0s). */
void mat44Zero(double m[4][4]){
    for (int i=0; i<4; i++) {
        for (int j=0; j<4; j++){ 
            m[i][j] = 0;
        }    
    }

}

/* Multiplies the transpose of the 3x3 matrix m by the 3x1 matrix v. To 
clarify, in math notation it computes M^T v. The output CANNOT safely alias the 
input. */
void mat331TransposeMultiply(
        const double m[3][3], const double v[3], double mTTimesV[3]){
           double tr[3][3];
           mat33Transpose(m, tr);
           mat331Multiply(tr, v, mTTimesV);

        }

/* Builds a 4x4 matrix for a viewport with lower left (0, 0) and upper right 
(width, height). This matrix maps a projected viewing volume 
[-1, 1] x [-1, 1] x [-1, 1] to screen [0, w] x [0, h] x [0, 1] (each interval 
in that order). */
void mat44Viewport(double width, double height, double view[4][4]){
    view[0][0] = width/2;
    view[0][1] = 0;
    view[0][2] = 0;
    view[0][3] = width/2;
    view[1][0] = 0;
    view[1][1] = height/2;
    view[1][2] = 0;
    view[1][3] = height/2;
    view[2][0] = 0;
    view[2][1] = 0;
    view[2][2] = 1/2;
    view[2][3] = 1/2;
    view[3][0] = 0;
    view[3][1] = 0;
    view[3][2] = 0;
    view[3][3] = 1;
}

/* Inverse to the matrix produced by mat44Viewport. */
void mat44InverseViewport(double width, double height, double view[4][4]){
    view[0][0] = 2/width;
    view[0][1] = 0;
    view[0][2] = 0;
    view[0][3] = -1;
    view[1][0] = 0;
    view[1][1] = 2/height;
    view[1][2] = 0;
    view[1][3] = -1;
    view[2][0] = 0;
    view[2][1] = 0;
    view[2][2] = 2;
    view[2][3] = -1;
    view[3][0] = 0;
    view[3][1] = 0;
    view[3][2] = 0;
    view[3][3] = 1;
}