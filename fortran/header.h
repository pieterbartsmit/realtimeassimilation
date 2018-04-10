//
// 
// Calculation Routines
void tracing( double *outpar,double *out,double *inpar,int *nRay, int *maxout, int *outstep, double *direction   );
void calc( double *matrix, double *xp,double *yp,int *np, double *freq, int *nfreq,double *angles,int *nang,int *nsubrays);

// Setup functions
void setnum( double *frac   , int *maxstep );
void setdom( double *inxlim , double *inylim );
void setbat( double *di , double *xi, double *yi, int *nxi, int *nyi );
void setcont( double *d);
void setbndlist( int *inbndlist , int *n);
