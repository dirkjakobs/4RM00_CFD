#define sqr(a)      (a*a)
#define max2(a,b)   (a > b ? a : b)
#define max3(a,b,c) (a > b ? max2(a,c) : max2(b,c))
#define min2(a,b)   (a < b ? a : b)
#define min3(a,b,c) (a < b ? min2(a,c) : min2(b,c))

extern void readInput (char *name);

extern void grid(void);
extern void gridbound(void);
extern void init(void);

extern void bound(void);
extern void globcont(void);
extern void derivatives(void);

extern void ucoeff(double **aE, double **aW, double **aN, double **aS, double **aP, double **b);
extern void vcoeff(double **aE, double **aW, double **aN, double **aS, double **aP, double **b);
extern void pccoeff(double **aE, double **aW, double **aN, double **aS, double **aP, double **b);
extern void Tcoeff(double **aE, double **aW, double **aN, double **aS, double **aP, double **b);
extern void epscoeff(double **aE, double **aW, double **aN, double **aS, double **aP, double **b);
extern void kcoeff(double **aE, double **aW, double **aN, double **aS, double **aP, double **b);

extern void conv(void);
extern void solve   (double **fi, double **b, double **aE, double **aW, double **aN, double **aS, double **aP);
extern void solveGS (double **fi, double **b, double **aE, double **aW, double **aN, double **aS, double **aP);
extern void solveSOR_ADI(double **fi, double **b, double **aE, double **aW, double **aN, double **aS, double **aP);
extern void solveSOR(double **fi, double **b, double **aE, double **aW, double **aN, double **aS, double **aP);
extern void solveSIP(double **fi, double **b, double **aE, double **aW, double **aN, double **aS, double **aP);
extern void velcorr(void);
extern void storeresults(void);

extern void calculateuplus(void);
extern void viscosity(void);

extern void output(void);
extern void printConv(double time, int iter);

extern void memalloc(void);
extern int *int_1D_array(int np);
extern double *double_1D_array(int np);
extern double **double_2D_matrix(int nm, int np);
extern void free_double_2D_matrix (double **m, int nm);

// SELF ADDED FUNCTIONS
extern void calc_uplus(void);
extern void animation(void);
