

// ### VARIABLES IMPORTED FROM TEXT FILE ###

typedef enum { false, true } bool;
bool ***CONS;

double XMAX;
double YMAX;
int    NPI;
int    NPJ;
double relax_u;
double relax_T;
int    MAX_ITER;
int    U_ITER;
int    V_ITER;
int    PC_ITER;
int    T_ITER;
int    EPS_ITER;
int    K_ITER;
double rho_init;
double mu_init;
double Cp_init;
double K_init;
double Gamma_init;

double TZERO;
double TEMP;
double U_IN;

double *x;               /* x coordinate on pressure points [m] */
double *x_u;             /* x coordinate on u-velocity points [m] */
double *y;               /* y coordinate on pressure points [m] */
double *y_v;             /* y coordinate on v-velocity points [m] */
double Dt;

double **u;
double **v;
double **pc;
double **p;
double **T;
double **rho;
double **mu;
double **mut;
double **mueff;
double **E;
double **E2;
double **Gamma;
double **Cp;
double **k;
double **gammak;
double **eps;
double **gammaeps;
double **delta;
double **Prandtl;
double **uplus;
double **yplus;
double **yplus1;
double **yplus2;
double **tw;

double **u_old;
double **v_old;
double **pc_old;
double **T_old;
double **k_old;
double **eps_old;

double **dudx;
double **dudy;
double **dvdx;
double **dvdy;

double **aE;
double **aW;
double **aN;
double **aS;
double **aP;
double **b;

int    Istart;
int    Iend;
int    Jstart;
int    Jend;

int    iter;
int    iter_u, iter_v, iter_pc, iter_T, iter_k, iter_eps;
double relax_u, relax_v, relax_pc, relax_T;

double SAVG;
double SMAX;

double **SP;
double **SP_u;
double **SP_v;
double **Su;
double **Su_u;
double **Su_v;

double **F_u;
double **F_v;

double **d_u;
double **d_v;

double m_in;
double m_out;

double *relax;

double omega;

// SELF ADDED VARIABLES
double **yplus_u;
double **yplus_v;
double **uplus_u;
double **uplus_v;
double **Tplus_u;
double **Tplus_v;
double **lambda;
double **p_abs;





