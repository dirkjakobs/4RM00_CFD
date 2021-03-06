/***** Solves: Unsteady, compressible convection-diffusion problems.

****** Description:
****** This program solves unsteady convection-diffusion problems	
****** using the transient simple algorithm described in ch. 8.7.1 in "Computational 
****** Fluid Dynamics" by H.K. Versteeg and W. Malalasekera. Symbols and
****** variables follow exactly the notations in this reference, and all 
****** equations cited are from this reference unless mentioned otherwise.

****** References: 1. Computational Fluid Dynamics, H.K. Versteeg and W. 
******			    Malalasekera, Longman Group Ltd, 1995
******/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "variables.h"
#include "constants.h"
#include "functions.h"


/* ################################################################# */
int main(int argc, char *argv[])
/* ################################################################# */
{
	int    iter_u, iter_v, iter_pc, iter_T, iter_eps, iter_k, Count;
	double du, dv, time;
	
	readInput("constraints.dat");

	init();
	bound(); /* apply boundary conditions */
	
	Count = 0;
//	animation(Count);
	
	for (time = Dt; time <= TOTAL_TIME; time += Dt) {
		iter = 0; 
		while (iter < MAX_ITER && SMAX > SMAXneeded && SAVG > SAVGneeded) { /* outer iteration loop */

			derivatives();
			ucoeff(aE, aW, aN, aS, aP, b);
			for (iter_u = 0; iter_u < U_ITER; iter_u++)
				solve(u, b, aE, aW, aN, aS, aP);

			vcoeff(aE, aW, aN, aS, aP, b);
			for (iter_v = 0; iter_v < V_ITER; iter_v++)
				solve(v, b, aE, aW, aN, aS, aP);

			bound();
			pccoeff(aE, aW, aN, aS, aP, b);
			for (iter_pc = 0; iter_pc < PC_ITER; iter_pc++)
				solve(pc, b, aE, aW, aN, aS, aP);

			velcorr(); /* Correct pressure and velocity */

			kcoeff(aE, aW, aN, aS, aP, b);
			for (iter_k = 0; iter_k < K_ITER; iter_k++)
				solve(k, b, aE, aW, aN, aS, aP);

			epscoeff(aE, aW, aN, aS, aP, b);
			for (iter_eps = 0; iter_eps < EPS_ITER; iter_eps++)
				solve(eps, b, aE, aW, aN, aS, aP);

			Tcoeff(aE, aW, aN, aS, aP, b);
			for (iter_T = 0; iter_T < T_ITER; iter_T++)
				solve(T, b, aE, aW, aN, aS, aP);

			properties(); 
						
			bound();
			storeresults(); /* Store data at current time level in arrays for "old" data*/

			iter++;
		} /* for outer iteration loop */

		printConv(time,iter); /* print convergence to the screen */

		/* reset SMAX and SAVG */
		SMAX = LARGE;
		SAVG = LARGE;
	
		Count++;
//		animation(Count);
		
	} /* for Dt */
	output();

	return 0;

} /* main */


/* ################################################################# */
void grid(void)
/* ################################################################# */
{
/***** Purpose: Defining the geometrical variables ******/
/*****          See fig. 6.2-6.4 in ref. 1 ******/
	int    I, J, i, j;
	double Dx, Dy;

	/* Length of volume element */

	Dx = XMAX/NPI;
	Dy = YMAX/NPJ;

	/* Length variable for the scalar points in the x direction */

	x[0] = 0.;
	x[1] = 0.5*Dx;

	for (I = 2; I <= NPI; I++)
		x[I] = x[I-1] + Dx;

	x[NPI+1] = x[NPI] + 0.5*Dx;

	/* Length variable for the scalar points fi[i][j] in the y direction */

	y[0] = 0.;
	y[1] = 0.5*Dy;

	for (J = 2; J <= NPJ; J++)
		y[J] = y[J-1] + Dy;

	y[NPJ+1] = y[NPJ] + 0.5*Dy;

	/* Length variable for the velocity components u[i][j] in the x direction */

	x_u[0] = 0.;
	x_u[1] = 0.;

	for (i = 2; i <= NPI + 1; i++)
		x_u[i] = x_u[i-1] + Dx;

	/* Length variable for the velocity components v[i][j] in the y direction */

	y_v[0] = 0.;
	y_v[1] = 0.;
	for (j = 2; j <= NPJ + 1; j++)
		y_v[j] = y_v[j-1] + Dy;

} /* grid */

/* ################################################################# */
void init(void)
/* ################################################################# */
{
/***** Purpose: To initilise all parameters. ******/
	int    I, J, i, j;
	double Dx, Dy;
	
	memalloc();
	grid();
	
	// Initialize for animation
	timestore = 1;//1/Dt;
	nloop = 0;
	
	/* Initialising all variables  */

	omega = 1.0; /* Over-relaxation factor for SOR solver ( 1 < omega < 2 ) */

	/* Initialize convergence parameters at large values */

	SMAX = LARGE;
	SAVG = LARGE;
		
	// PIM: make Dx and Dy global!
	Dx = XMAX/NPI;
	Dy = YMAX/NPJ;

	m_in  = 1.;
	m_out = 1.; 

	for (I = 0; I <= NPI + 1; I++) {
		i = I;
		for (J = 0; J <= NPJ + 1; J++) {
			j = J;
			u      [i][J] = U_IN;						/* Guess velocity profile in x-direction */
			v      [I][j] = 0.;       					/* Velocity in y-direction */
			p      [I][J] = 0.;       					/* Relative pressure */
			T      [I][J] = TZERO;   					/* Temperature, obtained from text file*/ 
			k      [I][J] = 1e-3;     					/* k */
			eps    [I][J] = 1e-4;     					/* epsilon */
			tw     [I][J] = 5.;                         /* tw */
			rho    [I][J] = rho_init;      				/* Density */
			mu     [I][J] = mu_init;    				/* Viscosity */
			Cp     [I][J] = Cp_init;     				/* J/(K*kg) Heat capacity - assumed constant for this problem */
			lambda [I][J] = lambda_init;				/* Thermal conductivity k [W/(mK)] */
			Gamma  [I][J] = lambda[I][J]/Cp[I][J];      /* Thermal conductivity divided by heat capacity - assumed laminair in initial state */
			Prandtl[I][J] = mu[I][J]/Gamma[I][J]; 		/* laminar (or moleculair) Prandtl number mu*Cp/K (eq. 3.50) */
			u_old  [i][J] = u[i][J];  					/* Velocity in x-direction old timestep */
			v_old  [I][j] = v[I][j];  					/* Velocity in y-direction old timestep */
			pc_old [I][J] = pc[I][J]; 					/* Pressure correction old timestep */
			T_old  [I][J] = T[I][J];  					/* Temperature old timestep */
			eps_old[I][J] = eps[I][J];  				/* epsilon old timestep*/
			k_old  [I][J] = k[I][J];    				/* k old timestep*/
			
			// Set yplus and uplus to 1
			yplus [I][J] = 1.;
			xplus [I][J] = 1.;
			uplus [I][J] = 1.;
			vplus [I][J] = 1.;
			
			// Guess yplus near object/wall		
	        if (CONS[I][J][YPLUS] == true) {
				yplus [I][J] = sqrt(rho[I][J] * u[I][J] / mu[I][J]) * (0.5*Dy);
	        } /* if */
	        
	        // Guess xplus near object/wall	
	        if (CONS[I][J][XPLUS] == true) {
				xplus [I][J] = sqrt(rho[I][J] * v[I][J] / mu[I][J]) * (0.5*Dx);
	        } /* if */		
			
		} /* for J */
	} /* for I */

	/* Setting the relaxation parameters */

	/* relax_u is set by ReadInput, See eq. 6.36 */
	/* relax_T is set by ReadInput, Relaxation factor for temperature */
	relax_v   = relax_u;       /* See eq. 6.37 */
	relax_pc  = 1.1 - relax_u; /* See eq. 6.33 */
	
} /* init */

/* ################################################################# */
void bound(void)
/* ################################################################# */
{
/***** Purpose: Specify boundary conditions for a calculation ******/
	int    I, J, i, j;

	/* Fixed temperature at the upper and lower wall */

	for (J = 0; J <= NPJ + 1; J++) {
		/* Temperature at the walls in Kelvin */
		u[1][J] = U_IN; /* inlet */
	} /* for J */

	for (I = 0; I <= NPI + 1; I++) {
		/* Temperature at the walls in Kelvin */
		T[I][0] = TZERO; /* bottom wall */
		T[I][NPJ+1] = TZERO; /* top wall */
	} /* for J */

	globcont();

	/* Velocity and temperature gradient at outlet = zero: */

	for (J = 1; J <= NPJ; J++) {
		j = J;
		/* Correction factor m_in/m_out is used to satisfy global continuity */
		u[NPI+1][J] = u[NPI][J]*m_in/m_out;
		v[NPI+1][j] = v[NPI][j];
		k[NPI+1][J] = k[NPI][J];
		eps[NPI+1][J] = eps[NPI][J];
	} /* for J */

	for (J = 0; J <= NPJ+1; J++) {
		// Set temperature gradient to zero
		T[NPI+1][J] = T[NPI][J];
		
		// Set relative pressure outlet to zero
		p[NPI+1][J] = 0;
	} /* for J */

	for (J = 0; J <= NPJ + 1; J++) {
		k[0][J] = 2./3.*sqr(U_IN*Ti); /* inlet */
		eps[0][J] = pow(Cmu,0.75)*pow(k[0][J],1.5)/(0.07*YMAX*0.5); /* inlet */
	} /* for J */
		
	/* Set Wall Boundary values if applicable
	   Boundary conditions need to be set if there is NO wall*/
	
	for (I = 1; I <= NPI + 1; I++) {
		i = I;
		// Lower wall gradient boundary conditions
		if (CONS[I][0][0] != true) {
			u[i][0] = u[i][1];
			v[I][0] = 0;				// velocity at walls in y-direction should be zero.
			k[I][0] = k[I][1];
			eps[I][0] = eps[I][1];
			T[I][0] = T[I][1];
		}
		// Upper wall gradient boundary conditionstrue
		if (CONS[I][NPJ+1][0] != true) {
			u[i][NPJ+1] = u[i][NPJ];
			v[I][NPJ+1] = 0;			// velocity at walls in y-direction should be zero.
			k[I][NPJ+1] = k[I][NPJ];
			eps[I][NPJ+1] = eps[I][NPJ];
			T[I][NPJ+1] = T[I][NPJ];
		}

	} /* for I */

} /* bound */

/* ################################################################# */
void globcont(void)
/* ################################################################# */
{
/***** Purpose: Calculate mass in and out of the calculation domain to ******/
/*****          correct for the continuity at outlet. ******/
	int    J, j;
	double AREAw;

	conv();

	m_in = 0.;
	m_out = 0.;

	for (J = 1; J <= NPJ; J++) {
		j = J;
		AREAw = y_v[j+1] - y_v[j]; /* See fig. 6.3 */
		m_in  += F_u[1][J]*AREAw;
		m_out += F_u[NPI][J]*AREAw;
	} /* for J */

} /* globcont */

/* ################################################################# */
void derivatives(void)
/* ################################################################# */
{
/***** Purpose: To calculate derivatives ******/
	int    I, J, i, j;
	
	for (I = 1; I <= NPI; I++) {
		i = I;
		for (J = 1; J <= NPJ; J++) {
			j = J;
			dudx[I][J] = (u[i+1][J] - u[i][J])   / (x_u[i+1] - x_u[i]);
			dudy[i][j] = (u[i][J]   - u[i][J-1]) / (y[J]     - y[J-1]);
			dvdx[i][j] = (v[I][j]   - v[I-1][j]) / (x[I]     - x[I-1]);
			dvdy[I][J] = (v[I][j+1] - v[I][j])   / (y_v[j+1] - y_v[j]);
			E [i][j] = sqr(dudy[i][j]) + sqr(dvdx[i][j]) + 2 * dudy[i][j] * dvdx[i][j];
         } /* for J */
   }  /* for I */                 

	for (I = 1; I <= NPI; I++) {
		i = I;
		for (J = 1; J <= NPJ; J++) {
			j = J;
			E2[I][J] = sqr(dudx[I][J]) + sqr(dvdy[I][J]) + 0.25*(E[i][j] + E[i+1][j] + E[i][j+1] + E[i+1][j+1]);
         } /* for J */
   }  /* for I */                 

} /* derivatives */

/* ################################################################# */
void solve(double **fi, double **b, double **aE, double **aW, double **aN, double **aS, double **aP)
/* ################################################################# */
{
/***** Purpose: To solve the algebraic equation 7.7. ******/
	int    I, J, space;
	double *Ari, *Cmri, Cri;

/* TDMA along a horizontal row from west to east. equation to solve: */

	/* - aW*fiW + aP*fiP - aE*fiE = aS*fiS + aN*fiN + b */

	/* equivalences with variables in eq. 7.1-7.6: */
	/* BETA = aW[I][J] Def. in eq. 7.2 */
	/* D    = aP[I][J] Def. in eq. 7.2 */
	/* ALFA = aE[I][J] Def. in eq. 7.2 */
	/* A    = Ari[I]	 Def. in eq. 7.6b */
	/* C    = Cri	 The right side assumed temporarily known (see eq. 7.8) */
	/* C�   = Cmri[I]  Def. in eq. 7.6c */
	/* b    = b[I][J]	 Def. in eq. 7.7 */

	space = max2((Iend - Istart + 3),(Jend - Jstart + 3));
	Ari   = double_1D_array(space);
	Cmri  = double_1D_array(space);

	/* Solving the (e-w) lines from the south */

	for (J = Jstart; J <= Jend; J++) {
		/* At the inlet boundary: */
		Ari [Istart-1] = 0;
		Cmri[Istart-1] = fi[Istart-1][J];

		for (I = Istart; I <= Iend; I++) { /* Forward substitution */
			Ari[I]  = aE[I][J]/(aP[I][J] - aW[I][J]*Ari[I-1]); /* eq. 7.6b */
			Cri     = aN[I][J]*fi[I][J+1] + aS[I][J]*fi[I][J-1] + b[I][J];
			Cmri[I] = (aW[I][J]*Cmri[I-1] + Cri)/(aP[I][J] - aW[I][J]*Ari[I-1]); /* eq. 7.6c */
		}

		for (I = Iend; I >= Istart; I--)  /* Back substitution */
			fi[I][J] = Ari[I]*fi[I+1][J] + Cmri[I]; /* eq. 7.6a */
	} /* for J from south */

	/* Solving the (e-w) lines from the north */

	for (J = Jend-1; J >= Jstart; J--) {
		/* At the inlet boundary: */
		Ari [Istart-1] = 0;
		Cmri[Istart-1] = fi[Istart-1][J];

		for (I = Istart; I <= Iend; I++) { /* Forward substitution */
			Ari[I]  = aE[I][J]/(aP[I][J] - aW[I][J]*Ari[I-1]); /* eq. 7.6b */
			Cri     = aN[I][J]*fi[I][J+1] + aS[I][J]*fi[I][J-1] + b[I][J];
			Cmri[I] = (aW[I][J]*Cmri[I-1] + Cri)/(aP[I][J] - aW[I][J]*Ari[I-1]);  /* eq. 7.6c */
		}

		for (I = Iend; I >= Istart; I--)  /* Back substitution */
			fi[I][J] = Ari[I]*fi[I+1][J] + Cmri[I]; /* eq. 7.6a */
	} /* for J from north */

/* TDMA along a vertical column from south to north. equation to solve: */

	/* - aS*fiW + aP*fiP - aN*fiE = aW*fiS + aE*fiN + b (eq. 7.8) */

	/* equivalences with variables in eq. 7.1-7.6: */
	/* BETA = aS[I][J] Def. in eq. 7.2 */
	/* D    = aP[I][J] Def. in eq. 7.2 */
	/* ALFA = aN[I][J] Def. in eq. 7.2 */
	/* A    = Ari[I]	 Def. in eq. 7.6b */
	/* C    = Cri      The right side assumed temporarily known (see eq. 7.8) */
	/* C�   = Cmri[I]  Def. in eq. 7.6c */
	/* b    = b[I][J]	 Def. in eq. 7.7 */

	/* Solving (n-s) lines from the west */

	for (I = Istart; I <= Iend; I++) {
		/* At the bottom boundary: */
		Ari[Jstart-1] = 0;
		Cmri[Jstart-1] = fi[I][Jstart-1];

		for (J = Jstart; J <= Jend; J++) { /* Forward substitution */
			Ari[J]  = aN[I][J]/(aP[I][J] - aS[I][J]*Ari[J-1]); /* eq. 7.6b */
			Cri     = aE[I][J]*fi[I+1][J] + aW[I][J]*fi[I-1][J] + b[I][J];
			Cmri[J] = (aS[I][J]*Cmri[J-1] + Cri)/(aP[I][J] - aS[I][J]*Ari[J-1]); /* eq. 7.6c */
		}

		for (J = Jend; J >= Jstart; J--) /* Back substitution */
			fi[I][J] = Ari[J]*fi[I][J+1] + Cmri[J]; /* eq. 7.6a */
	} /* for I from west */

	/* Solving (n-s) lines from the east */

	for (I = Iend - 1; I >= Istart; I--) {
		/* At the bottom boundary: */
		Ari[Jstart-1] = 0;
		Cmri[Jstart-1] = fi[I][Jstart-1];

		for (J = Jstart; J <= Jend; J++) { /* Forward substitution */
			Ari[J]  = aN[I][J]/(aP[I][J] - aS[I][J]*Ari[J-1]); /* eq. 7.6b */
			Cri     = aE[I][J]*fi[I+1][J] + aW[I][J]*fi[I-1][J] + b[I][J];
			Cmri[J] = (aS[I][J]*Cmri[J-1] + Cri)/(aP[I][J] - aS[I][J]*Ari[J-1]); /* eq. 7.6c */
		}

		for (J = Jend; J >= Jstart; J--) /* Back substitution */
			fi[I][J] = Ari[J]*fi[I][J+1] + Cmri[J]; /* eq. 7.6a */
	} /* for I from east */
	free(Ari);
	free(Cmri);

}

/* ################################################################# */
void solveGS(double **fi, double **b, double **aE, double **aW, double **aN, double **aS, double **aP)
/* ################################################################# */
{
/***** Purpose: To solve the algebraic equation 7.7. using Gauss-Seidel ******/
	int    I, J;

	for (I = Istart; I <= Iend; I++)
		for (J = Jstart; J <= Jend; J++)
			fi[I][J] = ( aE[I][J]*fi[I+1][J  ] + aW[I][J]*fi[I-1][J  ]
			           + aN[I][J]*fi[I  ][J+1] + aS[I][J]*fi[I  ][J-1] + b[I][J])
			           /aP[I][J];
} /* solveGS */

/* ################################################################# */
void solveSOR(double **fi, double **b, double **aE, double **aW, double **aN, double **aS, double **aP)
/* ################################################################# */
{
/***** Purpose: To solve the algebraic equation 7.7. using Successive over-relaxation ******/
	int    I, J;

	for (I = Istart; I <= Iend; I++)
		for (J = Jstart; J <= Jend; J++)
			fi[I][J] = ( aE[I][J]*fi[I+1][J  ] + aW[I][J]*fi[I-1][J  ]
			           + aN[I][J]*fi[I  ][J+1] + aS[I][J]*fi[I  ][J-1] + b[I][J])
			           /aP[I][J] * omega - (omega - 1)*fi[I][J];
} /* solveSOR */

/* ################################################################# */
void conv(void)
/* ################################################################# */
{
/***** Purpose: To calculate the convective mass flux component pr. unit ******/
/*****          area defined in eq. 5.7 ******/
	int    I, J, i, j;

	for (I = 1; I <= NPI + 1; I++) {
		i = I;
		for (J = 1; J <= NPJ + 1; J++) {
			j = J;
			F_u[i][J] = (rho[I-1][J  ]*(x[I] - x_u[i]) + rho[I][J]*(x_u[i] - x[I-1]))*u[i][J]/(x[I] - x[I-1]); /* = F(i, J) */
			F_v[I][j] = (rho[I  ][J-1]*(y[J] - y_v[j]) + rho[I][J]*(y_v[j] - y[J-1]))*v[I][j]/(y[J] - y[J-1]); /* = F(I, j) */
		} /* for J */
	} /* for I */

} /* conv */

/* ################################################################# */
void ucoeff(double **aE, double **aW, double **aN, double **aS, double **aP, double **b)
/* ################################################################# */
{
/***** Purpose: To calculate the coefficients for the u equation. ******/
	int    i, j, I, J;
	double Fw, Fe, Fs, Fn, 
	       Dw, De, Ds, Dn, 
	       AREAw, AREAe, AREAs, AREAn,
	       aPold, mun, mus;

	Istart = 2;
	Iend   = NPI;
	Jstart = 1;
	Jend   = NPJ;

	conv();

	for (I = Istart; I <= Iend; I++) {
		i = I;
		for (J = Jstart; J <= Jend; J++) {
			j = J;
			
			/* Geometrical parameters */
			/* Areas of the cell faces */

			AREAw = y_v[j+1] - y_v[j]; /* See fig. 6.3 */
			AREAe = AREAw;
			AREAs = x[I] - x[I-1];
			AREAn = AREAs;

			/* eq. 6.9a-6.9d - the convective mass flux defined in eq. 5.8a  */
			/* note:  F = rho*u but Fw = (rho*u)w = rho*u*AREAw per definition. */

			Fw = 0.5*(F_u[i  ][J  ] + F_u[i-1][J  ])*AREAw;
			Fe = 0.5*(F_u[i+1][J  ] + F_u[i  ][J  ])*AREAe;
			Fs = 0.5*(F_v[I  ][j  ] + F_v[I-1][j  ])*AREAs;
			Fn = 0.5*(F_v[I  ][j+1] + F_v[I-1][j+1])*AREAn;

			/* eq. 6.9e-6.9h - the transport by diffusion defined in eq. 5.8b  */
			/* note: D = mu/Dx but Dw = (mu/Dx)*AREAw per definition */

			Dw = mueff[I-1][J]/(x_u[i  ] - x_u[i-1])*AREAw;
			De = mueff[I  ][J]/(x_u[i+1] - x_u[i  ])*AREAe;
			Ds = 0.25*(mueff[I-1][J  ] + mueff[I][J  ] + mueff[I-1][J-1] + mueff[I][J-1])/(y[J  ] - y[J-1])*AREAs;
			Dn = 0.25*(mueff[I-1][J+1] + mueff[I][J+1] + mueff[I-1][J  ] + mueff[I][J  ])/(y[J+1] - y[J  ])*AREAn;

			/* The source terms */

			mus = 0.25*(mueff[I][J] + mueff[I-1][J] + mueff[I][J-1] + mueff[I-1][J-1]);
			mun = 0.25*(mueff[I][J] + mueff[I-1][J] + mueff[I][J+1] + mueff[I-1][J+1]);
			
			// Calculate sourceterm in u-direction:
			if (CONS[I][J][FIXED] == true) {
				SP[i][J] = - LARGE;
			}	
			else if (CONS[I][J][YPLUS] == true) {
				if(yplus[I][J] < 11.63)
					SP[i][J]= -mu[I][J]*AREAs/(0.5*AREAw);
				else
					SP[i][J]= -rho[I][J] * pow(Cmu, 0.25) * sqrt(k[I][J]) / uplus[I][J] * AREAs;
			}
			else
				SP[i][J] = 0.;  /* Source term Su, see page 451. */			

			Su[i][J] = (mueff[I][J]*dudx[I][J] - mueff[I-1][J]*dudx[I-1][J]) / (x[I] - x[I-1]) + 
			           (mun        *dvdx[i][j+1] - mus        *dvdx[i][j]) / (y_v[j+1] - y_v[j]) -
                       2./3. * (rho[I][J]*k[I][J] - rho[I-1][J]*k[I-1][J])/(x[I] - x[I-1]);
			Su[I][j] *= AREAw*AREAs;
			
			/* The coefficients (hybrid differencing sheme) */

			aW[i][J] = max3( Fw, Dw + 0.5*Fw, 0.);
			aE[i][J] = max3(-Fe, De - 0.5*Fe, 0.);
			
			/* aS, check current position and for wall to the south (J-1) */
			if (CONS[I][J][YPLUS] == true && CONS[I][J-1][FIXED] == true) aS[i][J]=0.;
			else      aS[i][J] = max3( Fs, Ds + 0.5*Fs, 0.);
            
			/* aN, check current position and for wall to the north (J+1) */
			if (CONS[I][J][YPLUS] == true && CONS[I][J+1][FIXED] == true) aN[i][J] =0.;
			else        aN[i][J] = max3(-Fn, Dn - 0.5*Fn, 0.);
            
			aPold    = 0.5*(rho[I-1][J] + rho[I][J])*AREAe*AREAn/Dt;

			/* eq. 8.31 without time dependent terms (see also eq. 5.14): */

			aP[i][J] = aW[i][J] + aE[i][J] + aS[i][J] + aN[i][J] + Fe - Fw + Fn - Fs - SP[I][J] + aPold;

			/* Calculation of d[i][J] = d_u[i][J] defined in eq. 6.23 for use in the  */
			/* equation for pression correction (eq. 6.32). See subroutine pccoeff. */

			d_u[i][J] = AREAw*relax_u/aP[i][J];

			/* Putting the integrated pressure gradient into the source term b[i][J] */
			/* The reason is to get an equation on the generalised form  */
			/* (eq. 7.7 ) to be solved by the TDMA algorithm.  */
			/* note: In reality b = a0p*fiP + Su = 0.  */

			b[i][J] = (p[I-1][J] - p[I][J])*AREAw + Su[I][J] + aPold*u_old[i][J];

			/* Introducing relaxation by eq. 6.36 . and putting also the last  */
			/* term on the right side into the source term b[i][J] */

			aP[i][J] /= relax_u;
			b [i][J] += (1 - relax_u)*aP[i][J]*u[i][J];

			/* now we have implemented eq. 6.36 in the form of eq. 7.7 */
			/* and the TDMA algorithm can be called to solve it. This is done  */
			/* in the next step of the main program. */

			} /* for j */
		} /* for i */
	
} /* ucoeff */

/* ################################################################# */
void vcoeff(double **aE, double **aW, double **aN, double **aS, double **aP, double **b)
/* ################################################################# */
{
/***** Purpose: To calculate the coefficients for the v equation. ******/
	int    i, j, I, J;
	double Fw, Fe, Fs, Fn, 
	       Dw, De, Ds, Dn, 
	       AREAw, AREAe, AREAs, AREAn,
	       aPold, mue, muw;

	Istart = 1;
	Iend   = NPI;
	Jstart = 2;
	Jend   = NPJ;

	conv();

	for (I = Istart; I <= Iend; I++) {
		i = I;
		for (J = Jstart; J <= Jend; J++) {
			j = J;

			/* Geometrical parameters */
			/* Areas of the cell faces */

			AREAw = y[J] - y[J-1]; /* See fig. 6.4 */
			AREAe = AREAw;
			AREAs = x_u[i+1] - x_u[i];
			AREAn = AREAs;

			/* eq. 6.11a-6.11d - the convective mass flux defined in eq. 5.8a  */
			/* note:  F = rho*u but Fw = (rho*u)w = rho*u*AREAw per definition. */

			Fw = 0.5*(F_u[i  ][J] + F_u[i  ][J-1])*AREAw;
			Fe = 0.5*(F_u[i+1][J] + F_u[i+1][J-1])*AREAe;
			Fs = 0.5*(F_v[I  ][j] + F_v[I  ][j-1])*AREAs;
			Fn = 0.5*(F_v[I  ][j] + F_v[I  ][j+1])*AREAn;

			/* eq. 6.11e-6.11h - the transport by diffusion defined in eq. 5.8b */
			/* note: D = mu/Dx but Dw = (mu/Dx)*AREAw per definition */

			Dw = 0.25*(mueff[I-1][J-1] + mueff[I  ][J-1] + mueff[I-1][J] + mueff[I  ][J])/(x[I  ] - x[I-1])*AREAw;
			De = 0.25*(mueff[I  ][J-1] + mueff[I+1][J-1] + mueff[I  ][J] + mueff[I+1][J])/(x[I+1] - x[I  ])*AREAe;
			Ds = mueff[I][J-1]/(y_v[j  ] - y_v[j-1])*AREAs;
			Dn = mueff[I][J  ]/(y_v[j+1] - y_v[j  ])*AREAn;

			/* The source terms */

			muw = 0.25*(mueff[I][J] + mueff[I-1][J] + mueff[I][J-1] + mueff[I-1][J-1]);
			mue = 0.25*(mueff[I][J] + mueff[I+1][J] + mueff[I][J-1] + mueff[I+1][J-1]);

			// Calculate sourceterm in v-direction:
			if (CONS[I][J][FIXED] == true) {
				SP[I][j] = -LARGE;
			}
			else if (CONS[I][J][XPLUS] == true) {
				if(xplus[I][J] < 11.63)
					SP[i][J]  = -mu[I][J]*AREAw/(0.5*AREAs);
				else
					SP[i][J]  = -rho[I][J] * pow(Cmu, 0.25) * sqrt(k[I][J]) / vplus[I][J] * AREAw;
			}
			// Set source terms to 0 if no vertical boundaries:
			else
				SP[I][j] = 0.; 

			Su[I][j] = (mueff[I][J]*dvdy[I][J] - mueff[I][J-1]*dvdy[I][J-1])/(y[J] - y[J-1]) + 
			           (mue*dudy[i+1][j] - muw*dudy[i][j])/(x_u[i+1] - x_u[i]) - 
                       2./3. * (rho[I][J]*k[I][J] - rho[I][J-1]*k[I][J-1])/(y[J] - y[J-1]); 

			Su[I][j] *= AREAw*AREAs;

			/* The coefficients (hybrid differencing scheme) */
			
			/* aW, check current position and for wall to the west (I-1) */
			if (CONS[I][J][XPLUS] == true && CONS[I-1][J][FIXED] == true) aW[I][j]=0.;
			else      aW[I][j] = max3( Fw, Dw + 0.5*Fw, 0.);
            
			/* aE, check current position and for wall to the east (I+1) */
			if (CONS[I][J][XPLUS] == true && CONS[I+1][J][FIXED] == true) aE[I][j]=0.;
			else      aE[I][j] = max3(-Fe, De - 0.5*Fe, 0.);

			aS[I][j] = max3( Fs, Ds + 0.5*Fs, 0.);
			aN[I][j] = max3(-Fn, Dn - 0.5*Fn, 0.);
			
			aPold    = 0.5*(rho[I][J-1] + rho[I][J])*AREAe*AREAn/Dt;

			/* eq. 8.31 without time dependent terms (see also eq. 5.14): */

			aP[I][j] = aW[I][j] + aE[I][j] + aS[I][j] + aN[I][j] + Fe - Fw + Fn - Fs - SP[I][J] + aPold;

			/* Calculation of d[I][j] = d_v[I][j] defined in eq. 6.23 for use in the */
			/* equation for pression correction (eq. 6.32) (see subroutine pccoeff). */

			d_v[I][j] = AREAs*relax_v/aP[I][j];

			/* Putting the integrated pressure gradient into the source term b[I][j] */
			/* The reason is to get an equation on the generalised form */
			/* (eq. 7.7 ) to be solved by the TDMA algorithm. */
			/* note: In reality b = a0p*fiP + Su = 0. */

			b[I][j] = (p[I][J-1] - p[I][J])*AREAs + Su[I][j] + aPold*v_old[I][j];

			/* Introducing relaxation by eq. 6.37 . and putting also the last */
			/* term on the right side into the source term b[i][J] */

			aP[I][j] /= relax_v;
			b [I][j] += (1 - relax_v)*aP[I][j]*v[I][j];

			/* now we have implemented eq. 6.37 in the form of eq. 7.7 */
			/* and the TDMA algorithm can be called to solve it. This is done */
			/* in the next step of the main program. */

			} /* for j */
		} /* for i */

} /* vcoeff */

/* ################################################################# */
void pccoeff(double **aE, double **aW, double **aN, double **aS, double **aP, double **b)
/* ################################################################# */
{
/***** Purpose: To calculate the coefficients for the pc equation. ******/
	int    i, j, I, J;
	double AREAw, AREAe, AREAs, AREAn;
	double SSUM;

	SMAX = 0.;
	SSUM = 0.;
	SAVG = 0.;

	Istart = 1;
	Iend   = NPI;
	Jstart = 1;
	Jend   = NPJ;

	conv();

	for (I = Istart; I <= Iend; I++) {
		i = I;
		for (J = Jstart; J <= Jend; J++) {
			j = J;

			/* Geometrical parameters */	
			/* Areas of the cell faces */

			AREAw = y_v[j+1] - y_v[j]; /* = A[i][J] See fig. 6.2 or fig. 6.5 */
			AREAe = AREAw;
			AREAs = x_u[i+1] - x_u[i]; /* = A[I][j] */
			AREAn = AREAs;

			/* The constant b� in eq. 6.32 */

			b[I][J] = F_u[i][J]*AREAw - F_u[i+1][J]*AREAe + F_v[I][j]*AREAs - F_v[I][j+1]*AREAn;

			SP[I][J] = 0.;
			Su[I][J] = 0.;
						
			b[I][J] += Su[I][J];

			SMAX     = max2(SMAX,fabs(b[I][J]));
			SSUM    += fabs(b[I][J]);
			
			/* The coefficients */

			aE[I][J] = 0.5*(rho[I  ][J  ] + rho[I+1][J  ])*d_u[i+1][J  ]*AREAe;
			aW[I][J] = 0.5*(rho[I-1][J  ] + rho[I  ][J  ])*d_u[i  ][J  ]*AREAw;
			aN[I][J] = 0.5*(rho[I  ][J  ] + rho[I  ][J+1])*d_v[I  ][j+1]*AREAn;
			aS[I][J] = 0.5*(rho[I  ][J-1] + rho[I  ][J  ])*d_v[I  ][j  ]*AREAs;

			aP[I][J] = aE[I][J] + aW[I][J] + aN[I][J] + aS[I][J] - SP[I][J];

			pc[I][J] = 0.;

			/* note: At the points nearest the boundaries, some coefficients are */
			/* necessarily zero. For instance at I = 1 and J = 1, the coefficients */
			/* aS and aW are zero since they are on the outside of the calculation */
			/* domain. This is automatically satisfied through the initialisation */
			/* where d_u[i][J] and d_v[I][j] are set to zero at these points. */

			} /* for J */
		} /* for I */

		/* Average error in mass balance is summed error devided by */
		/* number of internal grid points */
	      SAVG = SSUM/((Iend - Istart)*(Jend - Jstart));

} /* pccoeff */

/* ################################################################# */
void storeresults(void)
/* ################################################################# */
{
/***** To newly calculated variables are stored in the arrays ******/
/***** for old variables, which can be used in the next timestep ******/
	int I, J, i, j;

	for(i = 2; i <= NPI; i++)
		for(J = 1; J <= NPJ; J++)
			u_old[i][J] = u[i][J];

	for(I = 1; I <= NPI; I++)
		for(j = 2; j <= NPJ; j++)
			v_old[I][j] = v[I][j];

	for(I = 1; I <= NPI; I++)
		for(J = 1; J <= NPJ; J++) {
			pc_old[I][J]  = pc[I][J];
			T_old[I][J]   = T[I][J];
			eps_old[I][J] = eps[I][J];
			k_old[I][J]   = k[I][J];				
      }

} /* storeresults */


/* ################################################################# */
void velcorr(void)
/* ################################################################# */
{
/***** To correct the pressure and the velocities by eq. 6.24, 6.25 ******/
/*****  and a modified version of eq. 6.33. ******/
	int I, J, i, j;

	for (I = 1; I <= NPI; I++) {
		i = I;
		for (J = 1; J <= NPJ; J++) {
			j = J;

			p[I][J] += relax_pc*pc[I][J]; /* equation 6.33 */

			/* Velocity correction */
			/* Note: the relaxation factors for u and v are included  */
			/* in the d_u and d_v terms (see page 146) */

			if (i != 1)
				u[i][J] += d_u[i][J]*(pc[I-1][J  ] - pc[I][J]); /* eq. 6.24 */

			if (j != 1)
				v[I][j] += d_v[I][j]*(pc[I  ][J-1] - pc[I][J]); /* eq. 6.25 */

		} /* for J */
	} /* for I */

} /* velcorr */

/* ################################################################# */
void Tcoeff(double **aE, double **aW, double **aN, double **aS, double **aP, double **b)
/* ################################################################# */
{
/***** Purpose: To calculate the coefficients for the T equation. ******/
	int    i, j, I, J;
	double Fw, Fe, Fs, Fn, 
	       Dw, De, Ds, Dn, 
	       AREAw, AREAe, AREAs, AREAn,
	       aPold;

	Istart = 1;
	Iend   = NPI;
	Jstart = 1;
	Jend   = NPJ;

	conv();
	
	//################BEGIN SELF ADDED CODE################//
    properties(); 
	calc_wall_coeff();
	//#################END SELF ADDED CODE#################//

	for (I = Istart; I <= Iend; I++) {
		i = I;
		for (J = Jstart; J <= Jend; J++) {
			j = J;
			/* Geometrical parameters */
			/* Areas of the cell faces */

			AREAw = y_v[j+1] - y_v[j]; /* = A[i][J] See fig. 6.2 or fig. 6.5 */
			AREAe = AREAw;
			AREAs = x_u[i+1] - x_u[i]; /* = A[I][j] */
			AREAn = AREAs;

			/* The convective mass flux defined in eq. 5.8a */
			/* note:  F = rho*u but Fw = (rho*u)w = rho*u*AREAw per definition. */

			Fw = F_u[i  ][J  ]*AREAw;
			Fe = F_u[i+1][J  ]*AREAe;
			Fs = F_v[I  ][j  ]*AREAs;
			Fn = F_v[I  ][j+1]*AREAn;

			/* The transport by diffusion defined in eq. 5.8b */
			/* note: D = mu/Dx but Dw = (mu/Dx)*AREAw per definition */

			/* The conductivity, Gamma, at the interface is calculated */
			/* with the use of a harmonic mean. */

			Dw = Gamma[I-1][J  ]*Gamma[I  ][J  ]/(Gamma[I-1][J  ]*(x[I  ] - x_u[i  ]) + Gamma[I  ][J  ]*(x_u[i  ] - x[I-1]))*AREAw;
			De = Gamma[I  ][J  ]*Gamma[I+1][J  ]/(Gamma[I  ][J  ]*(x[I+1] - x_u[i+1]) + Gamma[I+1][J  ]*(x_u[i+1] - x[I  ]))*AREAe;
			Ds = Gamma[I  ][J-1]*Gamma[I  ][J  ]/(Gamma[I  ][J-1]*(y[J  ] - y_v[j  ]) + Gamma[I  ][J  ]*(y_v[j  ] - y[J-1]))*AREAs;
			Dn = Gamma[I  ][J  ]*Gamma[I  ][J+1]/(Gamma[I  ][J  ]*(y[J+1] - y_v[j+1]) + Gamma[I  ][J+1]*(y_v[j+1] - y[J  ]))*AREAn;

			/* The coefficients (hybrid differencing scheme) */
			aW[I][J] = max3( Fw, Dw + 0.5*Fw, 0.);
			aE[I][J] = max3(-Fe, De - 0.5*Fe, 0.);
			aS[I][J] = max3( Fs, Ds + 0.5*Fs, 0.);
			aN[I][J] = max3(-Fn, Dn - 0.5*Fn, 0.);

			/* The source terms, page 278 */

			/* On a wall, fix temperature */
			if (CONS[I][J][FIXED]*CONS[I][J][HOT] == true) {
				SP[I][J] = - LARGE;
				Su[I][J] = LARGE*TEMP;
			}
			
			/* Check if one of the source terms is valid */
			else if (CONS[I][J][YPLUS]*CONS[I][J][HOT] == true || CONS[I][J][XPLUS]*CONS[I][J][HOT] == true) {

				/* Calculate sourceterm in u-direction: */
				if (CONS[I][J][YPLUS]*CONS[I][J][HOT] == true) {
					
					if(yplus[I][J] < 11.63) {	/* laminar flow, eq. 9.13 */
						SP_u[I][J] = -mu[I][J]/Prandtl[I][J]*AREAs/(0.5*AREAw);
     		
					}
					else {						/* Turbulent flow, eq. 9.24 */
						SP_u[I][J] = -rho[I][J] * pow(Cmu, 0.25) * sqrt(k[I][J])* AREAs / Tplus_u[I][J];
					
						/* Coefficient aS, check current position and for wall to the south (J-1) */
						if      (CONS[I][J-1][FIXED] == true) aS[I][J] = 0.;
						/* Coefficient aN, check current position and for wall to the north (J+1) */
						else if (CONS[I][J+1][FIXED] == true) aN[I][J] = 0.; 
					}
					/* Source term Su */
					Su_u[I][J] = -SP_u[I][J]*TEMP;

					/* Store source terms */
					SP[I][J] = SP_u[I][J];
					Su[I][J] = Su_u[I][J];

				}
				
				/* Calculate sourceterm in v-direction: */
				if (CONS[I][J][XPLUS]*CONS[I][J][HOT] == true) {

					if(xplus[I][J] < 11.63) { 	/* laminar flow, eq. 9.13 */
						SP_v[I][J] = -mu[I][J]/Prandtl[I][J]*AREAw/(0.5*AREAs);
					}
					else {						/* Turbulent flow, eq. 9.24 */
						SP_v[I][J]  = -rho[I][J] * pow(Cmu, 0.25) * sqrt(k[I][J]) / Tplus_v[I][J] * AREAw;
					
						/* Coefficient aW, check current position and for wall to the west (I-1) */
						if      (CONS[I-1][J][FIXED] == true) aW[I][J] = 0.;
						/* Coefficient aE, check current position and for wall to the east (I+1) */
						else if (CONS[I+1][J][FIXED] == true) aE[I][J] = 0.;
					}

					/* Source term Su */
					Su_v[I][J] = -SP_v[I][J]*TEMP;

					/* Store source terms */
					SP[I][J] = SP_v[I][J];
					Su[I][J] = Su_v[I][J];
				}
				
				/* Calculate sourceterm in u&v-direction: */
				if (CONS[I][J][YPLUS]*CONS[I][J][XPLUS]*CONS[I][J][HOT] == true) {
					
					SP[I][J] = -mag(SP_u[I][J], SP_v[I][J]);	/* SP should always be negative */
					Su[I][J] = mag(Su_u[I][J], Su_v[I][J]);
				}
			}
			else {	/* for adiabatic walls, or places without wall*/
				SP[i][J] = 0.;
				Su[I][J] = 0.;
			}
			
			// Multiply by cell volume
			Su[I][J] *= AREAw*AREAs;
			SP[I][J] *= AREAw*AREAs;

			aPold    = rho[I][J]*AREAe*AREAn/Dt;

			/* eq. 8.31 with time dependent terms (see also eq. 5.14): */

			aP[I][J] = aW[I][J] + aE[I][J] + aS[I][J] + aN[I][J] + Fe - Fw + Fn - Fs - SP[I][J] + aPold;

			/* Setting the source term equal to b */

			b[I][J] = Su[I][J] + aPold*T_old[I][J];

			/* Introducing relaxation by eq. 6.37 . and putting also the last */
			/* term on the right side into the source term b[i][J] */

			aP[I][J] /= relax_T;
			b [I][J] += (1 - relax_T)*aP[I][J]*T[I][J];

			/* now the TDMA algorithm can be called to solve the equation. */
			/* This is done in the next step of the main program. */

			} /* for J */
		} /* for I */

} /* Tcoeff */

/* ################################################################# */
void epscoeff(double **aE, double **aW, double **aN, double **aS, double **aP, double **b)
/* ################################################################# */
{
/***** Purpose: Calculate the epsilon *****/
/*****           *****/
	int    i, j, I, J;
	double Fw, Fe, Fs, Fn, 
	       Dw, De, Ds, Dn, 
	       AREAw, AREAe, AREAs, AREAn,
	       aPold;

	Istart = 1;
	Iend   = NPI;
	Jstart = 1;
	Jend   = NPJ;

	conv();
    properties(); 
    
	for (I = Istart; I <= Iend; I++) {
		i = I;
		for (J = Jstart; J <= Jend; J++) {
			j = J;
			/* Geometrical parameters */
			/* Areas of the cell faces */

			AREAw = y_v[j+1] - y_v[j]; /* = A[i][J] See fig. 6.2 or fig. 6.5 */
			AREAe = AREAw;
			AREAs = x_u[i+1] - x_u[i]; /* = A[I][j] */
			AREAn = AREAs;

			/* The convective mass flux defined in eq. 5.8a */
			/* note:  F = rho*u but Fw = (rho*u)w = rho*u*AREAw per definition. */

			Fw = F_u[i  ][J  ]*AREAw;
			Fe = F_u[i+1][J  ]*AREAe;
			Fs = F_v[I  ][j  ]*AREAs;
			Fn = F_v[I  ][j+1]*AREAn;

			/* The transport by diffusion defined in eq. 5.8b */
			/* note: D = mu/Dx but Dw = (mu/Dx)*AREAw per definition */

			/* The conductivity, Gamma, at the interface is calculated */
			/* with the use of a harmonic mean. */

			Dw = mut[I-1][J  ]*mut[I  ][J  ]/sigmaeps/(mut[I-1][J  ]*(x[I  ] - x_u[i  ]) + mut[I  ][J  ]*(x_u[i  ] - x[I-1]))*AREAw;
			De = mut[I  ][J  ]*mut[I+1][J  ]/sigmaeps/(mut[I  ][J  ]*(x[I+1] - x_u[i+1]) + mut[I+1][J  ]*(x_u[i+1] - x[I  ]))*AREAe;
			Ds = mut[I  ][J-1]*mut[I  ][J  ]/sigmaeps/(mut[I  ][J-1]*(y[J  ] - y_v[j  ]) + mut[I  ][J  ]*(y_v[j  ] - y[J-1]))*AREAs;
			Dn = mut[I  ][J  ]*mut[I  ][J+1]/sigmaeps/(mut[I  ][J  ]*(y[J+1] - y_v[j+1]) + mut[I  ][J+1]*(y_v[j+1] - y[J  ]))*AREAn;

			/* The source terms */
			
			// Set sourceterm to zero on walls or objects:
			if (CONS[I][J][FIXED] == true) {
				Su[I][J]  = 0;
				SP[I][J]  = 0;			
			}
			/* Calculate sourceterm in u&v-direction: */
			else if (CONS[I][J][YPLUS]*CONS[I][J][XPLUS] == true) {
				SP[I][J] = -LARGE;
				Su[I][J] = pow(Cmu,0.75)*pow(k[I][J],1.5)/(kappa*mag(0.5*AREAw,0.5*AREAe))*LARGE;
			}
			/* Calculate sourceterm in u-direction: */
			else if (CONS[I][J][YPLUS] == true) {
				SP[I][J] = -LARGE;
				Su[I][J] = pow(Cmu,0.75)*pow(k[I][J],1.5)/(kappa*0.5*AREAw)*LARGE;
			}
			/* Calculate sourceterm in v-direction: */
			else if (CONS[I][J][XPLUS] == true) {
				SP[I][J] = -LARGE;
				Su[I][J] = pow(Cmu,0.75)*pow(k[I][J],1.5)/(kappa*0.5*AREAe)*LARGE;

			}
			/* Normal grid cell, not on or near object.*/
			else {
				Su[I][J] = C1eps * eps[I][J] / k[I][J] * 2. * mut[I][J] * E2[I][J];
				SP[I][J] = -C2eps * rho[I][J] * eps[I][J] / (k[I][J] + SMALL);
			}			
			                                 
			Su[I][J] *= AREAw*AREAs;
			SP[I][J] *= AREAw*AREAs;

			/* The coefficients (hybrid differencing scheme) */

			aW[I][J] = max3( Fw, Dw + 0.5*Fw, 0.);
			aE[I][J] = max3(-Fe, De - 0.5*Fe, 0.);
			aS[I][J] = max3( Fs, Ds + 0.5*Fs, 0.);
			aN[I][J] = max3(-Fn, Dn - 0.5*Fn, 0.);
			
			aPold    = rho[I][J]*AREAe*AREAn/Dt;

			/* eq. 8.31 with time dependent terms (see also eq. 5.14): */

			aP[I][J] = aW[I][J] + aE[I][J] + aS[I][J] + aN[I][J] + Fe - Fw + Fn - Fs - SP[I][J] + aPold;

			/* Setting the source term equal to b */

			b[I][J] = Su[I][J] + aPold*eps_old[I][J];

			/* now the TDMA algorithm can be called to solve the equation. */
			/* This is done in the next step of the main program. */

			} /* for J */
		} /* for I */
	
} /*epscoef*/

/* ################################################################# */
void kcoeff(double **aE, double **aW, double **aN, double **aS, double **aP, double **b)
/* ################################################################# */
{
/***** Purpose: Calculate the epsilon *****/
/*****           *****/
	int    i, j, I, J;
	double Fw, Fe, Fs, Fn, 
	       Dw, De, Ds, Dn, 
	       AREAw, AREAe, AREAs, AREAn,
	       aPold;

	Istart = 1;
	Iend   = NPI;
	Jstart = 1;
	Jend   = NPJ;

	conv();
   
    properties(); 
    
	calc_wall_coeff(); 
    
	for (I = Istart; I <= Iend; I++) {
		i = I;
		for (J = Jstart; J <= Jend; J++) {
			j = J;
			/* Geometrical parameters */
			/* Areas of the cell faces */

			AREAw = y_v[j+1] - y_v[j]; /* = A[i][J] See fig. 6.2 or fig. 6.5 */
			AREAe = AREAw;
			AREAs = x_u[i+1] - x_u[i]; /* = A[I][j] */
			AREAn = AREAs;

			/* The convective mass flux defined in eq. 5.8a */
			/* note:  F = rho*u but Fw = (rho*u)w = rho*u*AREAw per definition. */

			Fw = F_u[i  ][J  ]*AREAw;
			Fe = F_u[i+1][J  ]*AREAe;
			Fs = F_v[I  ][j  ]*AREAs;
			Fn = F_v[I  ][j+1]*AREAn;

			/* The transport by diffusion defined in eq. 5.8b */
			/* note: D = mu/Dx but Dw = (muw/Dx)*AREAw per definition */

			/* The conductivity, Gamma, at the interface is calculated */
			/* with the use of a harmonic mean. */

			Dw = mut[I-1][J  ]*mut[I  ][J  ]/sigmak/(mut[I-1][J  ]*(x[I  ] - x_u[i  ]) + mut[I  ][J  ]*(x_u[i  ] - x[I-1]))*AREAw;
			De = mut[I  ][J  ]*mut[I+1][J  ]/sigmak/(mut[I  ][J  ]*(x[I+1] - x_u[i+1]) + mut[I+1][J  ]*(x_u[i+1] - x[I  ]))*AREAe;
			Ds = mut[I  ][J-1]*mut[I  ][J  ]/sigmak/(mut[I  ][J-1]*(y[J  ] - y_v[j  ]) + mut[I  ][J  ]*(y_v[j  ] - y[J-1]))*AREAs;
			Dn = mut[I  ][J  ]*mut[I  ][J+1]/sigmak/(mut[I  ][J  ]*(y[J+1] - y_v[j+1]) + mut[I  ][J+1]*(y_v[j+1] - y[J  ]))*AREAn;
			
			/* The source terms */
			
			// Set sourceterm to zero on walls or objects:
			if (CONS[I][J][FIXED] == true) {
				Su[I][J]  = 0;
				SP[I][J]  = 0;			
			}
			/* Check if one of the source terms is valid */
			else if (CONS[I][J][YPLUS] == true || CONS[I][J][XPLUS] == true) {

				// Calculate sourceterm in u-direction:
				if (CONS[I][J][YPLUS] == true) {
					/* store in SP_u and SP, to calc magnitude when both xplus and yplus are active */
					SP_u[I][J] = -rho[I][J] * pow(Cmu,0.75) * sqrt(k[I][J]) * uplus[I][J]/(0.5*AREAw) * AREAs * AREAw;
					SP[I][J] = SP_u[I][J];
					Su_u[I][J] = tw[I][J] * 0.5 * (u[i][J] + u[i+1][J])/(0.5*AREAw) * AREAs * AREAw;
					Su[I][J] = Su_u[I][J];
				}
				// Calculate sourceterm in v-direction:
				if (CONS[I][J][XPLUS] == true) {
					SP_v[I][J] = -rho[I][J] * pow(Cmu,0.75) * sqrt(k[I][J]) * vplus[I][J]/(0.5*AREAs) * AREAs * AREAw;
					SP[I][J] = SP_v[I][J];
					Su_v[I][J] = tw[I][J] * 0.5 * (v[I][j] + v[I][j+1])/(0.5*AREAs) * AREAs * AREAw;
					Su[I][J] = Su_v[I][J];
				}
				
				// Calculate magnitude of source terms if both valid: k gaat uit van de wall shear stress. Je zou deze kunnen schrijven als de lengte van de vector (du/dy, dv/dx).
				if (CONS[I][J][YPLUS]*CONS[I][J][XPLUS] == true) {
					SP[I][J] = -rho[I][J] * pow(Cmu,0.75) * sqrt(k[I][J]) * mag(vplus[I][J], uplus[I][J])/(mag(0.5*AREAs,0.5*AREAw)) * AREAs * AREAw;
					Su[I][J] = tw[I][J] * 0.5 * mag((u[i][J] + u[i+1][J]),(v[I][j] + v[I][j+1]))/mag(0.5*AREAs,0.5*AREAw) * AREAs * AREAw;
//					SP[I][J] = mag(SP_u[I][J], SP_v[I][J]);
//					Su[I][J] = mag(Su_u[I][J], Su_v[I][J]);
				}
			}
			else {
				Su[I][J]  = 2. * mut[I][J] * E2[I][J];
				SP[I][J]  = -rho[I][J] * eps[I][J] / k[I][J];
			}                
			
 			Su[I][j] *= AREAw*AREAs;
			SP[I][j] *= AREAw*AREAs;

			/* The coefficients (hybrid differencing sheme) */

			/* aW, check current position and for wall to the west (I-1) */
			if (CONS[I][J][XPLUS] == true && CONS[I-1][J][FIXED] == true) aW[I][j]=0.;
			else      aW[I][j] = max3( Fw, Dw + 0.5*Fw, 0.);
            
			/* aE, check current position and for wall to the east (I+1) */
			if (CONS[I][J][XPLUS] == true && CONS[I+1][J][FIXED] == true) aE[I][j]=0.;
			else      aE[I][j] = max3(-Fe, De - 0.5*Fe, 0.);

			/* aS, check current position and for wall to the south (J-1) */
			if (CONS[I][J][YPLUS] == true && CONS[I][J-1][FIXED] == true) aS[i][J]=0.;
			else      aS[i][J] = max3( Fs, Ds + 0.5*Fs, 0.);
            
			/* aN, check current position and for wall to the north (J+1) */
			if (CONS[I][J][YPLUS] == true && CONS[I][J+1][FIXED] == true) aN[i][J] =0.;
			else      aN[i][J] = max3(-Fn, Dn - 0.5*Fn, 0.);
            
            aPold    = rho[I][J]*AREAe*AREAn/Dt;

			/* eq. 8.31 with time dependent terms (see also eq. 5.14): */

			aP[I][J] = aW[I][J] + aE[I][J] + aS[I][J] + aN[I][J] + Fe - Fw + Fn - Fs - SP[I][J] + aPold;

			/* Setting the source term equal to b */

			b[I][J] = Su[I][J] + aPold*k_old[I][J];

			/* now the TDMA algorithm can be called to solve the equation. */
			/* This is done in the next step of the main program. */

			} /* for J */
		} /* for I */
	
} /*kcoeff*/

/* ################################################################# */
void calc_wall_coeff(void)
/* ################################################################# */
{
/***** Purpose: Calculate uplus, yplus and tw  ******/
// Pim: Calculation of uplus is needed for the source terms Su and Sp 

	int    i, j, I, J;
	double Dx, Dy;
	
	properties(); 
	
	Dx = XMAX/NPI;
	Dy = YMAX/NPJ;
	
	for (I = Istart; I <= Iend; I++) {
		i = I;
		for (J = Jstart; J <= Jend; J++) {
			j = J;
	    	
	    	// Calculate: yplus
	        if (CONS[I][J][YPLUS] == true) {
				        
	        	if (yplus[I][J] < 11.63) { // PIM: Initialise yplus instead of yplus1!
                	tw[I][J]       = mu[I][J] * 0.5 * (u[i][J]+u[i+1][J]) / (0.5*Dy); // PIM: in general holds: 0.5*Dy = y[1] -y[0]
                  	yplus[I][J]  = sqrt(rho[I][J] * fabs(tw[I][J])) * (0.5*Dy) / mu[I][J];
                  	uplus[I][J]  = yplus[I][J];
                  	
            	}/* if */
            	else {
                  	tw[I][J]       = rho[I][J] * pow(Cmu,0.25) * sqrt(k[I][J]) * 0.5 * (u[i][J]+u[i+1][J]) / uplus[I][J];
                  	yplus [I][J] = sqrt(rho[I][J] * fabs(tw[I][J])) * (0.5*Dy) / mu[I][J];
                  	uplus [I][J] = log(ERough*yplus[I][J])/kappa;
					Tplus_u [I][J] = turb_Prandtl*(uplus[I][J] + Pee(turb_Prandtl, Prandtl[I][J])); /* Turbulet T+, eq. 3.50 */
					
					// Test Pee_u (never used)
                  	Pee_u   [I][J] = Pee(turb_Prandtl, Prandtl[I][J]);
            	}/* else */
	        } /* if */
	        
	        // Calculate: xplus
	        if (CONS[I][J][XPLUS] == true) {
				        
	        	if (xplus[I][J] < 11.63) { 
                	tw[I][J]       = mu[I][J] * 0.5 * (v[I][j]+v[I][j+1]) / (0.5*Dx); 
                  	xplus[I][J]  = sqrt(rho[I][J] * fabs(tw[I][J])) * (0.5*Dx) / mu[I][J];
                  	vplus[I][J]  = xplus[I][J];
                  	
            	}/* if */
            	else {
                  	tw[I][J]       = rho[I][J] * pow(Cmu,0.25) * sqrt(k[I][J]) * 0.5 * (v[I][j]+v[I][j+1]) / vplus[I][J];
                  	xplus [I][J] = sqrt(rho[I][J] * fabs(tw[I][J])) * (0.5*Dx) / mu[I][J];
                  	vplus [I][J] = log(ERough*xplus[I][J])/kappa;
					Tplus_v [I][J] = turb_Prandtl*(vplus[I][J] + Pee(turb_Prandtl, Prandtl[I][J])); /* Turbulet T+, eq. 3.50 */
                  	
                  	// Test Pee_v (never used)
                  	Pee_v   [I][J] = Pee(turb_Prandtl, Prandtl[I][J]);
            	}/* else */
	        } /* if */

        } /* for */
	} /* for */
} /* cuplus */

/* ################################################################# */
void properties(void)
/* ################################################################# */
{
/***** Purpose: Calculate the properties in the fluid as a function of temperature *****/
	int   I, J;

	for (I = 0; I <= NPI; I++) { // PIM: why not NPI+1?
		for (J = 1; J <= NPJ + 1; J++) { // PIM: why not J=0?
			// Calculate the density of air as function of pressure and temperature (does not converge)		
//			p_abs  [I][J] = P_ATM + p[I][J]; /* Absolute pressure [Pa]*/
//			rho    [I][J] = 28.964*p_abs[I][J]/T[I][J]/GAS_CONS;      /* Density, temperature and pressure dependent */
//			rho    [I][J] = 28.964*P_ATM/T[I][J]/GAS_CONS;      /* Density, temperature and constant pressure */
			
			// Calculate properties of air as function of temperature
//			mu     [I][J] = 352.975/T[I][J]*0.000001*(0.0264*pow(((T[I][J]-273.15)+50),1.24)+10);    /* Viscosity */
//			Cp     [I][J] = 1031.311-0.2028999*T[I][J]+0.0004005271*T[I][J]*T[I][J];     /* J/(K*kg) Heat capacity - assumed constant for this problem */
//			lambda [I][J] = 0.005+T[I][J]/13944;     /* W/(K*m) Thermal conductivity */

			// Calculate the viscosity in the fluid as a function of temperature
            mut[I][J] = rho[I][J]*Cmu*sqr(k[I][J])/(eps[I][J]+SMALL);
			mueff[I][J] = mu[I][J] + mut[I][J];
			
			// Calculate the laminar, turbulent and effective Gamma
			Gamma_L[I][J] = lambda[I][J]/Cp[I][J]; 			/* Thermal conductivity divided by heat capacity for laminar case*/
			Gamma_T[I][J] = mut[I][J]/turb_Prandtl;			/* Thermal conductivity divided by heat capacity for turbulent case*/
			Gamma  [I][J] = Gamma_L[I][J]+Gamma_T[I][J];	/* Thermal conductivity divided by heat capacity (effective)*/
			
			// Calculate laminair (or moleculair) Prandtl Number
			Prandtl[I][J] = mu[I][J]/Gamma_L[I][J];
        } /* for */
    }
} /* properties */

/* ################################################################# */
void printConv(double time, int iter)
/* ################################################################# */
{
/***** Purpose: Creating result table ******/
    if (time == Dt)
		printf ("ITER  Time\t   u\t       v\t   T\t       SMAX\t   SAVG\n");

	printf ("%4d %11.3e %11.2e %11.2e %11.2e %11.2e %11.2e\n", 
             iter, time, u[NPI][NPJ/2], v[NPI][NPJ/2], T[NPI][NPJ/2], SMAX, SAVG);

} /* printConv */

/* ################################################################# */
void output(void)
/* ################################################################# */
{
/***** Purpose: Creating result table ******/
	int    I, J, i, j;
	double ugrid, vgrid,stream,vorticity;
	FILE   *fp, *str, *velu, *velv, *vort;

/* Plot all results in output.dat */

	fp = fopen("output.dat", "w");

	for (I = 0; I <= NPI; I++) {
		i = I;
		for (J = 1; J <= NPJ; J++) {
			j = J;
			ugrid = 0.5*(u[i][J]+u[i+1][J  ]);
			vgrid = 0.5*(v[I][j]+v[I  ][j+1]);

			/* On an object set pressure to zero */
			if (CONS[I][J][FIXED] == true) {	
				p[I][J] = 0;
			}
			fprintf(fp, "%11.5e\t%11.5e\t%11.5e\t%11.5e\t%11.5e\t%11.5e\t%11.5e\t%11.5e\t%11.5e\t%11.5e\t%11.5e\t%11.5e\t%11.5e\t%11.5e\t%11.5e\t%11.5e\t%11.5e\t%11.5e\t%11.5e\n",
			             x[I], y[J], ugrid, vgrid, p[I][J], T[I][J], rho[I][J], mu[I][J], Gamma[I][J], k[I][J], eps[I][J], Tplus_u[I][J], Tplus_v[I][J], yplus[I][J], xplus[I][J], uplus[I][J], vplus[I][J], Pee_u[I][J], Pee_v[I][J]);
//			             1     2     3      4      5        6        7          8         9            10       11         12             13             14             15             16             17              18           19 
		} /* for J */
		fprintf(fp, "\n");
	} /* for I */

	fclose(fp);

/* Plot vorticity in vort.dat */

	vort = fopen("vort.dat", "w");

	for (I = 1; I <= NPI; I++) {
		i = I;
		for (J = 1; J <= NPJ; J++) {
			j = J;
			vorticity = (u[i][J] - u[i][J-1]) / (y[J] - y[J-1]) - (v[I][j] - v[I-1][j]) / (x[I] - x[I-1]);
			fprintf(vort, "%11.5e\t%11.5e\t%11.5e\n",
			             x[I], y[J], vorticity);
		} /* for J */
		fprintf(vort, "\n");
	} /* for I */

	fclose(vort);

/* Plot streamlines in str.dat */

	str = fopen("str.dat", "w");

	for (I = 0; I <= NPI; I++) {
		i = I;
		for (J = 0; J <= NPJ; J++) {
			j = J;
			stream = -0.5*(v[I+1][j]+v[I][j])*(x[I+1]-x[I])+0.5*(u[i][J+1]+u[i][J])*(y[J+1]-y[J]);
			fprintf(str, "%11.5e\t%11.5e\t%11.5e\n",
			             x[I], y[J], stream);
		} /* for J */
		fprintf(str,  "\n");
	} /* for I */

	fclose(str);

/* Plot horizontal velocity components in velu.dat */

	velu = fopen("velu.dat", "w");

	for (I = 1; I <= NPI+1; I++) {
		i = I;
		for (J = 0; J <= NPJ+1; J++) {
			fprintf(velu, "%11.5e\t%11.5e\t%11.5e\n",
			             x_u[i], y[J], u[i][J]);
		} /* for J */
		fprintf(velu, "\n");
	} /* for I */

	fclose(velu);

/* Plot vertical velocity components in velv.dat */

	velv = fopen("velv.dat", "w");

	for (I = 0; I <= NPI+1; I++) {
		for (J = 1; J <= NPJ+1; J++) {
			j = J;
			fprintf(velv, "%11.5e\t%11.5e\t%11.5e\n",
			             x[I], y_v[j], v[I][j]);
		} /* for J */
		fprintf(velv, "\n");
	} /* for I */

	fclose(velv);

} /* output */

/* ################################################################# */
int *int_1D_array(int np)
/* ################################################################# */
{
/* create an 1D array with size [np] of type int */
	int *a;

	a = (int *) calloc(np, sizeof(int));

	return a;

} /* int_1D_array */

/* ################################################################# */
double *double_1D_array(int np)
/* ################################################################# */
{
/* create an 1D array with size [np] of type double */
	double *a;

	a = (double *) calloc(np, sizeof(double));

	return a;

} /* double_1D_array */

/* ################################################################# */
double **double_2D_matrix (int nm, int np)
/* ################################################################# */
{
/* create an 2D matrix with size [nm, np] of type double */
	int i;
	double **m;

	m = (double **) calloc(nm, sizeof(double *));
	for ( i = 0; i < nm; i++)
		m[i] = (double *) calloc(np, sizeof(double));

	return m;

} /* double_2D_matrix */

/* ################################################################# */
bool ***bool_3D_matrix (int nm, int np, int dpt)
/* ################################################################# */
{
/* create an 3D matrix with size [nm, np, dpt] of type bool */
 	int i,j;
 	bool ***m;

	m = (bool ***) calloc(nm, sizeof(bool**));
    for (i = 0; i< nm; i++) {
    	m[i] = (bool **) calloc(np, sizeof(bool *));
        for (j = 0; j < np; j++)
        	m[i][j] = (bool *)calloc(dpt, sizeof(bool));
    }
	
	return m;

} /* bool_3D_matrix */

/* ################################################################# */
void memalloc(void)
/* ################################################################# */
{
	x    = double_1D_array(NPI + 2);
	x_u  = double_1D_array(NPI + 2);
	y    = double_1D_array(NPJ + 2);
	y_v  = double_1D_array(NPJ + 2);

	u      = double_2D_matrix(NPI + 2, NPJ + 2);
	v      = double_2D_matrix(NPI + 2, NPJ + 2);
	pc     = double_2D_matrix(NPI + 2, NPJ + 2);
	p      = double_2D_matrix(NPI + 2, NPJ + 2);
	T      = double_2D_matrix(NPI + 2, NPJ + 2);
	rho    = double_2D_matrix(NPI + 2, NPJ + 2);
	mu     = double_2D_matrix(NPI + 2, NPJ + 2);
	mut    = double_2D_matrix(NPI + 2, NPJ + 2);
	mueff  = double_2D_matrix(NPI + 2, NPJ + 2);
	Gamma  = double_2D_matrix(NPI + 2, NPJ + 2);
	Cp     = double_2D_matrix(NPI + 2, NPJ + 2);
	k      = double_2D_matrix(NPI + 2, NPJ + 2);
	eps    = double_2D_matrix(NPI + 2, NPJ + 2);
    delta  = double_2D_matrix(NPI + 2, NPJ + 2);
	E      = double_2D_matrix(NPI + 2, NPJ + 2);
	E2     = double_2D_matrix(NPI + 2, NPJ + 2);
	yplus  = double_2D_matrix(NPI + 2, NPJ + 2);
	uplus  = double_2D_matrix(NPI + 2, NPJ + 2);
	tw     = double_2D_matrix(NPI + 2, NPJ + 2);

	u_old  = double_2D_matrix(NPI + 2, NPJ + 2);
	v_old  = double_2D_matrix(NPI + 2, NPJ + 2);
	pc_old = double_2D_matrix(NPI + 2, NPJ + 2);
	T_old  = double_2D_matrix(NPI + 2, NPJ + 2);
	k_old  = double_2D_matrix(NPI + 2, NPJ + 2);
	eps_old= double_2D_matrix(NPI + 2, NPJ + 2);

	dudx   = double_2D_matrix(NPI + 2, NPJ + 2);
	dudy   = double_2D_matrix(NPI + 2, NPJ + 2);
	dvdx   = double_2D_matrix(NPI + 2, NPJ + 2);
	dvdy   = double_2D_matrix(NPI + 2, NPJ + 2);

	aP     = double_2D_matrix(NPI + 2, NPJ + 2);
	aE     = double_2D_matrix(NPI + 2, NPJ + 2);
	aW     = double_2D_matrix(NPI + 2, NPJ + 2);
	aN     = double_2D_matrix(NPI + 2, NPJ + 2);
	aS     = double_2D_matrix(NPI + 2, NPJ + 2);
	b      = double_2D_matrix(NPI + 2, NPJ + 2);

	SP     = double_2D_matrix(NPI + 2, NPJ + 2);
	SP_u   = double_2D_matrix(NPI + 2, NPJ + 2);
	SP_v   = double_2D_matrix(NPI + 2, NPJ + 2);
	Su     = double_2D_matrix(NPI + 2, NPJ + 2);
	Su_u   = double_2D_matrix(NPI + 2, NPJ + 2);
	Su_v   = double_2D_matrix(NPI + 2, NPJ + 2);

	F_u    = double_2D_matrix(NPI + 2, NPJ + 2);
	F_v    = double_2D_matrix(NPI + 2, NPJ + 2);

	d_u    = double_2D_matrix(NPI + 2, NPJ + 2);
	d_v    = double_2D_matrix(NPI + 2, NPJ + 2);
	
	yplus = double_2D_matrix(NPI + 2, NPJ + 2);
	xplus = double_2D_matrix(NPI + 2, NPJ + 2);
	uplus = double_2D_matrix(NPI + 2, NPJ + 2);
	vplus = double_2D_matrix(NPI + 2, NPJ + 2);
	Tplus_u = double_2D_matrix(NPI + 2, NPJ + 2); 
	Tplus_v = double_2D_matrix(NPI + 2, NPJ + 2);
	Prandtl = double_2D_matrix(NPI + 2, NPJ + 2);
	Gamma_L = double_2D_matrix(NPI + 2, NPJ + 2);
	Gamma_T = double_2D_matrix(NPI + 2, NPJ + 2);
	
	// TEST	
	Pee_u   = double_2D_matrix(NPI + 2, NPJ + 2);
	Pee_v   = double_2D_matrix(NPI + 2, NPJ + 2);
	   
	/* for Properties: */
	lambda  = double_2D_matrix(NPI + 2, NPJ + 2);
	p_abs   = double_2D_matrix(NPI + 2, NPJ + 2);

} /* memalloc */

/* ################################################################# */
void readInput (char *name)
/* ################################################################# */
{

	FILE        *fp;
	int			ncons, icons;
	int			nyplus, iyplus;
	int			nxplus, ixplus;
	int 		I, J, ad;
	
	// open file to read from
	fp=fopen(name,"r");
	
	fscanf( fp, "%*s %lf", &XMAX );
	fscanf( fp, "%*s %lf", &YMAX );
	fscanf( fp, "%*s %d", &NPI );
	fscanf( fp, "%*s %d", &NPJ );
	fscanf( fp, "%*s %lf", &relax_u );
	fscanf( fp, "%*s %lf", &relax_T );
	fscanf( fp, "%*s %lf", &Dt );
	fscanf( fp, "%*s %lf", &TOTAL_TIME );
	fscanf( fp, "%*s %d", &MAX_ITER );
	fscanf( fp, "%*s %d", &U_ITER );
	fscanf( fp, "%*s %d", &V_ITER );
	fscanf( fp, "%*s %d", &PC_ITER );
	fscanf( fp, "%*s %d", &T_ITER );
	fscanf( fp, "%*s %d", &EPS_ITER );
	fscanf( fp, "%*s %d", &K_ITER );
	fscanf( fp, "%*s %lf", &TZERO );
	fscanf( fp, "%*s %lf", &TEMP );
	fscanf( fp, "%*s %lf", &U_IN );
	fscanf( fp, "%*s %lf", &rho_init );
	fscanf( fp, "%*s %lf", &mu_init );
	fscanf( fp, "%*s %lf", &Cp_init );
	fscanf( fp, "%*s %lf", &lambda_init );

	// print grid parameters
	printf("From text file:\nGrid:           XMAX = %5.2f [m]         YMAX =  %4.2f [m]\n",XMAX,YMAX);
	printf("                 NPI =   %3d              NPJ =   %3d\n",NPI,NPJ);
	printf("Solver:      relax_u =  %4.2f          relax_T =  %4.2f\n",relax_u,relax_T);
	printf("Time:             Dt =%6.4f [s]       t_max =  %5.2f [s]\n",Dt, TOTAL_TIME);
	printf("Iterations: MAX_ITER =  %4d           U_ITER =  %4d           V_ITER =  %4d\n",MAX_ITER,   U_ITER, V_ITER);
	printf("             PC_ITER =  %4d         EPS_ITER =  %4d           K_ITER =  %4d\n", PC_ITER, EPS_ITER, K_ITER);
	printf("Physics:       Temp. =   %3.0f [K]         U_IN = %5.1f [m/s]\n",TEMP,U_IN);
	printf("                 rho =   %3.0f [kg/m^3]      mu = %5.0e [Pa*s]\n",rho_init, mu_init);
	printf("                  Cp =  %4.0f [J/kg/K]       K = %5.3f [W/m/K]\n",Cp_init, lambda_init);  

	// Allocate memory to save constraints
	CONS  = bool_3D_matrix(NPI + 2, NPJ + 2, 4);
	
	// ###########################################
	// Get fully constrained points from text file
	fscanf( fp, "%*s %d", &ncons );

	// loop through items in file	
	for (icons = 0; icons < ncons; icons++) {
		fscanf( fp, " %d  %d  %d", &I, &J, &ad);
		/* position constrain in 0 */
		CONS[I][J][FIXED] = true;
		/* temperature constrain in 3 */
		CONS[I][J][HOT] = ad;
	}

	// Print results
	printf("Constrains:    Fixed =  %4d ",ncons);

	// ###########################################
	// Get yplus points from text file
	fscanf( fp, "%*s %d", &nyplus );

	// loop through items in file	
	for (iyplus = 0; iyplus < nyplus; iyplus++) {
		fscanf( fp, " %d  %d  %d", &I, &J, &ad);
		/* yplus constrain in 1 */
		CONS[I][J][YPLUS] = true;
		/* temperature constrain in 3 */
		CONS[I][J][HOT] = ad;
	}

	// Print results
	printf("           yplus =  %4d ",nyplus);

	// ###########################################
	// Get xplus points from text file
	fscanf( fp, "%*s %d", &nxplus );

	// loop through items in file	
	for (ixplus = 0; ixplus < nxplus; ixplus++) {
		fscanf( fp, " %d  %d  %d", &I, &J, &ad);
		/* yplus constrain in 2 */
		CONS[I][J][XPLUS] = true;
		/* temperature constrain in 3 */
		CONS[I][J][HOT] = ad;
	}

	// Print results
	printf("           xplus =  %4d \n",nxplus);
	printf("\n");
	
} /* readinput */

void animation(int time)
{
/***** Purpose: Creating result table ******/
	int    I, J, i, j;
	double count, ugrid, vgrid,stream,vorticity, Tavg;
	
	FILE   *f_p, *str, *velu, *velv, *vort;
/* Plot all results in output.dat */
	
	if(timestore==1){
	
	char index[10];
//    sprintf(index, "T%d_%d.dat", nloop, time);
    sprintf(index, "Animation/T%d.dat", nloop);
    
	f_p = fopen(index, "w");
	
	count = 1.;
	
	for (I = 0; I <= NPI; I++) {
		i = I;
		for (J = 1; J <= NPJ; J++) {
			j = J;
			ugrid = 0.5*(u[i][J]+u[i+1][J  ]); /* interpolated horizontal velocity */
			vgrid = 0.5*(v[I][j]+v[I  ][j+1]); /* interpolated vertical velocity */

			/* On an object set pressure to zero */
			if (CONS[I][J][FIXED] == true) {	
				p[I][J] = 0;
			}
			//################BEGIN SELF ADDED CODE################//
			fprintf(f_p, "%11.5e\t%11.5e\t%11.5e\t%11.5e\t%11.5e\t%11.5e\t%11.5e\t%11.5e\t%11.5e\t%11.5e\t%11.5e\t%11.5e\t%11.5e\t%11.5e\t%11.5e\t%11.5e\t%11.5e\t%11.5e\t%11.5e\n",
			             x[I], y[J], ugrid, vgrid, p[I][J], T[I][J], rho[I][J], mu[I][J], Gamma[I][J], k[I][J], eps[I][J], Tplus_u[I][J], Tplus_v[I][J], yplus[I][J], xplus[I][J], uplus[I][J], vplus[I][J], Pee_u[I][J], Pee_v[I][J]);
//			             1     2     3      4      5        6        7          8         9            10       11         12             13             14             15             16             17              18           19
			//#################END SELF ADDED CODE#################//  
		} /* for J */
		fprintf(f_p, "\n");
	} /* for I */

	fclose(f_p);
	timestore=0;
	nloop+=1;
	}
timestore+=1;
}

