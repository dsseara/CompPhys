#include <stdio.h>
#include <math.h>

double G = 1.4;  //Adiabatic heat index

double L = 1.0;  //Length of box
double H = 1.0;  //Height of box
int COLS = 200;  //Columns
int ROWS = 200;  //Rows
int   Nq = 4;   //Number of conserved quantities
int   Ng = 2;   //Number of ghost zones
double T = 0.2; //Total time integrated over


enum{PPP, RHO, VXX, VYY}; //pressure, mass density, velocity in x, velocity in y
enum{DEN, MOX, MOY, EEE}; //mass density, momentum in x, momentum in y, energy
enum{XXX, YYY};
enum{MM_FLUX, PX_FLUX, PY_FLUX, EE_FLUX};

double getx (int j){
	return( L*( ( (double)(j-Ng) + 0.5) / (double)COLS ));
}

double gety (int i){
	return( H*( ( (double)(i-Ng) + 0.5) / (double)ROWS) );
}

//given a U vector in zone i, get the corresponding
//primitive variables. prim = (p, rho, vx, vy)
consToPrim(double * u, double * prim){

	double gamma = G;

	//Conserved densities
	double mass    = u[DEN];
	double momentX = u[MOX];  //Momentum in x direction
	double momentY = u[MOY];  //Momentum in y direction
	double E       = u[EEE];
	//Primitive variables to find
	double p, rho, vX, vY;

	rho = mass;
	vX  = momentX / mass;
	vY  = momentY / mass;
	p   = (gamma - 1.0)*(E - 0.5*mass*(vX*vX + vY*vY));

	prim[PPP] = p;
	prim[RHO] = rho;
	prim[VXX] = vX;
	prim[VYY] = vY;
}

//Convert an array of primitive variables (pressure, rho, vx, vy)
//to the appropriate conserved variables (rho, momentumX, momentumY, energy)
primToCons(double * prim, double * u){
	
	double gamma = G;

	double P   = prim[PPP];
	double rho = prim[RHO];
	double vx  = prim[VXX];
	double vy  = prim[VYY];
	
	//Conserved quantities to find
	double density, momentX, momentY, Energy;

	density = rho;
	momentX = rho*vx;
	momentY = rho*vy;
	Energy  = P/(gamma - 1) + 0.5 * rho * (vx*vx + vy*vy);

	u[DEN] = density;
	u[MOX] = momentX;
	u[MOY] = momentY;
	u[EEE] = Energy;
}

void alpha(double * alphX, double * alphY,  double * prim){
	
	double gamma = 1.4;
	int i, j, q;
	double maxP = -1e4;
	double maxM = -1e4;
	double alphPset[3], alphMset[3];
	double pL[4], pR[4];
	double vL, vR, vB, vA; //velocities left, right, below, above
	
	//X Direction
	for(i=0; i<(ROWS+2*Ng); ++i){
		for(j=0; j<(COLS+2*Ng-1); ++j){
			//Get left and right cells at boundays
			for(q=0; q<Nq; ++q){
				int jiq = Nq*(COLS+2*Ng)*i+ Nq*j +q; //Gets to cell (j,i), and 
				pL[q] = prim[jiq];	      //quantity q in that cell

				pR[q] = prim[jiq + Nq]; //need Nq(COLS*i + j+1) + q
			}

			//Speeds of sound on left and right
			double csL = sqrt(gamma * pL[PPP]/pL[RHO]);
			double csR = sqrt(gamma * pR[PPP]/pR[RHO]);
			
			vL = pL[VXX];
			vR = pR[VXX];

			double lambPL = vL + csL;
			double lambML = vL - csL;
			double lambPR = vR + csR;
			double lambMR = vR - csR;
			
			alphPset[0] = 0.0;
			alphPset[1] = lambPL;
			alphPset[2] = lambPR;
			alphMset[0] = 0.0;
			alphMset[1] = -lambML;
			alphMset[2] = -lambMR;

			for(q = 0; q<3; ++q){
				if (alphPset[q] > maxP) maxP = alphPset[q];
				if (alphMset[q] > maxM) maxM = alphMset[q];
			}
			
			alphX[2*i*(COLS + 2*Ng - 1) + 2*j  ] = maxM;
			alphX[2*i*(COLS + 2*Ng - 1) + 2*j+1] = maxP;
		}
	}

	//Y Direction
	for(j=0; j<(COLS+2*Ng); ++j){
		for(i=0; i<(ROWS+2*Ng-1); ++i){
			//Get left and right cells at boundays
			for(q=0; q<Nq; ++q){
				int jiq = Nq*(COLS+2*Ng)*i + Nq*j + q; //Gets to cell (j,i), and 
				pL[q] = prim[jiq];	      //quantity q in that cell

				pR[q] = prim[jiq + Nq*(COLS+2*Ng)]; //need Nq(COLS*(i+1) + j) + q
			}
			//Speeds of sound on left and right
			double csL = sqrt(gamma * pL[PPP]/pL[RHO]);
			double csR = sqrt(gamma * pR[PPP]/pR[RHO]);
			
			vL = pL[VYY];
			vR = pR[VYY];

			double lambPL = vL + csL;
			double lambML = vL - csL;
			double lambPR = vR + csR;
			double lambMR = vR - csR;
			
			alphPset[0] = 0.0;
			alphPset[1] =  lambPL;
			alphPset[2] =  lambPR;
			alphMset[0] = 0.0;
			alphMset[1] = -lambML;
			alphMset[2] = -lambMR;

			for(q = 0; q<3; ++q){
				if (alphPset[q] > maxP) maxP = alphPset[q];
				if (alphMset[q] > maxM) maxM = alphMset[q];
			}
			
			alphY[2*j*(ROWS + 2*Ng -1) + 2*i  ] = maxM;
			alphY[2*j*(ROWS + 2*Ng -1) + 2*i+1] = maxP;
		}
	}
}


void set_initial(double * prim){
	int i, j, q;
	double gamma = G;

	/*Initial left half of shock tube****/
	double left[4] = {1.0, 1.0, 0.0, 0.0};
	/*Initial right half of shock tube***/
	double right[4] = {0.1, 0.125, 0.0, 0.0};

	
	for(i =0 ; i<(ROWS+2*Ng) ; i++){
		for(j=0; j<(COLS+2*Ng); j++){
			for(q=0; q<Nq; q++){

				double x = getx(j);
				double y = gety(i);
				int jiq = Nq*(COLS+2*Ng)*i + Nq*j + q;

				if(y+x < 0.5){
					prim[jiq] = left[q];
				}else{
					prim[jiq] = right[q];
				
				}

			}
		}
	}
}

//Find fluxes in x direction
void fhll(double * Fhll, double * pL, double * pR, double alphM, double alphP, int i, int j){
	int q;
	double gamma = G;

	double fR[Nq], fL[Nq];
	double uR[Nq], uL[Nq];
	primToCons(pL, uL);
	primToCons(pR, uR);

	fR[MM_FLUX] = pR[RHO] * pR[VXX];
	fL[MM_FLUX] = pL[RHO] * pL[VXX];
	
	fR[PX_FLUX] = pR[RHO] * pR[VXX] * pR[VXX] + pR[PPP];
	fL[PX_FLUX] = pL[RHO] * pL[VXX] * pL[VXX] + pL[PPP];
	
	fR[PY_FLUX] = pR[RHO] * pR[VXX] * pR[VYY];
	fL[PY_FLUX] = pL[RHO] * pL[VXX] * pL[VYY];
	
	fR[EE_FLUX] = pR[VXX] * (pR[PPP]/(gamma-1.0) + 0.5 * pR[RHO] * (pR[VXX]*pR[VXX] + pR[VYY]*pR[VYY]) + pR[PPP]);
	fL[EE_FLUX] = pL[VXX] * (pL[PPP]/(gamma-1.0) + 0.5 * pL[RHO] * (pL[VXX]*pL[VXX] + pL[VYY]*pL[VYY]) + pL[PPP]);


	for(q=0; q<Nq; ++q){
		Fhll[Nq*(COLS + 2*Ng -1)*i + Nq*j + q] = (alphP*fL[q] + alphM*fR[q] - alphM*alphP * (uR[q] - uL[q])) / (alphP + alphM);
	}
	
}

void ghll(double * Ghll, double * pB, double * pA, double alphM, double alphP, int i, int j){
	int q;
	double gamma = G;

	double gB[Nq], gA[Nq];
	double uB[Nq], uA[Nq];
	primToCons(pB, uB);
	primToCons(pA, uA);
	
	gA[MM_FLUX] = pA[RHO] * pA[VYY];
	gB[MM_FLUX] = pB[RHO] * pB[VYY];
	
	gA[PX_FLUX] = pA[RHO] * pA[VXX] * pA[VYY];
	gB[PX_FLUX] = pB[RHO] * pB[VXX] * pB[VYY];
	
	gA[PY_FLUX] = pA[RHO] * pA[VYY] * pA[VYY] + pA[PPP];
	gB[PY_FLUX] = pB[RHO] * pB[VYY] * pB[VYY] + pB[PPP];
	
	gA[EE_FLUX] = pA[VYY] * (pA[PPP]/(gamma-1.0) + 0.5 * pA[RHO] * (pA[VXX]*pA[VXX] + pA[VYY]*pA[VYY]) + pA[PPP]);
	gB[EE_FLUX] = pB[VYY] * (pB[PPP]/(gamma-1.0) + 0.5 * pB[RHO] * (pB[VXX]*pB[VXX] + pB[VYY]*pB[VYY]) + pB[PPP]);

	for(q=0; q<Nq; ++q){
		Ghll[Nq*(ROWS + 2*Ng -1)*j + Nq*i + q] = (alphP*gB[q] + alphM*gA[q] - alphM*alphP * (uA[q] - uB[q])) / (alphP + alphM);
	}
		
}

double advance(double * prim, double dx, double dy, double dt){
	int i, j, q;
	
	int xBounds = (ROWS+2*Ng)*(COLS + 2*Ng -1);
	int yBounds = (COLS+2*Ng)*(ROWS + 2*Ng -1);
	double uL[Nq], uR[Nq];
	double pL[Nq], pR[Nq];
	double alphX[2  * xBounds];
	double alphY[2  * yBounds];
	double  Fiph[Nq * xBounds];
	double  Giph[Nq * yBounds];
	alpha(alphX, alphY, prim);
	double pLi[Nq*(ROWS + 2*Ng)];
	double pRi[Nq*(ROWS + 2*Ng)];
	double pAi[Nq*(COLS + 2*Ng)];
	double pBi[Nq*(COLS + 2*Ng)];

	///////////////////////////
/*
	//Constant boundary conditions
	//Store left and right side
	for(i=0; i<(ROWS+2*Ng); ++i){
		for(q=0; q<Nq; ++q){
			pLi[Nq*i + q] = prim[Nq*(COLS+2*Ng)*i + Nq*0               + q];
			pRi[Nq*i + q] = prim[Nq*(COLS+2*Ng)*i + Nq*(COLS+2*Ng - 1) + q];
		}
	}

	//Store top and bottom
	for(j=0; j<(COLS+2*Ng); ++j){
		for(q=0; q<Nq; ++q){
			pBi[Nq*j + q] = prim[Nq*(COLS+2*Ng)*0             + Nq*j + q];
			pAi[Nq*j + q] = prim[Nq*(COLS+2*Ng)*(ROWS+2*Ng-1) + Nq*j + q];
		}
	}
	/////////////////////////////////
*/	
	
	//Enforce Courant condition, dt < dx/(MAX alpha)
	double alphMax = -1e-4;

	for(i=0; i<(xBounds*2); ++i){
		if(alphX[i] > alphMax) alphMax = alphX[i];
	}

	if(dt > dx/alphMax) dt = 0.5*dx/alphMax;
	
	for(i=0; i<(yBounds*2); ++i){
		if(alphY[i] > alphMax) alphMax = alphY[i];
	}

	if(dt > dy/alphMax) dt = 0.5*dy/alphMax;
	////////////////////////////////////////////////
	
	//Find fluxes in X direction
	for(i=0; i<(ROWS+2*Ng); ++i){
		for(j=0; j<(COLS+2*Ng-1); ++j){
			for(q=0; q<Nq; ++q){
				int jiq = Nq*(COLS+2*Ng)*i + Nq*j + q;
				pL[q] = prim[jiq];
				pR[q] = prim[jiq + Nq];
			}
			fhll(Fiph, pL, pR, alphX[2*i*(COLS+2*Ng-1) + 2*j], alphX[2*i*(COLS+2*Ng-1) + 2*j+1], i, j);
		}
	}

	//Add stuff in X direction
	for(i=0; i<(ROWS+2*Ng); ++i){
		for(j=0; j<(COLS+2*Ng-1); ++j){
			for(q=0; q<Nq; ++q){
				int jiq = Nq*(COLS+2*Ng)*i + Nq*j + q;
				pL[q] = prim[jiq];
				pR[q] = prim[jiq + Nq];
			}
			primToCons(pL, uL);
			primToCons(pR, uR);
			
			for(q=0; q<Nq; ++q){
				uL[q] -= Fiph[Nq*(COLS+2*Ng-1)*i + Nq*j + q]*dt/dx;
				uR[q] += Fiph[Nq*(COLS+2*Ng-1)*i + Nq*j + q]*dt/dx;
			}

			consToPrim(uL, pL);
			consToPrim(uR, pR);
			
			for(q=0; q<Nq; ++q){
				int jiq = Nq*(COLS+2*Ng)*i + Nq*j + q;
				prim[jiq] = pL[q];
				prim[jiq + Nq] = pR[q];
			}
		}
	}

	//Find fluxes in Y direction
	for(j=0; j<(COLS+2*Ng); ++j){
		for(i=0; i<(ROWS+2*Ng-1); ++i){
			for(q=0; q<Nq; ++q){
				int jiq = Nq*(COLS+2*Ng)*i + Nq*j + q;
				pL[q] = prim[jiq];
				pR[q] = prim[jiq + Nq*(COLS+2*Ng)];
			}
			ghll(Giph, pL, pR, alphY[2*j*(ROWS+2*Ng-1) + 2*i], alphY[2*j*(ROWS+2*Ng-1) + 2*i+1], i, j);
		}
	}

	//Add stuff in Y direction
	for(j=0; j<(COLS+2*Ng); ++j){
		for(i=0; i<(ROWS+2*Ng-1); ++i){
			for(q=0; q<Nq; ++q){
				int jiq = Nq*(COLS+2*Ng)*i + Nq*j + q;
				pL[q] = prim[jiq];
				pR[q] = prim[jiq + Nq*(COLS+2*Ng)];
			}
			primToCons(pL, uL);
			primToCons(pR, uR);
			
			for(q=0; q<Nq; ++q){
				uL[q] -= Giph[Nq*(ROWS+2*Ng-1)*j + Nq*i + q]*dt/dy;
				uR[q] += Giph[Nq*(ROWS+2*Ng-1)*j + Nq*i + q]*dt/dy;
			}

			consToPrim(uL, pL);
			consToPrim(uR, pR);
			
			for(q=0; q<Nq; ++q){
				int jiq = Nq*(COLS+2*Ng)*i + Nq*j + q;
				prim[jiq] = pL[q];
				prim[jiq + Nq*(COLS+2*Ng)] = pR[q];
			}
		}
	}

	//Outflow in all boundaries
	
	//Boundary conditions for X
	for(i=0; i<(ROWS+2*Ng); ++i){
		for(q=0; q<Nq; ++q){
			//Left Side
			prim[Nq*(COLS+2*Ng)*i + Nq*0              + q] = prim[Nq*(COLS+2*Ng)*i + Nq*Ng       + q];
			prim[Nq*(COLS+2*Ng)*i + Nq*1              + q] = prim[Nq*(COLS+2*Ng)*i + Nq*Ng       + q];
			//Right Side
			prim[Nq*(COLS+2*Ng)*i + Nq*(COLS+2*Ng -1) + q] = prim[Nq*(COLS+2*Ng)*i + Nq*(COLS+1) + q];
			prim[Nq*(COLS+2*Ng)*i + Nq*(COLS+2*Ng -2) + q] = prim[Nq*(COLS+2*Ng)*i + Nq*(COLS+1) + q];

		}
	}

	//Boundary conditions for Y
	for(j=0; j<(COLS+2*Ng); ++j){
		for(q=0; q<Nq; ++q){
			
			//Bottom
			prim[Nq*(COLS+2*Ng)*0               + Nq*j + q] = prim[Nq*(COLS+2*Ng)*Ng     + Nq*j + q];
			prim[Nq*(COLS+2*Ng)*1               + Nq*j + q] = prim[Nq*(COLS+2*Ng)*Ng     + Nq*j + q];
			
			//Top
			prim[Nq*(COLS+2*Ng)*(ROWS+2*Ng - 2) + Nq*j + q] = prim[Nq*(COLS+2*Ng)*(ROWS+1) + Nq*j + q];
			prim[Nq*(COLS+2*Ng)*(ROWS+2*Ng - 1) + Nq*j + q] = prim[Nq*(COLS+2*Ng)*(ROWS+1) + Nq*j + q];

		}
	}

/*
	//Reinput boundary conditions in X
	for(i=0; i<(ROWS+2*Ng); ++i){
		for(q=0; q<Nq; ++q){
			prim[Nq*(COLS+2*Ng)*i + Nq*0               + q] = pLi[Nq*i + q];
			prim[Nq*(COLS+2*Ng)*i + Nq*(COLS+2*Ng - 1) + q] = pRi[Nq*i + q];
		}
	}

	//Reinput boundary conditions in Y
	for(j=0; j<(COLS+2*Ng); ++j){
		for(q=0; q<Nq; ++q){
			prim[Nq*(COLS+2*Ng)*0             + Nq*j + q] = pBi[Nq*j + q];
			prim[Nq*(COLS+2*Ng)*(ROWS+2*Ng-1) + Nq*j + q] = pAi[Nq*j + q];
		}
	}
*/
	///////////////////////////////////
	
	return(dt);
}

int main(void){
	FILE *xShock, *yShock, *inFlow;

//	xShock = fopen("shockX.txt", "w");
//	yShock = fopen("shockY.txt", "w");
	inFlow = fopen("flowIn.txt", "w");

	int i, j, q, count;
	count = 0;
	double prim[Nq*(ROWS+2*Ng)*(COLS+2*Ng)];
	set_initial(prim);
/*
	for(i=0; i<(ROWS+2*Ng); ++i){
		for(j=0; j<(COLS+2*Ng); ++j){
			printf("Zone(%d %d): P = %.2e, rho = %.2e, vx = %.2e, vy = %.2e\n",
				       	j, i, prim[Nq*(COLS*i+j) + 0], prim[Nq*(COLS*i+j) + 1], prim[Nq*(COLS*i+j) + 2], prim[Nq*(COLS*i+j) + 3]);
		}
	}
*/
	double dx, dy, t, dt;
        dx = L/(double)COLS;
	dy = H/(double)ROWS;
	t = 0.0;
	dt = T;

/*
	dt = advance(prim, dx, dy, dt);
	printf("%d dt = %e\n", count++, dt);
	dt = advance(prim, dx, dy, dt);
	printf("%d dt = %e\n", count++, dt);
	dt = advance(prim, dx, dy, dt);
	printf("%d dt = %e\n", count++, dt);
	dt = advance(prim, dx, dy, dt);
	printf("%d dt = %e\n", count++ ,dt);
	dt = advance(prim, dx, dy, dt);
	printf("%d dt = %e\n", count++ ,dt);
	dt = advance(prim, dx, dy, dt);
	printf("%d dt = %e\n", count++ ,dt);
	for(i=0; i<(ROWS+2*Ng); ++i){
		for(j=0; j<(COLS+2*Ng); ++j){
			printf("Zone(%d %d): P = %.2e, rho = %.2e, vx = %.2e, vy = %.2e\n",
				       	j, i, prim[Nq*(COLS+2*Ng)*i + Nq*j + 0], prim[Nq*(COLS+2*Ng)*i + Nq*j + 1], prim[Nq*(COLS+2*Ng)*i + Nq*j + 2], prim[Nq*(COLS+2*Ng)*i + Nq*j + 3]);
		}
	}

*/
	while(t<T){
		dt = advance(prim, dx, dy, dt);
		t+=dt;
		printf("%d dt = %e\n", count++ ,dt);
		if (dt<1e-20) break;
	}
	
	for(i=0; i<(ROWS+2*Ng); ++i){
		for(j=0; j<(COLS+2*Ng); ++j){
			fprintf(inFlow, "%e %e %e\n",
				       	prim[Nq*(COLS+2*Ng)*i+Nq*j + 1], getx(j), gety(i));
		}
	}

//	fclose(xShock);
//	fclose(yShock);
	fclose(inFlow);

	return(0);
}
