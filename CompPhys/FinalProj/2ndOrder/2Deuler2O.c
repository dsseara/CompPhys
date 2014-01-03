#include <stdio.h>
#include <math.h>

double G = 1.4;  //Adiabatic heat index

double L = 0.75 ;  //Length of box
double H = 1.5;  //Height of box
int COLS = 200;  //Columns
int ROWS = 200;  //Rows
int   Nq = 4;   //Number of conserved quantities
int   Ng = 2;   //Number of ghost zones
double T = 0.0002; //Total time integrated over
double gravity = -0.2; //acceleration due to gravity
double theta = 1.5;  //theta for minmod

enum{PPP, RHO, VXX, VYY}; //pressure, mass density, velocity in x, velocity in y
enum{DEN, MOX, MOY, EEE}; //mass density, momentum in x, momentum in y, energy
enum{XXX, YYY};
enum{MM_FLUX, PX_FLUX, PY_FLUX, EE_FLUX};

double getx (int j){
	return( L*( ( (double)(j-Ng) + 0.5) / (double)COLS ));
}

double gety (int i){
	return( H*( ( (double)(i-Ng) + 0.5) / (double)ROWS ));
}

double minmod(double x, double y, double z){
	
	double xsgn, ysgn, zsgn;
	double min1, min2;

	if(x<0){
		xsgn = -1.0;
	}else if(x>0){
		xsgn = 1.0;
	}else{
		xsgn = 0.0;
	}

	if(y<0){
		ysgn = -1.0;
	}else if(y>0){
		ysgn = 1.0;
	}else{
		ysgn = 0.0;
	}
	
	if(z<0){
		zsgn = -1.0;
	}else if(z>0){
		zsgn = 1.0;
	}else{
		zsgn = 0.0;
	}
	min1 = fmin(fabs(x), fabs(y));
	min2 = fmin(fabs(min1), fabs(z));

	return(0.25*fabs(xsgn + ysgn)*(xsgn+zsgn)*min2);
}

//given a U vector in zone i, get the corresponding
//primitive variables. prim = (p, rho, vx, vy)
void consToPrim(double * u, double * prim){

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
void primToCons(double * prim, double * u){
	
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
	Energy  = P/(gamma - 1.0) + 0.5 * rho * (vx*vx + vy*vy);

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
	double vL, vR; //velocities left, right
	double slopeL, slopeR;
	
	//X Direction
	for(i=0; i<(ROWS+2*Ng); ++i){
		for(j=1; j<(COLS+Ng); ++j){
			//Get left and right cells at boundarys
			for(q=0; q<Nq; ++q){
				int jiq = Nq*(COLS+2*Ng)*i+ Nq*j +q; //Gets to cell (j,i), and 
				
				slopeL = minmod(theta*(prim[jiq]      - prim[jiq - Nq]),
					       	  0.5*(prim[jiq + Nq] - prim[jiq - Nq]),
					       	theta*(prim[jiq + Nq] - prim[jiq]));

				slopeR = minmod(theta*(prim[jiq +   Nq] - prim[jiq]),
					       	  0.5*(prim[jiq + 2*Nq] - prim[jiq]),
					       	theta*(prim[jiq + 2*Nq] - prim[jiq+Nq]));

				
				pL[q] = prim[jiq]      + 0.5*slopeL; //quantity q in that cell
				pR[q] = prim[jiq + Nq] - 0.5*slopeR; //need Nq(COLS*i + j+1) + q
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
			
			alphX[2*i*(COLS+1) + 2*(j-1)  ] = maxM;
			alphX[2*i*(COLS+1) + 2*(j-1)+1] = maxP;
		}
	}

	//Y Direction
	for(j=0; j<(COLS+2*Ng); ++j){
		for(i=1; i<(ROWS+Ng); ++i){
			//Get left and right cells at boundays
			for(q=0; q<Nq; ++q){
				int jiq = Nq*(COLS+2*Ng)*i + Nq*j + q; //Gets to cell (j,i), and 
				
				slopeL = minmod(theta*(prim[jiq]                  - prim[jiq - (COLS+2*Ng)*Nq] ),
					       	0.5  *(prim[jiq + (COLS+2*Ng)*Nq] - prim[jiq - (COLS+2*Ng)*Nq] ),
					       	theta*(prim[jiq + (COLS+2*Ng)*Nq] - prim[jiq])                 );

				slopeR = minmod(theta*(prim[jiq +   (COLS+2*Ng)*Nq] - prim[jiq]                  ),
					       	0.5  *(prim[jiq + 2*(COLS+2*Ng)*Nq] - prim[jiq]                  ),
					       	theta*(prim[jiq + 2*(COLS+2*Ng)*Nq] - prim[jiq + (COLS+2*Ng)*Nq]));


				
				pL[q] = prim[jiq]                  + 0.5*slopeL; //quantity q in that cell
				pR[q] = prim[jiq + Nq*(COLS+2*Ng)] - 0.5*slopeR; //need Nq(COLS*(i+1) + j) + q
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
			
			alphY[2*j*(ROWS+1) + 2*(i-1)  ] = maxM;
			alphY[2*j*(ROWS+1) + 2*(i-1)+1] = maxP;

		}
	}
}


void set_initial(double * prim){
	int i, j, q;
	double A     = 0.100; //amplitude of perturbation
	double p0    =  5.0; //arbitrary initial pressure
	double sigma =  0.2; //verticle extent of velocity perturbation

	/*Initial top half of fluid****/
	double top[4] = {0.0, 5.0, 0.0, 0.0};
	/*Initial bottom half of tube***/
	double bot[4] = {0.0, 1.0, 0.0, 0.0};
	
	for(i =0 ; i<(ROWS+2*Ng) ; ++i){
		for(j=0; j<(COLS+2*Ng); ++j){
			for(q=0; q<Nq; ++q){

				double x = getx(j);
				double y = gety(i);

				top[VYY] = A*(1 + sin(2*M_PI*x / L))*exp(-(y - 0.5*H)*(y-0.5*H)/(sigma*sigma));
				top[PPP] = p0 + bot[RHO]*gravity*0.5*H + top[RHO]*gravity*(y-0.5*H);


				bot[VYY] = A*(1 + sin(2*M_PI*x / L))*exp(-(y - 0.5*H)*(y-0.5*H)/(sigma*sigma));
				bot[PPP] = p0 + bot[RHO] * gravity * y;

				int jiq = Nq*(COLS+2*Ng)*i + Nq*j + q;

			//	if(y < (0.5*H + A*sin(M_PI*x/L))){
				if(y < 0.5*H){
					prim[jiq] = bot[q];
				}else{
					prim[jiq] = top[q];
				
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
		Fhll[Nq*(COLS + 1)*i + Nq*(j-1) + q] = (alphP*fL[q] + alphM*fR[q] - alphM*alphP * (uR[q] - uL[q])) / (alphP + alphM);
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
		Ghll[Nq*(ROWS + 1)*j + Nq*(i-1) + q] = (alphP*gB[q] + alphM*gA[q] - alphM*alphP * (uA[q] - uB[q])) / (alphP + alphM);
	}
		
}

double step1(double * prim0, double dx, double dy, double dt){
	int i, j, q;
	
	int xBounds = (ROWS+2*Ng)*(COLS+1);
	int yBounds = (COLS+2*Ng)*(ROWS+1);
	double uL[Nq], uR[Nq];
	double pL[Nq], pR[Nq];
	double alphX[2  * xBounds];
	double alphY[2  * yBounds];
	double  Fiph[Nq * xBounds];
	double  Giph[Nq * yBounds];

	alpha(alphX, alphY, prim0);
	

	double slopeL, slopeR;

	//Enforce Courant condition, dt < dx/(MAX alpha)
	double alphMax = -1e-4;

	for(i=0; i<(xBounds*2); ++i){
		if(alphX[i] > alphMax) alphMax = alphX[i];
	}

	if(dt > dx/alphMax) dt = 0.5*dx/alphMax;
	
	for(i=0; i<(yBounds*2); ++i){
		if(alphY[i] > alphMax) alphMax = alphY[i];
	}

	if(dt >     dy/alphMax) dt = 0.5*dy/alphMax;
//	if(dt < 0.1*dy/alphMax) dt = 0.5*dy/alphMax;
	////////////////////////////////////////////////

	//Find fluxes in X direction
	for(i=0; i<(ROWS+2*Ng); ++i){
		for(j=1; j<(COLS+Ng); ++j){
			for(q=0; q<Nq; ++q){
				int jiq = Nq*(COLS+2*Ng)*i + Nq*j + q;

				slopeL = minmod(theta*(prim0[jiq]      - prim0[jiq - Nq]),
					       	  0.5*(prim0[jiq + Nq] - prim0[jiq - Nq]),
					       	theta*(prim0[jiq + Nq] - prim0[jiq])    );

				slopeR = minmod(theta*(prim0[jiq +   Nq] - prim0[jiq]     ),
					       	  0.5*(prim0[jiq + 2*Nq] - prim0[jiq]     ),
					       	theta*(prim0[jiq + 2*Nq] - prim0[jiq+Nq]) );

				pL[q] = prim0[jiq]      + 0.5*slopeL;
				pR[q] = prim0[jiq + Nq] - 0.5*slopeR;
			}
			fhll(Fiph, pL, pR, alphX[2*i*(COLS+1) + 2*(j-1)], alphX[2*i*(COLS+1) + 2*(j-1)+1], i, j);
		}
	}

	//Add stuff in X direction
	for(i=0; i<(ROWS+2*Ng); ++i){
		for(j=1; j<(COLS+Ng); ++j){
			for(q=0; q<Nq; ++q){
				int jiq = Nq*(COLS+2*Ng)*i + Nq*j + q;
			
				slopeL = minmod(theta*(prim0[jiq]      - prim0[jiq - Nq]), 
						  0.5*(prim0[jiq + Nq] - prim0[jiq - Nq]),
					       	theta*(prim0[jiq + Nq] - prim0[jiq])    );

				slopeR = minmod(theta*(prim0[jiq +   Nq] - prim0[jiq]      ),
					          0.5*(prim0[jiq + 2*Nq] - prim0[jiq]      ),
					       	theta*(prim0[jiq + 2*Nq] - prim0[jiq + Nq]));

				pL[q] = prim0[jiq]      + 0.5*slopeL;
				pR[q] = prim0[jiq + Nq] - 0.5*slopeR;
			}

			primToCons(pL, uL);
			primToCons(pR, uR);
			
			for(q=0; q<Nq; ++q){
				uL[q] -= Fiph[Nq*(COLS+1)*i + Nq*(j-1) + q]*dt/dx;
				uR[q] += Fiph[Nq*(COLS+1)*i + Nq*(j-1) + q]*dt/dx;
			}

			consToPrim(uL, pL);
			consToPrim(uR, pR);
			
			for(q=0; q<Nq; ++q){
				int jiq = Nq*(COLS+2*Ng)*i + Nq*j + q;
				prim0[jiq]      = pL[q];
				prim0[jiq + Nq] = pR[q];
			}
		}
	}

	//Find fluxes in Y direction
	//with gravitational source term
	double sourceL[4] = {0.0, 0.0, 0.0, 0.0};
	double sourceR[4] = {0.0, 0.0, 0.0, 0.0};

	for(j=0; j<(COLS+2*Ng); ++j){
		for(i=1; i<(ROWS+Ng); ++i){
			for(q=0; q<Nq; ++q){
				int jiq = Nq*(COLS+2*Ng)*i + Nq*j + q;
				
				slopeL = minmod(theta*(prim0[jiq]                    - prim0[jiq - (COLS+2*Ng)*Nq] ),
					       	0.5  *(prim0[jiq +   (COLS+2*Ng)*Nq] - prim0[jiq - (COLS+2*Ng)*Nq] ),
					       	theta*(prim0[jiq +   (COLS+2*Ng)*Nq] - prim0[jiq])                 );

				slopeR = minmod(theta*(prim0[jiq +   (COLS+2*Ng)*Nq] - prim0[jiq]                  ),
					       	0.5  *(prim0[jiq + 2*(COLS+2*Ng)*Nq] - prim0[jiq]                  ),
					       	theta*(prim0[jiq + 2*(COLS+2*Ng)*Nq] - prim0[jiq + (COLS+2*Ng)*Nq]));

				pL[q] = prim0[jiq]                  + 0.5*slopeL;
				pR[q] = prim0[jiq + (COLS+2*Ng)*Nq] - 0.5*slopeR;
			}
			
			ghll(Giph, pL, pR, alphY[2*j*(ROWS+1) + 2*(i-1)], alphY[2*j*(ROWS+1) + 2*(i-1)+1], i, j);
		}
	}

	//Add stuff in Y direction
	for(j=0; j<(COLS+2*Ng); ++j){
		for(i=1; i<(ROWS+Ng); ++i){
			for(q=0; q<Nq; ++q){
				int jiq = Nq*(COLS+2*Ng)*i + Nq*j + q;
			
				slopeL = minmod(theta*(prim0[jiq]                  - prim0[jiq - (COLS+2*Ng)*Nq]),
					       	  0.5*(prim0[jiq + (COLS+2*Ng)*Nq] - prim0[jiq - (COLS+2*Ng)*Nq]),
					       	theta*(prim0[jiq + (COLS+2*Ng)*Nq] - prim0[jiq])                );

				slopeR = minmod(theta*(prim0[jiq +   (COLS+2*Ng) * Nq] - prim0[jiq]),
					       	0.5  *(prim0[jiq + 2*(COLS+2*Ng) * Nq] - prim0[jiq]),
					       	theta*(prim0[jiq + 2*(COLS+2*Ng) * Nq] - prim0[jiq+(COLS+2*Ng)*Nq]));

				pL[q] = prim0[jiq]                    + 0.5*slopeL;
				pR[q] = prim0[jiq + (COLS + 2*Ng)*Nq] - 0.5*slopeR;
			}

			sourceL[2] = pL[RHO]*gravity;
			sourceL[3] = pL[RHO]*pL[VYY]*gravity;
			sourceR[2] = pR[RHO]*gravity;
			sourceR[3] = pR[RHO]*pR[VYY]*gravity;

			primToCons(pL, uL);
			primToCons(pR, uR);
			
			for(q=0; q<Nq; ++q){
				uL[q] = uL[q] - Giph[Nq*(ROWS+1)*j + Nq*(i-1) + q]*dt/dy + dt*sourceL[q];
				uR[q] = uR[q] + Giph[Nq*(ROWS+1)*j + Nq*(i-1) + q]*dt/dy + dt*sourceR[q] ;
			}

			consToPrim(uL, pL);
			consToPrim(uR, pR);
			
			for(q=0; q<Nq; ++q){
				int jiq = Nq*(COLS+2*Ng)*i + Nq*j + q;
				prim0[jiq]                  = pL[q];
				prim0[jiq + Nq*(COLS+2*Ng)] = pR[q];
			}
		}
	}

	////////////////////////////////////
	//Boundary conditions for X
	//Periodic
	for(i=0; i<(ROWS+2*Ng); ++i){
		for(q=0; q<Nq; ++q){
			//Left Side
			prim0[Nq*(COLS+2*Ng)*i + Nq*0           + q] = prim0[Nq*(COLS+2*Ng)*i + Nq*(COLS+0) + q];
			prim0[Nq*(COLS+2*Ng)*i + Nq*1           + q] = prim0[Nq*(COLS+2*Ng)*i + Nq*(COLS+1) + q];
			//Right Side
			prim0[Nq*(COLS+2*Ng)*i + Nq*(COLS+Ng+1) + q] = prim0[Nq*(COLS+2*Ng)*i + Nq*(Ng+1)   + q];
			prim0[Nq*(COLS+2*Ng)*i + Nq*(COLS+Ng+0) + q] = prim0[Nq*(COLS+2*Ng)*i + Nq*(Ng+0)   + q];

		}
	}

	//Boundary conditions for Y
	//Reflective
	for(j=0; j<(COLS+2*Ng); ++j){
		for(q=0; q<Nq; ++q){

			if(q == 3){
				//Bottom
				prim0[Nq*(COLS+2*Ng)*0 + Nq*j + q] = -prim0[Nq*(COLS+2*Ng)*Ng + Nq*j + q];
				prim0[Nq*(COLS+2*Ng)*1 + Nq*j + q] = -prim0[Nq*(COLS+2*Ng)*Ng + Nq*j + q];
				//Top
				prim0[Nq*(COLS+2*Ng)*(ROWS+2*Ng - 2) + Nq*j + q] = -prim0[Nq*(COLS+2*Ng)*(ROWS+1) + Nq*j + q];
				prim0[Nq*(COLS+2*Ng)*(ROWS+2*Ng - 1) + Nq*j + q] = -prim0[Nq*(COLS+2*Ng)*(ROWS+1) + Nq*j + q];
			}else{
				//Bottom
				prim0[Nq*(COLS+2*Ng)*0 + Nq*j + q] = prim0[Nq*(COLS+2*Ng)*Ng + Nq*j + q];
				prim0[Nq*(COLS+2*Ng)*1 + Nq*j + q] = prim0[Nq*(COLS+2*Ng)*Ng + Nq*j + q];
				//Top
				prim0[Nq*(COLS+2*Ng)*(ROWS+2*Ng - 2) + Nq*j + q] = prim0[Nq*(COLS+2*Ng)*(ROWS+1) + Nq*j + q];
				prim0[Nq*(COLS+2*Ng)*(ROWS+2*Ng - 1) + Nq*j + q] = prim0[Nq*(COLS+2*Ng)*(ROWS+1) + Nq*j + q];

			}
		}
	}

	///////////////////////////////////
	return(dt);
}
/*
void step2(double * prim0, double * prim1, double * prim2, double dx, double dy, double dt){
	int i, j, q;
	
	int xBounds = (ROWS+2*Ng)*(COLS + 2*Ng -1);
	int yBounds = (COLS+2*Ng)*(ROWS + 2*Ng -1);
	double uL[Nq], uR[Nq];
	double pL[Nq], pR[Nq];
	double alphX[2  * xBounds];
	double alphY[2  * yBounds];
	double  Fiph[Nq * xBounds];
	double  Giph[Nq * yBounds];
	alpha(alphX, alphY, prim0);
	double slopeL, slopeR;

	//Find fluxes in X direction
	for(i=0; i<(ROWS+2*Ng); ++i){
		for(j=0; j<(COLS+2*Ng-1); ++j){
			for(q=0; q<Nq; ++q){
				int jiq = Nq*(COLS+2*Ng)*i + Nq*j + q;

				slopeL = minmod(theta*(prim1[jiq]      - prim1[jiq - Nq]), 
						  0.5*(prim1[jiq + Nq] - prim1[jiq - Nq]), 
						theta*(prim1[jiq + Nq] - prim1[jiq]));
				
				slopeR = minmod(theta*(prim1[jiq +   Nq] - prim1[jiq]),
					       	  0.5*(prim1[jiq + 2*Nq] - prim1[jiq]),
					       	theta*(prim1[jiq + 2*Nq] - prim1[jiq+Nq]));

				pL[q] = prim1[jiq]      + 0.5*slopeL;
				pR[q] = prim1[jiq + Nq] - 0.5*slopeR;
			}
			fhll(Fiph, pL, pR, alphX[2*i*(COLS+2*Ng-1) + 2*j], alphX[2*i*(COLS+2*Ng-1) + 2*j+1], i, j);
		}
	}

	//Add stuff in X direction
	for(i=0; i<(ROWS+2*Ng); ++i){
		for(j=0; j<(COLS+2*Ng-1); ++j){
			for(q=0; q<Nq; ++q){
				int jiq = Nq*(COLS+2*Ng)*i + Nq*j + q;
			
				slopeL = minmod(theta*(prim1[jiq]      - prim1[jiq - Nq]), 
						  0.5*(prim1[jiq + Nq] - prim1[jiq - Nq]), 
						theta*(prim1[jiq + Nq] - prim1[jiq]));
				
				slopeR = minmod(theta*(prim1[jiq +   Nq] - prim1[jiq]),
					       	  0.5*(prim1[jiq + 2*Nq] - prim1[jiq]),
					       	theta*(prim1[jiq + 2*Nq] - prim1[jiq+Nq]));


				pL[q] = prim1[jiq]      + 0.5*slopeL;
				pR[q] = prim1[jiq + Nq] - 0.5*slopeR;
			}

			primToCons(pL, uL);
			primToCons(pR, uR);
			
			for(q=0; q<Nq; ++q){
				uL[q] -= 0.25 * Fiph[Nq*(COLS+2*Ng-1)*i + Nq*j + q]*dt/dx;
				uR[q] += 0.25 * Fiph[Nq*(COLS+2*Ng-1)*i + Nq*j + q]*dt/dx;
			}

			consToPrim(uL, pL);
			consToPrim(uR, pR);
			
			for(q=0; q<Nq; ++q){
				int jiq = Nq*(COLS+2*Ng)*i + Nq*j + q;
				prim2[jiq]      = pL[q];
				prim2[jiq + Nq] = pR[q];
			}
		}
	}

	//Find fluxes in Y direction
	//with gravitational source term
	double sourceL[4] = {0.0, 0.0, 0.0, 0.0};
	double sourceR[4] = {0.0, 0.0, 0.0, 0.0};

	for(j=0; j<(COLS+2*Ng); ++j){
		for(i=0; i<(ROWS+2*Ng-1); ++i){
			for(q=0; q<Nq; ++q){
				int jiq = Nq*(COLS+2*Ng)*i + Nq*j + q;
				
				slopeL = minmod(theta*(prim1[jiq]                      - prim1[jiq - (COLS+2*Ng)*Nq]),
					       	0.5  *(prim1[jiq +   (COLS+2*Ng) * Nq] - prim1[jiq - (COLS+2*Ng)*Nq]),
					       	theta*(prim1[jiq +   (COLS+2*Ng) * Nq] - prim1[jiq])                );

				slopeR = minmod(theta*(prim1[jiq +   (COLS+2*Ng) * Nq] - prim1[jiq]                  ),
					       	0.5  *(prim1[jiq + 2*(COLS+2*Ng) * Nq] - prim1[jiq]                  ),
					       	theta*(prim1[jiq + 2*(COLS+2*Ng) * Nq] - prim1[jiq + (COLS+2*Ng)*Nq]));

				pL[q] = prim1[jiq]                  + 0.5*slopeL;
				pR[q] = prim1[jiq + (COLS+2*Ng)*Nq] - 0.5*slopeR;
			}
			
			ghll(Giph, pL, pR, alphY[2*j*(ROWS+2*Ng-1) + 2*i], alphY[2*j*(ROWS+2*Ng-1) + 2*i+1], i, j);
		}
	}

	//Add stuff in Y direction
	for(j=0; j<(COLS+2*Ng); ++j){
		for(i=0; i<(ROWS+2*Ng-1); ++i){
			for(q=0; q<Nq; ++q){
				int jiq = Nq*(COLS+2*Ng)*i + Nq*j + q;
			
				slopeL = minmod(theta*(prim1[jiq]                      - prim1[jiq - (COLS+2*Ng) * Nq]),
					          0.5*(prim1[jiq +   (COLS+2*Ng) * Nq] - prim1[jiq - (COLS+2*Ng) * Nq]),
					       	theta*(prim1[jiq +   (COLS+2*Ng) * Nq] - prim1[jiq])                  );

				slopeR = minmod(theta*(prim1[jiq +   (COLS+2*Ng) * Nq] - prim1[jiq]                   ),
					       	0.5  *(prim1[jiq + 2*(COLS+2*Ng) * Nq] - prim1[jiq]                   ),
					       	theta*(prim1[jiq + 2*(COLS+2*Ng) * Nq] - prim1[jiq + (COLS+2*Ng) * Nq]));

				pL[q] = prim1[jiq]                    + 0.5*slopeL;
				pR[q] = prim1[jiq + (COLS+2*Ng) * Nq] - 0.5*slopeR;
			}

			sourceL[2] = pL[RHO]*gravity;
			sourceL[3] = pL[RHO]*pL[VYY]*gravity;
			sourceR[2] = pR[RHO]*gravity;
			sourceR[3] = pR[RHO]*pR[VYY]*gravity;

			primToCons(pL, uL);
			primToCons(pR, uR);
			
			for(q=0; q<Nq; ++q){
				uL[q] -= (Giph[Nq*(ROWS+2*Ng-1)*j + Nq*i + q]*dt/dy + dt*sourceL[q]);
				uR[q] +=  Giph[Nq*(ROWS+2*Ng-1)*j + Nq*i + q]*dt/dy + dt*sourceR[q] ;
			}

			consToPrim(uL, pL);
			consToPrim(uR, pR);
			
			for(q=0; q<Nq; ++q){
				int jiq = Nq*(COLS+2*Ng)*i + Nq*j + q;
				prim2[jiq]                  = pL[q];
				prim2[jiq + Nq*(COLS+2*Ng)] = pR[q];
			}
		}
	}


	////////////////////////////////////
	//Boundary conditions for X
	//Periodic
	for(i=0; i<(ROWS+2*Ng); ++i){
		for(q=0; q<Nq; ++q){
			//Left Side
			prim[Nq*(COLS+2*Ng)*i + Nq*0           + q] = prim[Nq*(COLS+2*Ng)*i + Nq*(COLS+0) + q];
			prim[Nq*(COLS+2*Ng)*i + Nq*1           + q] = prim[Nq*(COLS+2*Ng)*i + Nq*(COLS+1) + q];
			//Right Side
			prim[Nq*(COLS+2*Ng)*i + Nq*(COLS+Ng+1) + q] = prim[Nq*(COLS+2*Ng)*i + Nq*(Ng+1)   + q];
			prim[Nq*(COLS+2*Ng)*i + Nq*(COLS+Ng+0) + q] = prim[Nq*(COLS+2*Ng)*i + Nq*(Ng+0)   + q];

		}
	}
	//Boundary conditions for Y
	//Reflective
	for(j=0; j<(COLS+2*Ng); ++j){
		for(q=0; q<Nq; ++q){

			if(q == 3){
				//Bottom
				prim[Nq*(COLS+2*Ng)*0 + Nq*j + q] = -prim[Nq*(COLS+2*Ng)*Ng + Nq*j + q];
				prim[Nq*(COLS+2*Ng)*1 + Nq*j + q] = -prim[Nq*(COLS+2*Ng)*Ng + Nq*j + q];
				//Top
				prim[Nq*(COLS+2*Ng)*(ROWS+2*Ng - 2) + Nq*j + VYY] = -prim[Nq*(COLS+2*Ng)*(ROWS+1) + Nq*j + VYY];
				prim[Nq*(COLS+2*Ng)*(ROWS+2*Ng - 1) + Nq*j + VYY] = -prim[Nq*(COLS+2*Ng)*(ROWS+1) + Nq*j + VYY];
			}else{
				//Bottom
				prim[Nq*(COLS+2*Ng)*0 + Nq*j + q] = prim[Nq*(COLS+2*Ng)*Ng + Nq*j + q];
				prim[Nq*(COLS+2*Ng)*1 + Nq*j + q] = prim[Nq*(COLS+2*Ng)*Ng + Nq*j + q];
				//Top
				prim[Nq*(COLS+2*Ng)*(ROWS+2*Ng - 2) + Nq*j + VYY] = prim[Nq*(COLS+2*Ng)*(ROWS+1) + Nq*j + VYY];
				prim[Nq*(COLS+2*Ng)*(ROWS+2*Ng - 1) + Nq*j + VYY] = prim[Nq*(COLS+2*Ng)*(ROWS+1) + Nq*j + VYY];

			}
		}
	}

	///////////////////////////////////
}

void advance(double * prim0, double * prim1, double * prim2, double dx, double dy, double dt){
	int i, j, q;
	
	int xBounds = (ROWS+2*Ng)*(COLS + 2*Ng -1);
	int yBounds = (COLS+2*Ng)*(ROWS + 2*Ng -1);
	double uL[Nq], uR[Nq];
	double pL[Nq], pR[Nq];
	double alphX[2  * xBounds];
	double alphY[2  * yBounds];
	double  Fiph[Nq * xBounds];
	double  Giph[Nq * yBounds];
	alpha(alphX, alphY, prim0);
	double slopeL, slopeR;

	//Enforce Courant condition, dt < dx/(MAX alpha)
	double alphMax = -1e-4;

	for(i=0; i<(xBounds*2); ++i){
		if(alphX[i] > alphMax) alphMax = alphX[i];
	}

	if(dt > dx/alphMax) dt = 0.5*dx/alphMax;
	
	for(i=0; i<(yBounds*2); ++i){
		if(alphY[i] > alphMax) alphMax = alphY[i];
	}

	if(dt >     dy/alphMax) dt = 0.8*dy/alphMax;
	if(dt < 0.1*dy/alphMax) dt = 0.5*dy/alphMax;
	////////////////////////////////////////////////
	
	//Find fluxes in X direction
	for(i=0; i<(ROWS+2*Ng); ++i){
		for(j=0; j<(COLS+2*Ng-1); ++j){
			for(q=0; q<Nq; ++q){
				int jiq = Nq*(COLS+2*Ng)*i + Nq*j + q;

				slopeL = minmod(theta*(prim1[jiq]      - prim1[jiq - Nq]), 
						  0.5*(prim1[jiq + Nq] - prim1[jiq - Nq]), 
						theta*(prim1[jiq + Nq] - prim1[jiq]));
				
				slopeR = minmod(theta*(prim1[jiq +   Nq] - prim1[jiq]),
					       	  0.5*(prim1[jiq + 2*Nq] - prim1[jiq]),
					       	theta*(prim1[jiq + 2*Nq] - prim1[jiq+Nq]));

				pL[q] = prim1[jiq]      + 0.5*slopeL;
				pR[q] = prim1[jiq + Nq] - 0.5*slopeR;
			}
			fhll(Fiph, pL, pR, alphX[2*i*(COLS+2*Ng-1) + 2*j], alphX[2*i*(COLS+2*Ng-1) + 2*j+1], i, j);
		}
	}

	//Add stuff in X direction
	for(i=0; i<(ROWS+2*Ng); ++i){
		for(j=0; j<(COLS+2*Ng-1); ++j){
			for(q=0; q<Nq; ++q){
				int jiq = Nq*(COLS+2*Ng)*i + Nq*j + q;
			
				slopeL = minmod(theta*(prim1[jiq]      - prim1[jiq - Nq]), 
						  0.5*(prim1[jiq + Nq] - prim1[jiq - Nq]), 
						theta*(prim1[jiq + Nq] - prim1[jiq]));
				
				slopeR = minmod(theta*(prim1[jiq +   Nq] - prim1[jiq]),
					       	  0.5*(prim1[jiq + 2*Nq] - prim1[jiq]),
					       	theta*(prim1[jiq + 2*Nq] - prim1[jiq+Nq]));


				pL[q] = prim1[jiq]      + 0.5*slopeL;
				pR[q] = prim1[jiq + Nq] - 0.5*slopeR;
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
				prim2[jiq]      = pL[q];
				prim2[jiq + Nq] = pR[q];
			}
		}
	}

	//Find fluxes in Y direction
	//with gravitational source term
	double sourceL[4] = {0.0, 0.0, 0.0, 0.0};
	double sourceR[4] = {0.0, 0.0, 0.0, 0.0};

	for(j=0; j<(COLS+2*Ng); ++j){
		for(i=0; i<(ROWS+2*Ng-1); ++i){
			for(q=0; q<Nq; ++q){
				int jiq = Nq*(COLS+2*Ng)*i + Nq*j + q;
				
				slopeL = minmod(theta*(prim1[jiq]                      - prim1[jiq - (COLS+2*Ng)*Nq]),
					       	0.5  *(prim1[jiq +   (COLS+2*Ng) * Nq] - prim1[jiq - (COLS+2*Ng)*Nq]),
					       	theta*(prim1[jiq +   (COLS+2*Ng) * Nq] - prim1[jiq])                );

				slopeR = minmod(theta*(prim1[jiq +   (COLS+2*Ng) * Nq] - prim1[jiq]                  ),
					       	0.5  *(prim1[jiq + 2*(COLS+2*Ng) * Nq] - prim1[jiq]                  ),
					       	theta*(prim1[jiq + 2*(COLS+2*Ng) * Nq] - prim1[jiq + (COLS+2*Ng)*Nq]));

				pL[q] = prim1[jiq]                  + 0.5*slopeL;
				pR[q] = prim1[jiq + (COLS+2*Ng)*Nq] - 0.5*slopeR;
			}
			
			ghll(Giph, pL, pR, alphY[2*j*(ROWS+2*Ng-1) + 2*i], alphY[2*j*(ROWS+2*Ng-1) + 2*i+1], i, j);
		}
	}

	//Add stuff in Y direction
	for(j=0; j<(COLS+2*Ng); ++j){
		for(i=0; i<(ROWS+2*Ng-1); ++i){
			for(q=0; q<Nq; ++q){
				int jiq = Nq*(COLS+2*Ng)*i + Nq*j + q;
			
				slopeL = minmod(theta*(prim1[jiq]                      - prim1[jiq - (COLS+2*Ng) * Nq]),
					          0.5*(prim1[jiq +   (COLS+2*Ng) * Nq] - prim1[jiq - (COLS+2*Ng) * Nq]),
					       	theta*(prim1[jiq +   (COLS+2*Ng) * Nq] - prim1[jiq])                  );

				slopeR = minmod(theta*(prim1[jiq +   (COLS+2*Ng) * Nq] - prim1[jiq]                   ),
					       	0.5  *(prim1[jiq + 2*(COLS+2*Ng) * Nq] - prim1[jiq]                   ),
					       	theta*(prim1[jiq + 2*(COLS+2*Ng) * Nq] - prim1[jiq + (COLS+2*Ng) * Nq]));

				pL[q] = prim1[jiq]                    + 0.5*slopeL;
				pR[q] = prim1[jiq + (COLS+2*Ng) * Nq] - 0.5*slopeR;
			}

			sourceL[2] = pL[RHO]*gravity;
			sourceL[3] = pL[RHO]*pL[VYY]*gravity;
			sourceR[2] = pR[RHO]*gravity;
			sourceR[3] = pR[RHO]*pR[VYY]*gravity;

			primToCons(pL, uL);
			primToCons(pR, uR);
			
			for(q=0; q<Nq; ++q){
				uL[q] -= (Giph[Nq*(ROWS+2*Ng-1)*j + Nq*i + q]*dt/dy + dt*sourceL[q]);
				uR[q] +=  Giph[Nq*(ROWS+2*Ng-1)*j + Nq*i + q]*dt/dy + dt*sourceR[q] ;
			}

			consToPrim(uL, pL);
			consToPrim(uR, pR);
			
			for(q=0; q<Nq; ++q){
				int jiq = Nq*(COLS+2*Ng)*i + Nq*j + q;
				prim2[jiq]                  = pL[q];
				prim2[jiq + Nq*(COLS+2*Ng)] = pR[q];
			}
		}
	}

	////////////////////////////////////
	//Boundary conditions for X
	//Periodic
	for(i=0; i<(ROWS+2*Ng); ++i){
		for(q=0; q<Nq; ++q){
			//Left Side
			prim[Nq*(COLS+2*Ng)*i + Nq*0           + q] = prim[Nq*(COLS+2*Ng)*i + Nq*(COLS+0) + q];
			prim[Nq*(COLS+2*Ng)*i + Nq*1           + q] = prim[Nq*(COLS+2*Ng)*i + Nq*(COLS+1) + q];
			//Right Side
			prim[Nq*(COLS+2*Ng)*i + Nq*(COLS+Ng+1) + q] = prim[Nq*(COLS+2*Ng)*i + Nq*(Ng+1)   + q];
			prim[Nq*(COLS+2*Ng)*i + Nq*(COLS+Ng+0) + q] = prim[Nq*(COLS+2*Ng)*i + Nq*(Ng+0)   + q];

		}
	}

	//Boundary conditions for Y
	//Reflective
	for(j=0; j<(COLS+2*Ng); ++j){
		for(q=0; q<Nq; ++q){

			if(q == 3){
				//Bottom
				prim[Nq*(COLS+2*Ng)*0 + Nq*j + q] = -prim[Nq*(COLS+2*Ng)*Ng + Nq*j + q];
				prim[Nq*(COLS+2*Ng)*1 + Nq*j + q] = -prim[Nq*(COLS+2*Ng)*Ng + Nq*j + q];
				//Top
				prim[Nq*(COLS+2*Ng)*(ROWS+2*Ng - 2) + Nq*j + VYY] = -prim[Nq*(COLS+2*Ng)*(ROWS+1) + Nq*j + VYY];
				prim[Nq*(COLS+2*Ng)*(ROWS+2*Ng - 1) + Nq*j + VYY] = -prim[Nq*(COLS+2*Ng)*(ROWS+1) + Nq*j + VYY];
			}else{
				//Bottom
				prim[Nq*(COLS+2*Ng)*0 + Nq*j + q] = prim[Nq*(COLS+2*Ng)*Ng + Nq*j + q];
				prim[Nq*(COLS+2*Ng)*1 + Nq*j + q] = prim[Nq*(COLS+2*Ng)*Ng + Nq*j + q];
				//Top
				prim[Nq*(COLS+2*Ng)*(ROWS+2*Ng - 2) + Nq*j + VYY] = prim[Nq*(COLS+2*Ng)*(ROWS+1) + Nq*j + VYY];
				prim[Nq*(COLS+2*Ng)*(ROWS+2*Ng - 1) + Nq*j + VYY] = prim[Nq*(COLS+2*Ng)*(ROWS+1) + Nq*j + VYY];

			}
		}
	}

	///////////////////////////////////
	return(dt);
}
*/
int main(void){
	FILE *rtI;

	rtI = fopen("/home/danny/CompPhys/FinalProj/2ndOrder/2ndOrderData/rt0.txt", "w");

	int i, j, count;
	count = 0;
	double prim[Nq*(ROWS+2*Ng)*(COLS+2*Ng)];
	set_initial(prim);

	
	double dx, dy, dt;
        dx = L/(double)COLS;
	dy = H/(double)ROWS;
	double t = 0.0;
	dt = 1e-5;
	

/*	
	for(i=0; i<(ROWS+2*Ng); ++i){
		for(j=0; j<(COLS+2*Ng); ++j){
			printf("Zone(%d %d): P = %.2e, rho = %.2e, vx = %.2e, vy = %.2e\n",
				       	j, i, prim[Nq*(COLS+2*Ng)*i + Nq*j + 0], prim[Nq*(COLS+2*Ng)*i + Nq*j + 1], prim[Nq*(COLS+2*Ng)*i + Nq*j + 2], prim[Nq*(COLS+2*Ng)*i + Nq*j + 3]);
		}
	}

*/
	while(t<T){

		if (count == 1000) break;
		dt = step1(prim, dx, dy, dt);
		t+=dt;
		printf("%d dt = %e\n", count++ ,dt);
		if (dt<1e-20) break;
	}


	for(i=Ng; i<(ROWS+Ng); ++i){
		for(j=Ng; j<(COLS+Ng); ++j){
			fprintf(rtI, "%e %e %e %e %e %e\n", prim[Nq*(COLS+2*Ng)*i+Nq*j + 0], prim[Nq*(COLS+2*Ng)*i+Nq*j + 1], prim[Nq*(COLS+2*Ng)*i + Nq*j + 2], prim[Nq*(COLS+2*Ng)*i + Nq*j + 3], getx(j), gety(i));
		}
	}

	fclose(rtI);

	return(0);
}
