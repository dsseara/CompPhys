/* Daniel Seara   Andrew MacFadyen    Paul Duffel
Computational Physics Assignment 1: Particle in Gravity Field.
This code simulates a point particle in a 2 dimensional gravity well.
The equation of motion is simply:
	x'' = -G*M/r^2
where GM characterizes the strength of the field, and r is the distance
between the particle and source. These equations can be reduced to 
two couple first order ODEs:
	x' = v
	v' = -G*M/r^2
The code integrates each variable separately and writes the data out to a
file which can be used for visualization. All three methods are run using the
same initial conditions.

As well as integrate the trajectory, the "integrate" function also calculates
the error in terms the energy and angular momentum loss as a function of time.
*/

#include <stdio.h>
#include <math.h>

#define GM 1.0

/*Method = 0 -> Euler
 *Method = 1 -> RK2
 *Method = 2 -> RK4
 *T is the final time, N is the number of steps, and dt is the size of the time step
 */

int method;
double T = 50;
int N = 10;
double dt;

double integrate(void);
void   euler(double, double, double, double, double, double*, double*, double*, double*);
void   rk2(double, double, double, double, double, double*, double*, double*, double*);
void   rk4(double, double, double, double, double, double*, double*, double*, double*);
double dx_dt(double, double, double, double ,double);
double dy_dt(double, double, double, double, double);
double dvx_dt(double, double, double);
double dvy_dt(double, double, double);



int main(){
	FILE *l2Norm;
	int i;
	int j, jmax;
	double l2;
	jmax = 55;
	int n;
	/* Use this section to get the L2 norm for a circular orbit.
	 * For each method, run through all values of n, the number of time steps, then do it 
	 * again for the other methods.*/
	for (i=0; i<3; ++i){
		
		n = N;	
		method = i;
		if (method == 0){
			l2Norm = fopen("l2NormE.txt", "w");
		} else if (method==1){
			l2Norm = fopen("l2NormRK2.txt", "w");
		} else if (method==2){
			l2Norm = fopen("l2NormRK4.txt", "w");
		}

		for (j = 0 ; j<jmax; ++j){

			dt = T/n;
			l2 = integrate();
		/*	printf("l2 for method %d with n = %d: %25.22E\n", i, n, l2); */

			fprintf(l2Norm, "%25.16f %25.16f\n", dt, l2);
			n = n*2;
		}
	}  
 	/*Just for getting trajectories. 
	 *Also uncomment the fopen command in integrate */
/*	for (i=0; i<3; ++i){
		dt = T/N;
		method = i;
		integrate();
	}  
*/
	return 0;
}

double integrate(void){
	FILE *out;
	
	/*Initial Values*/
	double  t = 0.0;
	double vx = 0.0;
	double vy = -1.0;
	double  x = 1.0;
	double  y = 0.0;

	double x_new, y_new, vx_new, vy_new;
	double E0, E, dE;
	double L0, L, dL;
	double l2sq = 0.0;

	int i;

	if (method == 0){
		out = fopen("gravityEuler.txt", "w");
	} else if (method == 1){
		out = fopen("gravityRK2.txt", "w");
	} else if (method == 2){
		out = fopen("gravityRK4.txt", "w");
	}
	
	/* Initial values of energy (E =0.5 * m (vx^2 + vy^2) + GM/r) and angular momentum
	 * (L = r x p = x*m*vy - y*m*vx). Note that the mass of the particle m=1. */
	E0 = -0.5; 
	L0 = 1.0;

	/* Now begins the actual integrating. Depending on the method called, this function
	 * uses any of the 3 methods mentioned above, as well as calculating the energy and
	 * angular momentum at each time step to be used in calculating the error over time.*/
	for (i = 0; i < N; i++){
	
		E = 0.5*(pow(vx,2) + pow(vy,2)) - GM/( pow( pow(x,2) + pow(y,2), 0.5));
		dE = fabs((E - E0) / E0);

		L = x * vy - y * vx;
		dL = fabs((L - L0) / L0);

	/*	printf("%f %f %f %f %f %f %f %d\n",t, x, y, vx, vy, dE, dL, i);  */
		fprintf(out, "%f %25.16f %25.16f %25.16f %25.16f %25.16f %25.16f %d\n",t, x, y, vx, vy, dE, dL, i);

		if (method == 0){
			euler(t, x, y, vx, vy, &x_new, &y_new, &vx_new, &vy_new);
		} else if (method == 1){
			rk2(t, x, y, vx, vy, &x_new, &y_new, &vx_new, &vy_new);
		} else if (method == 2){
			rk4(t, x, y, vx, vy, &x_new, &y_new, &vx_new, &vy_new);
		}
		
		if (t > T){
			break;
		}
		
		/*Find the L2 norm as a function of step size, h.
		 * L2 = sqrt(int_(0)^(T) || r_num(t) - r_real(t)||^2 dt)
		 * r_real(t) = 1 for circular orbit with initial conditions
		 * r0 = (1,0), v0 = (0,1).
		 */  

		double r_num =  pow( (pow(x,2) + pow(y,2)) , 0.5);
		l2sq += fabs( pow((r_num - 1), 2)) *dt; 
	

		t  = t+dt;
		x  = x_new;
		y  = y_new;
		vx = vx_new;
		vy = vy_new;
	}
	fclose(out);
	return pow(l2sq, 0.5);
}


double dvx_dt(double t, double x, double y){
	double dsq = (x*x + y*y);
	double gm = GM;
	double force = -1*gm / dsq;
	double theta = atan(y/x);
	double ax = force*cos(theta);
	if (x<0) {
		ax = fabs(ax);
	} else {
		ax = fabs(ax)*-1;
	}
	return ax;
}

double dvy_dt(double t, double x, double y){
	double dsq = (x*x + y*y);
	double gm = GM;
	double force = -1*gm / dsq;
	double theta = atan(y/x);
	double ay = force*sin(theta);
	if (y<0) {
		ay = fabs(ay);
	} else {
		ay = fabs(ay)*-1;
	}	
	return ay;
}

double dx_dt(double t, double x, double y, double vx, double vy){
	return vx;
}

double dy_dt(double t, double x, double y, double vx, double vy){
	return vy;
}

void euler(double t, double x, double y, double vx, double vy, double* x_new, double* y_new, double* vx_new, double* vy_new){

	/*Forward Euler Method */	
	double h = dt;
		
	*x_new  =  x + h * dx_dt(t, x, y, vx, vy);	
	*y_new  =  y + h * dy_dt(t, x, y, vx, vy);
	*vx_new = vx + h * dvx_dt(t, x, y);
	*vy_new = vy + h * dvy_dt(t, x, y);
}

void rk2(double t, double x, double y, double vx, double vy, double* x_new, double* y_new, double* vx_new, double* vy_new){

	/* Mid-point Method */
	double h = dt;

	/*k's are changes in velocity, l's for position*/
	double k1x, k1y, l1x, l1y;

	k1x = h * dvx_dt(t, x, y);
	k1y = h * dvy_dt(t, x, y);
	l1x = h *  dx_dt(t, x, y, vx, vy);
	l1y = h *  dy_dt(t, x, y, vx, vy);

	*x_new  =  x + h * dx_dt (t + 0.5*h, x + 0.5*l1x, y + 0.5*l1y, vx + 0.5 * k1x, vy + 0.5 * k1y);
	*y_new  =  y + h * dy_dt (t + 0.5*h, x + 0.5*l1x, y + 0.5*l1y, vx + 0.5 * k1x, vy + 0.5 * k1y);
	*vx_new = vx + h * dvx_dt(t + 0.5*h, x + 0.5*l1x, y + 0.5*l1y);
	*vy_new = vy + h * dvy_dt(t + 0.5*h, x + 0.5*l1x, y + 0.5*l1y);

}

void rk4(double t, double x, double y, double vx, double vy, double* x_new, double* y_new, double* vx_new, double* vy_new){
	/*RK4*/

	double h = dt;
	
	/*k's are changes in velocity, l's for position*/
	double l1x, l2x, l3x, l4x;
	double l1y, l2y, l3y, l4y;
	double k1x, k2x, k3x, k4x;
	double k1y, k2y, k3y, k4y;
		
	l1x = h *  dx_dt(t, x, y, vx, vy);
	l1y = h *  dy_dt(t, x, y, vx, vy);
	k1x = h * dvx_dt(t, x, y);
	k1y = h * dvy_dt(t, x, y);

	l2x = h *  dx_dt(t+0.5*dt, x+0.5*l1x, y+0.5*l1y, vx+0.5*k1x, vy+0.5*k1y);
	l2y = h *  dy_dt(t+0.5*dt, x+0.5*l1x, y+0.5*l1y, vx+0.5*k1x, vy+0.5*k1y);
	k2x = h * dvx_dt(t+0.5*dt, x+0.5*l1x, y+0.5*l1y);
	k2y = h * dvy_dt(t+0.5*dt, x+0.5*l1x, y+0.5*l1y);

	l3x = h *  dx_dt(t+0.5*dt, x+0.5*l2x, y+0.5*l2y, vx+0.5*k2x, vy+0.5*k2y);
	l3y = h *  dy_dt(t+0.5*dt, x+0.5*l2x, y+0.5*l2y, vx+0.5*k2x, vy+0.5*k2y);
	k3x = h * dvx_dt(t+0.5*dt, x+0.5*l2x, y+0.5*l2y);
	k3y = h * dvy_dt(t+0.5*dt, x+0.5*l2x, y+0.5*l2y);

	l4x = h *  dx_dt(t+0.5*dt, x+l3x, y+l3y, vx+k3x, vy+k3y);
	l4y = h *  dy_dt(t+0.5*dt, x+l3x, y+l3y, vx+k3x, vy+k3y);	
	k4x = h * dvx_dt(t+0.5*dt, x+l3x, y+l3y);
	k4y = h * dvy_dt(t+0.5*dt, x+l3x, y+l3y);
		
	*x_new  =  x + (l1x + 2.*l2x + 2.*l3x + l4x)/6.0;
	*y_new  =  y + (l1y + 2.*l2y + 2.*l3y + l4y)/6.0;
	*vx_new = vx + (k1x + 2.*k2x + 2.*k3x + k4x)/6.0;
     	*vy_new = vy + (k1y + 2.*k2y + 2.*k3y + k4y)/6.0;
}
