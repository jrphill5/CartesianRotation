#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <signal.h>
#include <unistd.h>

using namespace std;

void top();

void systemA(   double a,                     double dx, double dy, double dz, const char* color );
void systemAB(  double a, double b,           double dx, double dy, double dz, const char* color );
void systemABC( double a, double b, double c, double dx, double dy, double dz, const char* color );

bool initplot();
void plot( const char* command);

void printv( double* v );
void printm( double* m );
void drawv( int id, const char* color, double* o, double* v );

double dotv( double* v1, double* v2 );
double* unitv( double* v );
double magv( double* v );
void zerov( double* v );

double* transform( double* m, double* v );
double* transpose( double* m );

const double pi = 3.141592653589793;

void sighandler(int sig);
FILE* gnuplot = NULL;

double r  = 1.0;
double th = pi/6.0;
double ph = pi/4.0;

double a, a0 = 0;
double b, b0 = pi/8.0;
double c, c0 = 0;

double* vector;
double* origin;

double* trunk;
double* tier1;
double* tier2;

double xhat[3] = { 1.0, 0.0, 0.0 };
double yhat[3] = { 0.0, 1.0, 0.0 };
double zhat[3] = { 0.0, 0.0, 1.0 };

int main()
{

	if ( !initplot() )
	{

		printf("Plot initialization failed!\n");
		return 1;

	}

	vector = (double*) calloc(3, sizeof(double));
	origin = (double*) calloc(3, sizeof(double));

	trunk = (double*) calloc(3, sizeof(double));
	tier1 = (double*) calloc(3, sizeof(double));
	tier2 = (double*) calloc(3, sizeof(double));

	zerov( origin );

	trunk[0] = origin[0];
	trunk[1] = origin[1];
	trunk[2] = origin[2] + r;

	drawv( 0, "#000000", origin, trunk );

	double rotAB[]  = {  cos(ph)*cos(th), -sin(ph),  cos(ph)*sin(th),
	                     sin(ph)*cos(th),  cos(ph),  sin(ph)*sin(th),
	                    -sin(th),          0.0,      cos(th)  };

	int max = 3;
	for ( int i = 0; i < max; i++ )
	{

		double r1 = r;
		double a1 = 2.0*pi/max * i;
		double b1 = th;

		tier1[0] = 0.0;
		tier1[1] = 0.0;
		tier1[2] = r1;

		double dx = trunk[0];
		double dy = trunk[1];
		double dz = trunk[2];

		double rotAB1[]  = {  cos(a1)*cos(b1), -sin(a1),  cos(a1)*sin(b1),
		                      sin(a1)*cos(b1),  cos(a1),  sin(a1)*sin(b1),
		                     -sin(b1),          0.0,      cos(b1)  };

		systemAB( a1, b1, dx, dy, dz, "#00FF00" );

		vector = transform( rotAB1, tier1 );
		vector[0] += dx; vector[1] += dy; vector[2] += dz;
		drawv( 100*(i + 1), "#000000", trunk, vector );
		plot( "replot" );

		tier1[0] = vector[0];
		tier1[1] = vector[1];
		tier1[2] = vector[2];

		for ( int j = 0; j < max; j++ )
		{

			double r2 = r;
			double a2 = 2.0*pi/max * j;
			double b2 = th;

			tier2[0] = r2*sin(b2)*cos(a2);
			tier2[1] = r2*sin(b2)*sin(a2);
			tier2[2] = r2*cos(b2);

			double dx = tier1[0];
			double dy = tier1[1];
			double dz = tier1[2];

			double rotAB2[]  = {  cos(a2)*cos(b2), -sin(a2),  cos(a2)*sin(b2),
			                      sin(a2)*cos(b2),  cos(a2),  sin(a2)*sin(b2),
			                     -sin(b2),          0.0,      cos(b2)          };

			systemAB( a1, b1, dx, dy, dz, "#FF0000" );

			vector = transform( rotAB1, tier2 );
			vector[0] += dx; vector[1] += dy; vector[2] += dz;
			drawv( 100*(i + 1)+10*(j + 1), "#0000FF", tier1, vector );
			plot( "replot" );

			tier2[0] = vector[0];
			tier2[1] = vector[1];
			tier2[2] = vector[2];

			for ( int k = 0; k < max; k++ )
			{

			}

		}			

	}
	
	while ( true )
		signal(2, &sighandler);

	double o[3] = {0.0,0.0,0.0};

	drawv( 1, "#000000", o, origin );
	plot( "replot" );
	top();

	while ( true )
		signal(2, &sighandler);

	return 0;

}

bool initplot()
{

	if ( gnuplot == NULL )
	{

		gnuplot = popen("gnuplot -geometry 600x600 - > /dev/null 2>&1","w");
		if ( gnuplot == NULL ) return false;
		//plot("set term postscript eps color");
		plot("set xlabel \"x\" offset first 3,0,0");
		plot("set xrange[-3:3]");
		plot("set ylabel \"y\" offset first 0,3,0");
		plot("set yrange[-3:3]");
		plot("set zlabel \"z\" offset first 0,0,1.5");
		plot("set zrange[0:3]");
		plot("set view equal_axes xyz");
		plot("set view 60,135");
		plot("unset key");
		plot("unset border");
		plot("set ticslevel 0");
		plot("set zeroaxis lw 2 lt 1 lc rgb \"#000000\"");
		plot("set xtics 1 axis nomirror");
		plot("set ytics 1 axis nomirror");
		plot("set ztics 1 axis nomirror");
		plot("f(x,y) = -1");
		plot("splot f(x,y)");
		plot("replot");
		return true;

	}
	else return false;

}

void top()
{

	a0 = 0.0;
	b0 = pi/6.0;
	c0 = 0.0;

	int speed = 0;

	while ( speed < 1 || speed > 100 )
	{

		printf("speed: ");
		scanf("%d",&speed);

	}

	int kmax = 100;
	int jmax = 1;
	int imax = 100;

	while ( true )
	{

		for ( int k = 0; k < kmax; k++ )
		{

			double a = a0 + 2.0*pi/kmax * k;
			systemA( a, 0,0,0, "#FF0000" );

			for ( int j = 0; j < jmax; j++ )
			{

				double b = b0 + 2.0*pi/jmax * j;
				systemAB( a, b, 0,0,0, "#00FF00" );

				for ( int i = 0; i < imax; i++ )
				{

					double c = c0 + 2.0*pi/imax * i;
					systemABC( a, b, c, 0,0,0, "#0000FF" );
					plot("replot");

					usleep( 750*(100 - speed) );
					signal(2, &sighandler);

				}

			}

		}

	}

}

void systemA( double a, double dx, double dy, double dz, const char* color )
{

	double origin[3] = { dx, dy, dz };

	double rotA[]   = {  cos(a), -sin(a),  0.0,
	                     sin(a),  cos(a),  0.0,
	                     0.0,     0.0,     1.0     };	

	vector = transform( rotA, xhat );
	vector[0] += dx; vector[1] += dy; vector[2] += dz;
	drawv( 11, color, origin, vector );

	vector = transform( rotA, yhat );
	vector[0] += dx; vector[1] += dy; vector[2] += dz;
	drawv( 12, color, origin, vector );

	vector = transform( rotA, zhat );
	vector[0] += dx; vector[1] += dy; vector[2] += dz;
	drawv( 13, color, origin, vector );

}

void systemAB( double a, double b, double dx, double dy, double dz, const char* color )
{

	double origin[3] = { dx, dy, dz };

	double rotAB[]  = {  cos(a)*cos(b), -sin(a),  cos(a)*sin(b),
	                     sin(a)*cos(b),  cos(a),  sin(a)*sin(b),
	                    -sin(b),         0.0,     cos(b)  };

	vector = transform( rotAB, xhat );
	vector[0] += dx; vector[1] += dy; vector[2] += dz;
//	drawv( 21, color, origin, vector );
	drawv( 0, color, origin, vector );

	vector = transform( rotAB, yhat );
	vector[0] += dx; vector[1] += dy; vector[2] += dz;
//	drawv( 22, color, origin, vector );
	drawv( 0, color, origin, vector );

	vector = transform( rotAB, zhat );
	vector[0] += dx; vector[1] += dy; vector[2] += dz;
//	drawv( 23, color, origin, vector );
	drawv( 0, color, origin, vector );

}

void systemABC( double a, double b, double c, double dx, double dy, double dz, const char* color )
{

	double origin[3] = { dx, dy, dz };

	double rotABC[] = {  cos(a)*cos(b)*cos(c) - sin(a)*sin(c), -sin(a)*cos(c) - cos(a)*cos(b)*sin(c),  cos(a)*sin(b),
	                     sin(a)*cos(b)*cos(c) + cos(a)*sin(c),  cos(a)*cos(c) - sin(a)*cos(b)*sin(c),  sin(a)*sin(b),
	                    -sin(b)*cos(c),                         sin(b)*sin(c)                       ,  cos(b)         };

	vector = transform( rotABC, xhat );
	vector[0] += dx; vector[1] += dy; vector[2] += dz;
	drawv( 31, color, origin, vector );

	vector = transform( rotABC, yhat );
	vector[0] += dx; vector[1] += dy; vector[2] += dz;
	drawv( 32, color, origin, vector );

	vector = transform( rotABC, zhat );
	vector[0] += dx; vector[1] += dy; vector[2] += dz;
	drawv( 33, color, origin, vector );

}

void sighandler(int sig)
{

	plot("exit");
	fclose(gnuplot);
	exit(1);

}

void plot( const char* command )
{

	fprintf(gnuplot, "%s\n", command);
	fflush(gnuplot);

}

void printv( double* v )
{

	for ( int i = 0; i < 3; i++ )
		printf("[ % 1f ]", v[i]);
	printf(" %f",magv(v));
	printf("\n");

}

void printm( double* m )
{

	for ( int j = 0; j < 3; j++ )
	{
		printf("[ ");
		for ( int i = 0; i < 3; i++ )
			printf("% 1f ", m[3*j+i]);
		printf("]\n");
	}
	printf("\n");

}

double dot( double* v1, double* v2 )
{

	return v1[0]*v2[0] + v1[1]*v2[1] + v1[2]*v2[2];

}

double* unitv( double* v )
{

	static double vector[3] = {0.0};
	double magnitude = magv( v );

	vector[0] /= magnitude;
	vector[1] /= magnitude;
	vector[2] /= magnitude;

	return vector;

}

double magv( double* v )
{

	return pow( dot(v, v), 0.5 );

}

void zerov( double* v )
{

	v[0] = 0.0;
	v[1] = 0.0;
	v[2] = 0.0;

}

double* transform( double* m, double* v )
{

	static double vector[3] = {0};

/*	for ( int j = 0; j < 3; j++ )
		for ( int i = 0; i < 3; i++ )
			vector[j] += m[3*j+i]*v[i];*/

	vector[0] = m[0]*v[0] + m[1]*v[1] + m[2]*v[2];
	vector[1] = m[3]*v[0] + m[4]*v[1] + m[5]*v[2];
	vector[2] = m[6]*v[0] + m[7]*v[1] + m[8]*v[2];

	return vector;

}

double* transpose( double* m )
{

	static double matrix[9] = {0};

	for ( int j = 0; j < 3; j++ )
        for ( int i = 0; i < 3; i++ )
            matrix[j] += m[3*i+j];

	return matrix;

}

void drawv( int id, const char* color, double* o, double* v )
{

	char command[128];
	if ( id <= 0 ) sprintf(command,"set arrow    from %f,%f,%f to %f,%f,%f lc rgb \"%s\"",     o[0],o[1],o[2], v[0],v[1],v[2], color);
	else           sprintf(command,"set arrow %d from %f,%f,%f to %f,%f,%f lc rgb \"%s\"", id, o[0],o[1],o[2], v[0],v[1],v[2], color);
	plot(command);

}
