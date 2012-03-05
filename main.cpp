#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <signal.h>
#include <unistd.h>

using namespace std;

bool initplot();
void top();
void systemA( double a );
void systemAB( double a, double b );
void systemABC( double a, double b, double c );
void plot( const char* command);
void printv( double* v );
void printm( double* m );
double dot( double* v1, double* v2 );
double* unitv( double* v );
double length( double* v );
double* transform( double* m, double* v );
double* transpose( double* m );
void drawv( int id, const char* color, double* o, double* v );

const double pi = 3.141592653589793;

void sighandler(int sig);
FILE* gnuplot = NULL;

double r = 1.0;
double th = pi/4.0;
double ph = pi/2.0;

double a0 = 0;
double b0 = pi/6.0;
double c0 = 0;

double* vector;
double* origin;
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

	origin[0] = 1.0;
	origin[1] = 0.0;
	origin[2] = 0.0;

	double vect[3] = { origin[0] + r*sin(th)*cos(ph),
	                   origin[1] + r*sin(th)*sin(ph),
	                   origin[2] + r*cos(th)          };

	drawv( 1, "#000000", origin, vect );
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

	origin[0] = 0.0;
	origin[1] = 0.0;
	origin[2] = 0.0;

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
			systemA( a );

			for ( int j = 0; j < jmax; j++ )
			{

				double b = b0 + 2.0*pi/jmax * j;
				systemAB( a, b );

				for ( int i = 0; i < imax; i++ )
				{

					double c = c0 + 2.0*pi/imax * i;
					systemABC( a, b, c );
					plot("replot");

					usleep( 750*(100 - speed) );
					signal(2, &sighandler);

				}

			}

		}

	}

}

void systemA( double a )
{

	double rotA[]   = {  cos(a), -sin(a),  0.0,
	                     sin(a),  cos(a),  0.0,
	                     0.0,     0.0,     1.0     };	

	vector = transform( rotA, xhat );
	drawv( 11, "#FF0000", origin, vector );

	vector = transform( rotA, yhat );
	drawv( 12, "#FF0000", origin, vector );

	vector = transform( rotA, zhat );
	drawv( 13, "#FF0000", origin, vector );

}

void systemAB( double a, double b )
{

	double rotAB[]  = {  cos(a)*cos(b), -sin(a),  cos(a)*sin(b),
	                     sin(a)*cos(b),  cos(a),  sin(a)*sin(b),
	                    -sin(b),         0.0,     cos(b)  };


	vector = transform( rotAB, xhat );
	drawv( 21, "#00FF00", origin, vector );

	vector = transform( rotAB, yhat );
	drawv( 22, "#00FF00", origin, vector );

	vector = transform( rotAB, zhat );
	drawv( 23, "#00FF00", origin, vector );

}

void systemABC( double a, double b, double c )
{

	double rotABC[] = {  cos(a)*cos(b)*cos(c) - sin(a)*sin(c), -sin(a)*cos(c) - cos(a)*cos(b)*sin(c),  cos(a)*sin(b),
	                     sin(a)*cos(b)*cos(c) + cos(a)*sin(c),  cos(a)*cos(c) - sin(a)*cos(b)*sin(c),  sin(a)*sin(b),
	                    -sin(b)*cos(c),                         sin(b)*sin(c)                       ,  cos(b)         };

	vector = transform( rotABC, xhat );
	drawv( 31, "#0000FF", origin, vector );

	vector = transform( rotABC, yhat );
	drawv( 32, "#0000FF", origin, vector );

	vector = transform( rotABC, zhat );
	drawv( 33, "#0000FF", origin, vector );

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
	printf(" %f",length(v));
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
	double magnitude = length( v );

	vector[0] /= magnitude;
	vector[1] /= magnitude;
	vector[2] /= magnitude;

	return vector;

}

double length( double* v )
{

	return pow( dot(v, v), 0.5 );

}

double* transform( double* m, double* v )
{

	static double vector[3] = {0.0};

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
	sprintf(command,"set arrow %d from %f,%f,%f to %f,%f,%f lc rgb \"%s\"", id, o[0],o[1],o[2], v[0],v[1],v[2], color);
	plot(command);

}
