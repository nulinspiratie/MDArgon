#include <cstdlib>
#include <iostream>
#include <stdio.h>
#include <math.h>
#include <fstream>
#include <stdlib.h>
//#include <plplot/plstream.h>
//#include <armadillo>


#define n 2 			//unit cells per direction
#define T 5			//Temperature (actually kT/m)
#define dt 0.001			//timestep
#define iterations 1000		//number of iterations
#define loopsperprint 10	//loops before another print is made


#define N 4*n*n*n		//number of particles
#define rho 0.2			//number density
#define L pow(N/rho,1/3.0)	//box dimension
#define pi 3.141592653589793

double pos[N][3]={};		//positions
double vel[N][3]={};		//velocities
double t=0;

void initialize();
void print();
double normalrand();
void displace(int p);
double mod(double num, double div);

using namespace std;
//using namespace arma;
int main()
{
	if (rho > sqrt(2))//Must find correct value
	{
		cout << "Density cannot be larger than sqrt(2)" << endl;
		return 1;
	}
	if (n%2==1)
	{
		cout << "Must be an even number of particles" << endl;
		return 1;
	}
	initialize();
	print();
	for (int loop = 0; loop < iterations; loop++)
	{
		if (loop%loopsperprint==0) print();
		for (int i = 0; i < N; i++)
		{
			displace(i);
		}
		t += dt;
	}
	return 0;
}

void displace(int p)
{
	for (int i=0;i<3;i++)
	{
		pos[p][i] = mod( pos[p][i] +vel[p][i] * dt, L);
	}
}

void print()
{
	static int filecount=1;
	char filename[100];
	sprintf(filename,"images/fcc-%04d.dat",filecount);
	ofstream file(filename);
	file << N << endl;
	file << "0 " << L << endl;
	file << "0 " << L << endl;
	file << "0 " << L << endl;
	for (int i=0;i<N;i++)
		file << pos[i][0] << " " << pos[i][1] << " " << pos[i][2] << " 1" << endl;
	file.close();
	filecount++;

}

void initialize()
{
	double d=pow(4/rho,1/3.0);
	int i=0;
	for (int z=0;z<2*n;z++)
		for (int y=0;y<n;y++)
			for (int x=0;x<n;x++)
			{

				pos[i][0] = (x * d + (d/2) * (z % 2));
				pos[i][1] = y * d;
				pos[i][2] = z * d/2;
				for (int j=0;j<3;j++)
					vel[i][j]=normalrand();
				i++;
				pos[i][0] = (x + 1/2.) * d + (d/2) * (z % 2);
				pos[i][1] = (y + 1/2.) * d;
				pos[i][2] = z * d / 2;
				for (int j=0;j<3;j++)
					vel[i][j]=normalrand();
				i++;
			}
}

double normalrand()
{
	static double second;
	static bool extra=0;

	if (extra==1)
	{
		extra=0;
		return second;
	}
	else
	{
		double x = drand48();
		double y = drand48();
		double u = T * sqrt(-2*log(x))*cos(2*pi*y);
		second = T * sqrt(-2*log(x))*sin(2*pi*y);
		extra = 1;
		return u;
	}
}

double mod(double num,double div)
{
	return num - floor(num/div) * div;
}

/*void print()
  {
  plclear();
  plcol0(1);//Sets color (red)
  plbox3("bnstu", "x", 0, 0, "bnstu", "y", 0, 0, "bcnmstuv", "z", 0, 0);//Draws box with axes
  plcol0(2);//Sets color

//Now the tricky part: getting the rows of pos into plpoin3
//cout << pos;
double x[N],y[N],z[N];
for (int i=0;i<N;i++)
{
x[i]=pos[i][0];
y[i]=pos[i][1];
z[i]=pos[i][2];
}
plpoin3(N,x,y,z,4);
plflush();

Below should be used to initialise
plsdev("xcairo");
plinit();
pladv(0);
plvpor(0,1,0,1);
plwind(-1, 1, -2/3, 4 / 3);
plw3d(1, 1, 1, 0, box[0], 0, box[1], 0, box[2], 30, -45);//Viewing angle and sizes
print();
cin.ignore();
plspause(false);
plend();



}*/


