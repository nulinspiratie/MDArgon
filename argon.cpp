#include <cstdlib>
#include <iostream>
#include <stdio.h>
#include <math.h>
#include <fstream>
#include <stdlib.h>
#include <ctime>
#include <cstring>
//#include <plplot/plstream.h>
//#include <armadillo>


#define rho 0.1			//number density
#define n 2 			//unit cells per direction
#define T 1			//Temperature (actually kT/m)
#define rc2 9			//Cutoff length squared

#define dt 0.00001		//timestep
#define iterations 1000000		//number of iterations
#define loopsperprint 10000	//loops before another print is made

#define notices 10

#define N 4*n*n*n		//number of particles
#define L pow(N/rho,1/3.0)	//box dimension
#define pi 3.141592653589793
#define ecut 4*(1/pow(rc2,6) - 1/pow(rc2,3))

double pos[N][3]={};		//positions
double vel[N][3]={};		//velocities
double force[N][3]={};		//forces
double energy=0;		//energy
double t=0;


void initialize();
void calcforce();
double dist(int i, int j,int coord=3);	//if coord=3, it returns the norm, else it returns the distance in coord dimension
void print();
double normalrand();
void displace();
double mod(double num, double div);

using namespace std;
//using namespace arma;
int main()
{
	srand48(2);//(long)time(NULL));
	if (rho > sqrt(2))//Must find correct value
	{
		cout << "Density cannot be larger than sqrt(2)" << endl;
		return 1;
	}
	/*if (n%2==1)
	{
		cout << "Must be an even number of particles" << endl;
		return 1;
	}*/
	initialize();
	print();
	clock_t begin = clock();
	for (int perc=0; perc < notices; perc++)
	{
		for (int loop = 0; loop < iterations/notices; loop++)
		{
			if (loop%loopsperprint==0) print();
			//calcforce();
			displace();
			t += dt;
		}
		cout << "Percent done: " << (perc+1) * 100 / notices << "%";
		cout << "\ttime simulated: " << double(clock() - begin) / CLOCKS_PER_SEC;
		cout << "\t time left: " << ((notices - perc) / (perc+1) * double(clock() - begin)) / CLOCKS_PER_SEC << endl;
	}
	return 0;
}


void displace()
{
	for (int p=0;p<N;p++)
	{
		for (int i=0;i<3;i++)
		{
			vel[p][i]+=force[p][i]*dt/2.;
			pos[p][i] = mod(pos[p][i] + vel[p][i] * dt, L);
		}
	}
	calcforce();
	for (int p=0;p<N;p++)
		for (int i=0;i<3;i++)
			vel[p][i] += force[p][i] * dt / 2.;
}




void calcforce()
{
	energy=0;
	memset(force,0,sizeof(force));
	for (int i=0; i<N; i++)
	{
		for (int j=i+1; j<N; j++)
		{
			double dr[3];
			dr[0]=dist(i,j,0);
			dr[1]=dist(i,j,1);
			dr[2]=dist(i,j,2);
			double dr2 = dr[0]*dr[0] + dr[1]*dr[1] + dr[2]*dr[2];
			if (dr2 < rc2)
			{
				double r2i = 1/dr2;
				double r6i = pow(r2i,3);
				double f=48*r2i*r6i*(r6i-0.5); 	//Lennard-Jones Potential
				for (int k=0;k<3;k++)
				{
					force[i][k] += f*dr[k];
					force[j][k] -= f*dr[k];
				}
				energy += 4*r6i *(r6i - 1) - ecut;
				//printf("{%f,%f,%f}, {%f,%f,%f}, dr={%f,%f,%f}\nf=(%f,%f,%f)\n",pos[i][0],pos[i][1],pos[i][2],pos[j][0],pos[j][1],pos[j][2],dr[0],dr[1],dr[2],force[i][0],force[i][1],force[i][2]);
			}
		}
	}
}

double dist(int i, int j, int coord)
{
	double dr=0;
	if (coord!=3)
	{
		dr = pos[i][coord] - pos[j][coord];
		dr = dr - round(dr / L) * L;
		return dr;
	}
	else
	{
		double dx;
		for (int k = 0; k < 3 ; k++)
		{
			dx = pos[i][k] - pos[j][k];
			dx = dx - round(dx / L) * L;
			dr += dx * dx;
		}
		return dr;
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
	file.open("output.dat",fstream::app);
	file << energy << endl;
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
	//FCC Lattice is created
	calcforce();
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


