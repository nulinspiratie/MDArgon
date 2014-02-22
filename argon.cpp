#include <cstdlib>
#include <iostream>
#include <stdio.h>
#include <cmath>
#include <fstream>
#include <stdlib.h>
#include <ctime>
#include <cstring>
//#include <plplot/plstream.h>
//#include <armadillo>


#define rho 0.1			//number density
#define n 4 			//unit cells per direction
#define T 0.1			//Temperature (actually kT/m)
#define rc2 9			//Cutoff length squared

#define dt 0.00001		//timestep
#define iterations 1000000	//number of iterations
#define loopsperprint 10000	//loops before another print is made
#define listinterval 100	//Update verlet list every listinterval loops

#define notices 10

#define N 4*n*n*n		//number of particles
#define L pow(N/rho,1/3.0)	//box dimension
#define pi 3.141592653589793
#define ecut 4*(1/pow(rc2,6) - 1/pow(rc2,3))
#define rv2 rc2			//Verlet list cutoff distance

double pos[N][3]={};		//positions
double vel[N][3]={};		//velocities
double force[N][3]={};		//forces
double epot=0;			//potential energy
double ekin=0;			//kinetic energy
double t=0;
int nlist[N]={};		//Verlet particle number list
double list[N][N]={};

void updatelist();
void initialize();
void calcforce();
double dist(int i, int j,int coord=3);	//if coord=3, it returns the norm squared, else it returns the distance in coord dimension
void print();
double normalrand();
void displace();
double mod(double num, double div);

using namespace std;
//using namespace arma;
int main()
{
	cout << "\n\n\n\n";
	cout << "Molecular Dynamics simulation of Argon\n\n";
	cout << "Number of particles: " << N << endl;
	cout << "Density: " << rho << endl;
	cout << "Box length: " << L << endl;
	cout << "Temperature: " << T << endl;
	cout << "dt: " << dt << "s\t iterations: " << iterations << endl;
	cout << "simulation time: " << dt * iterations << " seconds\n";
	
	srand48(2);//(long)time(NULL));
	initialize();
	cout << "Starting Simulation\n";
	clock_t begin = clock();
	updatelist();
	calcforce();
	print();
	for (int perc=0; perc < notices; perc++)
	{
		for (int loop = 0; loop < iterations/notices; loop++)
		{
			if (loop%listinterval==0) updatelist();
			if (loop%loopsperprint==0) print();
			//calcforce();
			displace();
			t += dt;
		}
		cout << "Percent done: " << (perc+1) * 100 / notices << "%";
		cout << "\ttime simulated: " << double(clock() - begin) / CLOCKS_PER_SEC << " s";
		double timeleft = (double(notices - (perc+1)) / (perc+1) * double(clock() - begin)) / CLOCKS_PER_SEC;
		cout << "\t time left: " << floor(timeleft/60) << "m " << round(mod(timeleft,60)) << "s\n";
		}
	return 0;
}

void updatelist()
{
	memset(nlist,0,sizeof(nlist));
	for (int i=0; i<N-1; i++)
		for (int j=i+1; j<N; j++)
			if( dist(i,j,3) < rc2)
			{
				//nlist[j]++;
				list[i][nlist[i]] = j;
				nlist[i]++;
				//list[j][nlist[j]] = i;
			}
}

void displace()
{
	ekin=0;
	for (int p=0;p<N;p++)
	{
		for (int i=0;i<3;i++)
		{
			vel[p][i] += force[p][i]*dt/2.;
			ekin += 0.5 * vel[p][i] * vel[p][i];
			pos[p][i] = mod(pos[p][i] + vel[p][i] * dt, L);
		}
	}
	calcforce();
	//cout << "epot=" << epot << endl;
	//cin.ignore();
	for (int p=0;p<N;p++)
		for (int i=0;i<3;i++)
			vel[p][i] += force[p][i] * dt / 2.;
}




void calcforce()
{
	memset(force,0,sizeof(force));
	epot = 0;
	for (int i=0;i<N;i++)
	{
		for (int pn=0; pn<nlist[i]; pn++)
		{
			int j=list[i][pn];//This is the correct particle number
			//cout << "i,j=" << i << "," << j << "\t";
			//cin.ignore();
			double dr[3];
			dr[0]=dist(i,j,0);
			dr[1]=dist(i,j,1);
			dr[2]=dist(i,j,2);
			double r2i = 1./(dr[0]*dr[0] + dr[1]*dr[1] + dr[2]*dr[2]);
			//cout << r2i;
			//cin.ignore();
			double r6i = r2i * r2i * r2i;
			double f=48*r2i*r6i*(r6i-0.5);       //Lennard-Jones Potential
			for (int k=0;k<3;k++)
			{
				force[i][k] += f*dr[k];
				force[j][k] -= f*dr[k];
			}
			epot += 4*r6i *(r6i - 1) - ecut;

		}
	}
}

double dist(int i, int j, int coord)
{
	double dr=0;
	if (coord!=3) //Returns distance (may be negative)
	{
		dr = pos[i][coord] - pos[j][coord];
		dr -= round(dr / L) * L;
		return dr;
	}
	else //Returns the distance squared (always positive)
	{
		double dx;
		for (int k = 0; k < 3 ; k++)
		{
			dx = pos[i][k] - pos[j][k];
			dx -= round(dx / L) * L;
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
	if (filecount==1)
		file.open("output.dat");
	else	
		file.open("output.dat",fstream::app);
	file << ekin + epot << " " << ekin << " " << epot << endl;
	file.close();

	/*Temporary check for maximum velocity
	  double maxvel=0;
	  for (int i=0;i<N;i++)
	  for (int coord=0;coord<3;coord++)
	  if (abs(vel[i][coord])>maxvel) maxvel = vel[i][coord];
	  cout << "Maximum velocity = " << maxvel << endl;
	  */
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



