#include <cstdlib>
#include <iostream>
#include <stdio.h>
#include <cmath>
#include <fstream>
#include <stdlib.h>
#include <ctime>
#include <cstring>
#include <sys/stat.h>
//#include <plplot/plstream.h>
//#include <armadillo>


#define rho 0.1			//number density
#define n 4 			//unit cells per direction
#define T 10			//Temperature (actually kT/m)
#define rc2 9			//Cutoff length squared

#define dt 0.00000000001		//timestep
#define iterations 500000	//number of iterations

#define thermiter 100000	//number of thermalization iterations
#define rescaleiter 1000	//number of iterations between temperature rescaling

#define loopsperdatastore 1000	//loops before data is stored
#define loopsperthermdatastore 100	//loops before thermalization data is stored
#define loopsperimage 100000	//loops before another printed image is made
#define listinterval 100	//Update verlet list every listinterval loops

#define notices 10

#define N 4*n*n*n		//number of particles
#define L pow(N/rho,1/3.0)	//box dimension
#define pi 3.141592653589793
#define ecut 4*(1/pow(rc2,6) - 1/pow(rc2,3))
#define rv2 rc2+0.5		//Verlet list cutoff distance

double pos[N][3]={};		//positions
double posinit[N][3]={};	//Positions right after thermalization (For diffusion)
double vel[N][3]={};		//velocities
double force[N][3]={};		//forces
double epot=0;			//potential energy
double ekin=0;			//kinetic energy
double einit=0;			//initial energy
double t=0;
double tmeasurestart=0;		//Actual measurement starting time
double diffusion=0;		//Diffusion
clock_t begin;
int nlist[N]={};		//Verlet particle number list
double list[N][N]={};
char folder[100];
char filename[100];

double instantaneoustemp();
void updatelist();
void initialize();
void calcforce();
void calcekin();			//To calculate the initial kinetic energy
double dist(int i, int j,int coord=3);	//if coord=3, it returns the norm squared, else it returns the distance in coord dimension
void printdata(bool thermalize=0,double time=0);
void printimage();
void output(int stage);
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
	cout << "Thermalization time: " << dt * thermiter << " seconds\n";
	cout << "simulation time: " << dt * iterations << " seconds\n";

	//Create folders and files
	sprintf(folder,"data/%g-%g-%d-%g/",(double)T,(double)rho,n,dt);
	mkdir(folder, S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);
	//sprintf(folder,"data/%g-%g-%d-%g/images",(double)T,(double)rho,n);
	sprintf(filename,"%simages",folder);
	mkdir(filename, S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);
	//sprintf(folder,"data/%g-%g-%d-%g/data.dat",(double)T,(double)rho,n);
	sprintf(filename,"%sdata.dat",folder);
	ofstream file(filename);
	file.close();
	//sprintf(folder,"data/%g-%g-%d-%g/thermdata.dat",(double)T,(double)rho,n);
	sprintf(filename,"%sthermdata.dat",folder);
	file.open(filename);
	file.close();
	
	//Check if not too many images will be made
	if ( (iterations + thermiter) / loopsperimage > 1000)
	{
		cout << "Too many images will be made, exiting...\n";
		return 1;
	}

	srand48((long)time(NULL));
	initialize();
	begin = clock();	//Start timer

	updatelist();		//Create first verlet list
	calcforce();		//Calculate initial force and potential energy
	calcekin();		//Calculate initial kinetic energy
	einit = epot + ekin;
	printimage();
	output(0);//Also write output to file
	cout << "Starting thermalization\n";
	for (int perc=0; perc < notices; perc++)
	{
		for (int loop = 0; loop < thermiter/notices; loop++)
		{
			if (loop%listinterval==0) updatelist();
			if (loop%loopsperimage==0) printimage();
			if (loop%loopsperthermdatastore==0) printdata(1);
			
			if (loop%rescaleiter==0)//Rescale velocities for correct temperature
			{
				double lambda = sqrt( (N-1) * 3 * T / (2 * N * ekin));
				for (int i=0;i<N;i++)
					for (int j=0;j<3;j++)
						vel[i][j] *= lambda;

			}

			displace();
			t += dt;
		}

		cout << "Percent done: " << (perc+1) * 100 / notices << "%";
		cout << "\ttime simulated: " << double(clock() - begin) / CLOCKS_PER_SEC << " s";
		double timeleft = (double(notices - (perc+1)) / (perc+1) * double(clock() - begin)) / CLOCKS_PER_SEC;
		cout << "\t time left: " << floor(timeleft/60) << "m " << round(mod(timeleft,60)) << "s\n";
	}
	output(1);	//Write thermalization time to file
	cout << "\nStarting Simulation\n";
	tmeasurestart=t;			//Measurement starts at this time
	for (int i=0;i<N;i++)			//Initial positions for diffusion
		for (int j=0;j<3;j++)
			posinit[i][j] = pos[i][j];
	cout << "Temperature equals " << instantaneoustemp() << endl;

	for (int perc=0; perc < notices; perc++)
	{
		for (int loop = 0; loop < iterations/notices; loop++)
		{
			if (loop%listinterval==0) updatelist();
			if (loop%loopsperimage==0) printimage();
			if (loop%loopsperdatastore==0) printdata(0,t-tmeasurestart);
			displace();
			t += dt;
		}
		cout << "Percent done: " << (perc+1) * 100 / notices << "%";
		cout << "\ttime simulated: " << double(clock() - begin) / CLOCKS_PER_SEC << " s";
		double total = double( (thermiter + iterations));
		double current = (double)(thermiter +  double(perc+1)/notices * iterations);
		double timeleft = ((total/current - 1)  * double(clock() - begin)) / CLOCKS_PER_SEC;
		cout << "\t time left: " << floor(timeleft/60) << "m " << round(mod(timeleft,60)) << "s\n";
	}
	diffusion=0;
	for (int i=0;i<N;i++)
		for (int j=0;j<3;j++)
			diffusion+=(pos[i][j] - posinit[i][j]) * (pos[i][j] - posinit[i][j]);
	diffusion/= 6. * N * (t-tmeasurestart);


	cout << "\nSimulation Finished\n";
	cout << "Initial energy: " << einit << "\tFinal energy: " << epot + ekin << "\tDifference(error): " << epot + ekin-einit << endl;
	cout << "Diffusion = " << diffusion << endl;
	output(2);
	return 0;
}

double instantaneoustemp()
{
	double sum=0;
	for (int i = 0; i < N; i++)
		for (int j = 0; j < 3; j++)
			sum+= vel[i][j] * vel[i][j];
	return sum / (3 * (N-1));
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
			pos[p][i] = pos[p][i] + vel[p][i] * dt;
		}
	}
	ekin /= N ;	//To make it per particle
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
			double dr2 = (dr[0]*dr[0] + dr[1]*dr[1] + dr[2]*dr[2]);
			if (dr2 < rc2)
			{
				//cout << r2i;
				//cin.ignore();
				double r2i = 1./dr2;
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
	epot /= N;
}

void calcekin()
{
	for (int p=0; p<N; p++)
		for (int i=0;i<3;i++)
			ekin += 0.5 * vel[p][i] * vel[p][i];
	ekin /= N;
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


void printdata(bool thermalize,double time)
{
	ofstream file;
	if (thermalize)
	{
		//sprintf(filename,"data/%g-%g-%d/thermdata.dat",(double)T,(double)rho,n);
		sprintf(filename,"%sthermdata.dat",folder);
		file.open(filename,fstream::app);
		file << ekin + epot << " " << ekin << " " << epot << endl;
	}
	else 
	{
		//sprintf(filename,"data/%g-%g-%d/data.dat",(double)T,(double)rho,n);
		sprintf(filename,"%sdata.dat",folder);
		file.open(filename,fstream::app);

		//Calculate diffusion
		diffusion=0;
		for (int i=0;i<N;i++)
			for (int j=0;j<3;j++)
				diffusion+=(pos[i][j] - posinit[i][j]) * (pos[i][j] - posinit[i][j]);
		diffusion/= 6. * N * time;

		file << ekin + epot << " " << ekin << " " << epot << " " << diffusion << endl;
	}
	file.close();
}

void printimage()
{
	static int filecount=1;
	char filename[100];
	sprintf(filename,"%simages/fcc-%04d.dat",folder,filecount);
	ofstream file(filename);
	file << N << endl;
	file << "0 " << L << endl;
	file << "0 " << L << endl;
	file << "0 " << L << endl;
	for (int i=0;i<N;i++)
		file << mod(pos[i][0],L) << " " << mod(pos[i][1],L) << " " << mod(pos[i][2],L) << " 1" << endl;
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
void output(int stage)
{
	//char filename[100];
	//sprintf(filename,"data/%g-%g-%d/output.txt",(double)T,(double)rho,n);
	sprintf(filename,"%soutput.txt",folder);
	ofstream file(filename,fstream::app);	
	if (stage ==0)
	{
	file << "\n\n\n\n";
	file << "Molecular Dynamics simulation of Argon\n\n";
	file << "Number of particles: " << N << endl;
	file << "Density: " << rho << endl;
	file << "Box length: " << L << endl;
	file << "Temperature: " << T << endl;
	file << "Thermalization iterations: " << thermiter << "\t measurement iterations" << iterations << endl;
	file << "dt: " << dt << "s" << endl;
	file << "Thermalization time: " << dt * thermiter << " seconds\n";
	file << "simulation time: " << dt * iterations << " seconds\n";
	file << "Initial energy: " << einit << "\n\n";
	}
	else if (stage ==1)
	{
		file << "Thermalization took " << double(clock() - begin) / CLOCKS_PER_SEC << " s\n";
		file << "Current energy: " << ekin+epot << "\tekin=" << ekin << "epot=" << epot << "\n\n";
	}
	else if (stage == 2)
	{
		file << "Simulation took " << double(clock() - begin) / CLOCKS_PER_SEC << " s\n";
		file << "Initial energy: " << einit << "\tFinal energy: " << epot + ekin << "\tDifference(error): " << epot + ekin-einit << endl;
		file << "Kinetic energy: " << ekin << "\tPotential energy: " << epot << endl;
		file << "Diffusion = " << diffusion << endl;

	}

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
		double u = sqrt(T) * sqrt(-2*log(x))*cos(2*pi*y);
		second = sqrt(T) * sqrt(-2*log(x))*sin(2*pi*y);
		extra = 1;
		return u;
	}
}

double mod(double num,double div)
{
	return num - floor(num/div) * div;
}



