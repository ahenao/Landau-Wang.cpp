/*
This Program perform a single spin flip Glauber Dynamics simulation of the Ising
Magnetic system, by the technique of Monte Carlo

Andres Henao Aristizabal
Physics Engineer- Universidad Nacional de Colombia
Studying Master in Computational Physics- Universitat Politecnica Catalunya, Unversitat Barcelona
ahenaoa@unal.edu.co 
*/




#include <iostream>
#include <fstream> //Contiene la función para escribir ficheros .txt
#include <cmath>
#include <ctime>                      // define time()
#include <cstring>
#include <cstdlib>

#include "randomc.h"    //Library with the Random Number Generators

/*To compile put in the same directory of the .cpp the library randomc.h, mersenne.o, userintf.o
and write: 
g++ 2DIsing_principal.cpp -o 2DIsing  mersenne.o usernintf.o

The codes mersenne.cpp and userintf.cpp are constructed as .o to avoid the warnings on the differences between C and C++*/

using namespace std;

void usage ( void )
{
  cout<<"\n--------------------------------------------------------------------------------"<<endl;
  cout<<"\t\t\tMonte Carlo Simulation of\n\t\t\t     the Ising Model\n";
  cout<<"\n\t\t\t Written by Andres Henao";
  cout<<"\n\t       Master in Computational and Applied Physics\n\t\t\t    UPC-UB June 2011\n\n";
  cout<<"./2DIsing usage:\n";
  cout<<">> 2DIsing [options]\n\n";
  cout<<"Options:\n";
  cout<<"\t-d\t[integer]\tDimension\n";
  cout<<"\t-L\t[integer]\tLinear Size\n";
  cout<<"\t-seed\t[integer]\tSeed of the random number generator\n";
  cout<<"\t-h  \t\t\tPrint this info\n";
  cout<<"--------------------------------------------------------------------------------\n"<<endl;
}

/********************Initialization of Variables****************/
int qq,kk;//The random spin to change
float r,p,J=1,h=0;//Random number, probability, exchange constant, external field
double Hd,delta;//To measure delta of energy
int d=2;
double L=20;
float T=2.0;
double nmcs=1e5;
int nmeas=1;
int seed=0;
double N;
int n;
double energy=0,energy_trial=0,magnet=0;
int meas_ints=217;//
double dos[218];//This is for each subinterval of energy
double hist[218];//
float f;
int flat;
double entries=0;//will tell us the number of entries in histogram

ofstream a11("g1");  //Put here the name of the densitiy of states
ofstream a12("H1");  //Put here the name of the histogram
/****************************************************************/


/******Here we parse the command line arguments; If you add an option, document it in the usage() function!****/
void read(int argc, char *argv[])  
{
  for (int i = 1; i < argc; i++) {
         if (!strcmp(argv[i],"-L")) L = atof(argv[++i]);
	 else if (!strcmp(argv[i],"-seed")) seed = atoi(argv[++i]);
	 else if (!strcmp(argv[i],"-h")) 
      {
      usage(); 
      exit(0);
      }
    else {
      fprintf(stderr,"Error: Command-line argument '%s' not recognized.\n",
	      argv[i]);
      exit(-1);
    }
  }
}  
/****************************************************************/
CRandomMersenne RanGen(seed);       // make instance of random number generator

/****************Initialize the Lattice***************************/
void initialize_square(int *spinij[],int n)
{
  for(int i=1;i<=n;i++)        //Initial Configuration
    {   
      for(int j=1;j<=n;j++)
	{
	  r= RanGen.Random();
	  if(r<0.5)
	    {
	      spinij[i][j]=1;
	      magnet+=1;
	    }
	  else
	    {
	      spinij[i][j]=-1;
	      magnet-=1;
	    }
	  if(i==1){spinij[n+1][j]=spinij[1][j];}//Boundary neighbours
	  if(i==n){spinij[0][j]=spinij[n][j];}
	  if(j==1){spinij[i][n+1]=spinij[i][1];}
	  if(j==n){spinij[i][0]=spinij[i][n];}
	}
    }
  for(int i = 1; i <= n; i++)       //Measuring Initial energy
    {
      for(int j = 1; j <= n; j++)
	{
	  energy=energy-h*spinij[i][j];
	  if(spinij[i][j]==spinij[i-1][j])
	    {energy=energy-J;}
	  else
	    {energy=energy+J;}
	  if(spinij[i][j]==spinij[i+1][j])
	    {energy=energy-J;}
	  else
	    {energy=energy+J;}
	  if(spinij[i][j]==spinij[i][j-1])
	    {energy=energy-J;}
	  else
	    {energy=energy+J;}
	  if(spinij[i][j]==spinij[i][j+1])
	    {energy=energy-J;}
	  else
	    {energy=energy+J;}	                     
	}
    }
  energy*=0.5;
  energy=energy/(n*n);

  //------------Set all the g(E)=1, or the same ln(g(E))=0----------//
  for(int i = 0; i < meas_ints; i++)
    {
      dos[i]=1;
      hist[i]=0;
    }
  //---------------------------------------------------------------//
}
/********************************************************************/


/********************Monte Carlo Glauber and Mteropolis Update******************************/
void Monte_Carlo(int *spinij[],int n,double N, int inter, float f)
{     
  for(int m = 1; m <= n; m++)//N attempted spin flips
    {
      for(int o = 1; o <= n; o++)
	{
  //------------------Chosing a random spin to flip--------------------------//     
	   qq=RanGen.IRandom(1,n);
	  kk=RanGen.IRandom(1,n);//Vamos a cambiar el spin qq,kk
	  //qq=m;
	  //kk=o;
      
      
      //--------------------Finding the new Hamiltonian--------------------------//
      Hd=0;delta=0;
      Hd=spinij[qq-1][kk]+spinij[qq+1][kk]+spinij[qq][kk-1]+spinij[qq][kk+1];//*Uncomment for Single spin flip
      delta=2*Hd*spinij[qq][kk]*J + 2*h*spinij[qq][kk];
      energy_trial=energy+delta/N;
      int bin1=0,bin2=0,in1=0,in2=0;//bin1 and bin2 choose the interval and in tells if the energy is in the interval
                       //otherwise we do nothing
      float change;
      change = -2.03+0.3*(inter-1);
      for(int i=0; i < meas_ints-1; i++)
       	{
      	  if(energy>=(change+i/100.0)&&energy<(change+(i+1)/100.0))
      	    {
      	    bin1=i;
      	    in1=1;
      	    }
      	  change = -2.03+0.3*(inter-1);
      	  if(energy_trial>(change+i/100.0)&&energy_trial<(change+(i+1)/100.0))
      	    {
      	    bin2=i;
      	    in2=1;
      	    }
      	}


      if(in1*in2==1)//both the energy and energy trial are within the interval
      	{
      	  if(bin1!=bin2)
      	    {
      	      if(log(dos[bin2])>=log(dos[bin1]))
      		{
      		  p=exp(log(dos[bin1])-log(dos[bin2]));
      		  r= RanGen.Random();
      		  if(r<=p)
      		    {
      		      spinij[qq][kk]=-1*spinij[qq][kk];
      		      magnet+=spinij[qq][kk];
      		      energy=energy_trial;
      		      if(qq==1&&kk!=1&&kk!=n){spinij[n+1][kk]=-spinij[n+1][kk];}
      		      if(qq==n&&kk!=1&&kk!=n){spinij[0][kk]=-spinij[0][kk];}
      		      if(qq==1&&kk==1){spinij[n+1][n+1]=-spinij[n+1][n+1];}
      		      if(qq==n&&kk==n){spinij[0][0]=-spinij[0][0];}
      		      if(qq==1&&kk==n){spinij[n+1][0]=-spinij[n+1][0];}
      		      if(qq==n&&kk==1){spinij[0][n+1]=-spinij[0][n+1];}
      		      if(kk==1&&qq!=1&&qq!=n){spinij[qq][n+1]=-spinij[qq][n+1];}
      		      if(kk==n&&qq!=1&&qq!=n){spinij[qq][0]=-spinij[qq][0];}
      		      dos[bin2]=dos[bin2]+log(f);
      		      hist[bin2]++;
      		      entries++;
      		    } 
		  else
		    {
		      dos[bin1]=dos[bin1]+log(f);
		      hist[bin1]++;
		      entries++;
		    }
      		}
	      else
		{
		  spinij[qq][kk]=-1*spinij[qq][kk];
		  magnet+=spinij[qq][kk];
		  energy=energy_trial;
		  if(qq==1&&kk!=1&&kk!=n){spinij[n+1][kk]=-spinij[n+1][kk];}
		  if(qq==n&&kk!=1&&kk!=n){spinij[0][kk]=-spinij[0][kk];}
		  if(qq==1&&kk==1){spinij[n+1][n+1]=-spinij[n+1][n+1];}
		  if(qq==n&&kk==n){spinij[0][0]=-spinij[0][0];}
		  if(qq==1&&kk==n){spinij[n+1][0]=-spinij[n+1][0];}
		  if(qq==n&&kk==1){spinij[0][n+1]=-spinij[0][n+1];}
		  if(kk==1&&qq!=1&&qq!=n){spinij[qq][n+1]=-spinij[qq][n+1];}
		  if(kk==n&&qq!=1&&qq!=n){spinij[qq][0]=-spinij[qq][0];}
		  dos[bin2]=dos[bin2]+log(f);
		  hist[bin2]++;
		  entries++;
		}
      	    } 
      	}
      else if(in1==1&&in2==0){ dos[bin1]=dos[bin1]+log(f);}//The trial state is not in the energy interval
      else if(in1==0&&in2==1)//This is for the first time that the energy enters to the interval
      	{
      	  spinij[qq][kk]=-1*spinij[qq][kk];
      	  energy=energy_trial;
      	  if(qq==1&&kk!=1&&kk!=n){spinij[n+1][kk]=-spinij[n+1][kk];}
      	  if(qq==n&&kk!=1&&kk!=n){spinij[0][kk]=-spinij[0][kk];}
      	  if(qq==1&&kk==1){spinij[n+1][n+1]=-spinij[n+1][n+1];}
      	  if(qq==n&&kk==n){spinij[0][0]=-spinij[0][0];}
      	  if(qq==1&&kk==n){spinij[n+1][0]=-spinij[n+1][0];}
      	  if(qq==n&&kk==1){spinij[0][n+1]=-spinij[0][n+1];}
      	  if(kk==1&&qq!=1&&qq!=n){spinij[qq][n+1]=-spinij[qq][n+1];}
      	  if(kk==n&&qq!=1&&qq!=n){spinij[qq][0]=-spinij[qq][0];}
      	  dos[bin2]=dos[bin2]+log(f);
      	  hist[bin2]++;
      	  entries++;
      	}
    } 
    }
}
/*******************************************************************/


/********************checks the flatness************************/
float histogram()
{
  float flatt=0,error=1,err=0;
  double sum=0,mean=0;
  int i;
  for(i = 0; i < meas_ints; i++)
    {
      sum+=hist[i];
    }
  mean=sum/meas_ints;
  i=0;
  while(error==1&&i<meas_ints)
    {
      err=(hist[i]-mean)/mean;
      if(abs(err)<0.1){error=1;flatt=1;i++;}//Here you control the flat parameter
      else{error=0;flat=0;}//If any of the entries has an error bigger than you defined, the histogram is not flat
    }
  return flatt;
}
/********************************************************************/


/****************Set hist=0**********************************/
void set_hist()
{
 for(int i = 0; i < meas_ints; i++)
    {
      hist[i]=0;
    }
}
/*************************************************************/


/************************Output Files***************************/
void write(int inter)
{
  float ch;
  ch = -2.03+0.3*(inter-1);
  for(int i = 0; i < meas_ints; i++)
    {
      a11<<ch<<"\t"<<dos[i]<<endl;
      a12<<ch<<"\t"<<hist[i]<<endl;
      ch = ch+0.01;
    }
}
/**************************************************************/




int main( int argc, char * argv[] )
{
  read(argc,&argv[0]);//Reads the Input
  N=pow(L,d);         
  n=int(L);
  int interval;
  cout<<"Choose the interval:\n";
  cout<<"\t [1] [-2.00 , -1.67]\n";
  cout<<"\t [2] [-1.73 , -1.37]\n";
  cout<<"\t [3] [-1.43 , -1.07]\n";
  cout<<"\t [4] [-1.13 , -0.77]\n";
  cout<<"\t [5] [-0.83 , -0.47]\n";
  cout<<"\t [6] [-0.53 , -0.17]\n";
  cout<<"\t [7] [-0.23 ,  0.13]\n";
  cin>>interval;
  int **spinij;	
  spinij=new int*[n+2];
  for(int i = 0; i <= n+2; i++)
    {
      spinij[i]=new int[n+2];
    }
  initialize_square(spinij,n); //Initialize the Lattice and measures the Initial Conditions  
  f=1.2;
  nmcs=1;
  flat=0;
 
  while(f>1.0001)
    {
      while(flat==0&&nmcs<1e3)
      	{
	  int merse;
	  nmcs++;
	      Monte_Carlo(spinij,n,N,interval,f);//Glauber and Metropolis
	      
	      if(int(nmcs)%int(1000)==0)
		{
		  flat=histogram();
		  cout<<"flatting\t"<<flat<<endl;
		}
      	}
      f=sqrt(f);
      cout<<nmcs<<"\t"<<flat<<"\t"<<f<<"\t"<<entries<<"\t"<<energy<<endl;
      nmcs=1;
      flat=0;
      if(f>1.0001){
	set_hist();   }
    }
  write(interval);


  
  //*****Release the memory*********//
  for(int i = 0; i <= n; i++){
    delete [] spinij[i];}
  delete [] spinij;
  //********************************//
  
  return 0;
}


