#include <fstream>
#include <iostream>
#include <string>
#include <cmath>

/*param filename :the name of the file to open
param rows       :number of rows in the file data
param cols       :number of columns in the files data
return float**   :a 2d float array holding the read in values
			or NULL on failure.*/
float** read_file(std::string filename,int rows,int cols)
{
	std::fstream file;//create a stream for the file
	file.open(filename.c_str(), std::ios::in);//open the file to read in
	if(!file.is_open()){return 0;}//if the file failed to open return NULL

        //float** is just a 2d array like float values[cols][rows]
	//create the column pointers
	float** floats = new float*[cols+1];
	
	//create the row pointers
	for(int i = 0; i <cols;++i){ floats[i] = new float[rows+1]; }

	//read each through each row
	for(int i =0;i<rows;++i)
	{
		//read the values in this row and push into the correct
		//column.floats is [cols][rows]
		for(int j =0;j<cols;++j)//push into the col
		{ file >>floats[j][i]; }
	}
	file.close();//close the file

	return floats;//return the 2d array
}

int main()
{
  int rows = int(114);//number of rows in the file
  int cols = 2;//number of columns in the file
  float** floats;//2d array
  
  /*the func read_file returns null on failure, so set floats to the 
    value returned (which should be a 2d array) and check for failure.
    if the function failed exit main.*/ 
  if( !(floats = read_file("dos",rows,cols) ) ){return 0;}
  double prom[800],sum2,count,nmcs,sigma;
  sum2=0;
  count=0;
  sigma=0;
  double a;
  int t;
  for (float T =0.01; T <= 8.0; T=T+0.01)
  {   
      t=0;
      for(int i = 0; i < 114; i++)
	{
	  a=((floats[1][i]))*exp(-(floats[0][i])/T);
	  count+=a;
	  sum2+=(floats[0][i])*a;
	  //a=floats[1][i]*pow(2,-int(floats[1][i]));
	  //std::cout<<floats[0][i]<<"\t"<<a<<"\n";
	}
      sigma=-T*log(count);
      prom[t]=sum2/count;
      std::cout<<T<<"\t"<<prom[t]<<"\t"<<sigma<<"\n";
      t++;
      }


  
  return 0;
}
