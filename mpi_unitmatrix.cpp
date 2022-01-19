#include <iostream>
#include <cmath>
#include <algorithm>
#include <cstdlib>
#include "mpi.h"

void f_matrix(double * & matrix,int nrows, int ncols, int pid, int np);
void p_matrix(const double * & matrix,int nrows,int ncols, int pid, int np);
void p_slab(const double * & matrix,int nrows, int ncols);

int main(int argc, char **argv[])
{
	MPI_Init(&argc &argv);
	int pid,np;
	MPI_Comm_rank(MPI_COMM_WORLD,&pid);
	MPI_Comm_size(MPI_COMM_WORLD,&np);

	int N=std::atoi(argv[1]);
	int Nlo= N/np;

	double * matrix = new double[Nlo*N]{0.0};
	f_matrix(matrix,Nlo,N,pid,np);
	p_matrix(matrix,Nlo,N,pid,np);
	delete [] matrix; 
	matrix =nullptr;

	MPI_Finalize();
		return 0;
	}


void p_slab(const double * matrix,int nrows, int ncols){
	for (int i = 0; i < nrows; ++i)
	{
		for (int j = 0; j < ncols; ++j)
		{
			std::cout << matrix[i*ncols+j]<<" ";
		}
		std::cout<<"\n";
	}
}

void p_matrix(const double * & matrix,int nrows,int ncols, int pid, int np){
	if (0==pid)
	{
		p_slab(matrix,nrows,ncols);
		double * buffer = new double[nrows*ncols]{0.0};
		for (int scr = 1; scr < np; ++scr)
		{
			MPI_Recv(&buffer[0], nrows*ncols, MPI_DOUBLE,scr,tag, MPI_COMM_WORLD,MPI_STATUS_IGNORE);
			p_slab(buffer,nrows,ncols);
		}
		delete[] buffer;
		buffer=nullptr;
	}
	else{
		MPI_Send(&matrix[0],nrows*ncols, MPI_DOUBLE,0,tag, MPI_COMM_WORLD);
	}
}

void f_matrix(double * matrix, int nrows,int ncols, int pid, int np){
	for (int i = 0; i < nrows; ++i)
	{
		int j=nrows*pid +i;
		matrix[i*ncols+j]=1.0;
	}
}