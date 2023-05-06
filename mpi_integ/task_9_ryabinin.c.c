//the program calculates the integral sin(x) from 1 to 2
#include <stdio.h>
#include <math.h>
#include <mpi.h>

int main(int argc, char* argv[])
{
	double a = 1, b = 2; //integration boundaries
	double start, end; //the variables for the timer
	int N = 600000000;//number of steps
	int size, rank;//number of processes and process number
	double sum = 0;//variable for the result
	double partialAmount = 0;// partial amount for each process
	double x;

	double n = 1.0 / (double)N;//step length

	MPI_Status status;

	MPI_Init(&argc, &argv);

	MPI_Comm_size(MPI_COMM_WORLD, &size);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);

	double partialSegment = (double)rank / (double)size;//the starting point for the process
	double partiakSegmentNext = (double)(rank + 1) / (double)size;//the starting point for the process with the number n+1



	if (rank == 0)
	{
		//THE START RECTANGLE METHOD
		start = MPI_Wtime();
		double s; //partial amount for 0 process
		double t;
		for (double i = a + partialSegment; i < a + partiakSegmentNext; i += n) //rectangle method
		{
			partialAmount += sin(i);
		}
		partialAmount *= n;
		sum += partialAmount;

		if (size != 1)
		{
			for (int i = 1; i < size; i++)//receiving messages
			{
				MPI_Recv(&s, 1, MPI_DOUBLE, i, 1, MPI_COMM_WORLD, &status);
				sum += s;
			}
		}

		if (N % size != 0)
		{
			MPI_Recv(&t, 1, MPI_DOUBLE, size - 1, 4, MPI_COMM_WORLD, &status);
			sum += t;
		}

		end = MPI_Wtime();
		printf("sum rectangle = %lf\n", sum);
		printf("time = %lf\n", end - start);


		//THE START TRAPEZIOID METHOD
		sum = 0;
		partialAmount = 0;

		start = MPI_Wtime();
		for (double i = a + partialSegment; i < a + partiakSegmentNext; i += n)
		{
			partialAmount += (sin(i) + sin(i + n));
		}
		partialAmount *= (n / 2);
		sum += partialAmount;

		if (size != 1)
		{
			for (int i = 1; i < size; i++)
			{
				MPI_Recv(&s, 1, MPI_DOUBLE, i, 2, MPI_COMM_WORLD, &status);
				sum += s;
			}
		}

		if (N % size != 0)
		{
			MPI_Recv(&t, 1, MPI_DOUBLE, size - 1, 5, MPI_COMM_WORLD, &status);
			sum += t;
		}
		end = MPI_Wtime();
		printf("sum trapezioid = %lf\n", sum);
		printf("time = %lf\n", end - start);


		//THE START SIMPSON METHOD
		sum = 0;
		partialAmount = 0;

		int curr = 0;
		start = MPI_Wtime();
		for (double i = a + partialSegment; i < a + partiakSegmentNext; i += n)
		{
			int coefficient = 2 + 2 * (curr % 2);
			partialAmount += (coefficient * sin(i + n));
			curr++;
		}
		partialAmount *= (n / 3);
		sum += partialAmount;
		if (size != 1)
		{
			for (int i = 1; i < size; i++)
			{
				MPI_Recv(&s, 1, MPI_DOUBLE, i, 3, MPI_COMM_WORLD, &status);
				sum += s;
			}
		}

		if (N % size != 0)
		{
			MPI_Recv(&t, 1, MPI_DOUBLE, size - 1, 6, MPI_COMM_WORLD, &status);
			sum += t;
		}
		end = MPI_Wtime();
		printf("sum simpson = %lf\n", sum);
		printf("time = %lf\n", end - start);

	}
	else
	{
		//THE START RECTANGLE METHOD
		for (double i = a + partialSegment; i < a + partiakSegmentNext; i += n)
		{
			partialAmount += sin(i);
		}
		partialAmount *= n;

		MPI_Send(&partialAmount, 1, MPI_DOUBLE, 0, 1, MPI_COMM_WORLD);



		//THE START TRAPEZIOID METHOD
		partialAmount = 0;

		for (double i = a + partialSegment; i < a + partiakSegmentNext; i += n)
		{
			partialAmount += (sin(i) + sin(i + n));
		}
		partialAmount *= (n / 2);
		MPI_Send(&partialAmount, 1, MPI_DOUBLE, 0, 2, MPI_COMM_WORLD);



		//THE START SIMPSON METHOD
		partialAmount = 0;

		int curr = 0;
		for (double i = a + partialSegment; i < a + partiakSegmentNext; i += n)
		{

			int coefficient = 2 + 2 * (curr % 2);
			partialAmount += (coefficient * sin(i + n));
			curr++;
		}
		partialAmount *= (n / 3);
		MPI_Send(&partialAmount, 1, MPI_DOUBLE, 0, 3, MPI_COMM_WORLD);

		partialAmount = 0;
	}
	if ((N % size != 0) && (rank = (size - 1)))
	{
		for (int i = 0; i < N % size; i++)
		{
			x = i * n + partialSegment;
			partialAmount += sin(x);
		}
		partialAmount *= n;
		MPI_Send(&partialAmount, 1, MPI_DOUBLE, 0, 4, MPI_COMM_WORLD);

		partialAmount = 0;

		for (int i = 0; i < N % size; i++)
		{
			x = i * n + partialSegment;
			partialAmount += sin(x) + sin(x + n);
		}
		partialAmount *= n / 2;
		MPI_Send(&partialAmount, 1, MPI_DOUBLE, 0, 5, MPI_COMM_WORLD);
		partialAmount = 0;
		for (int i = 0; i < N % size; i++)
		{
			x = i * n + partialSegment;
			int coefficient = 2 + 2 * (i % 2);
			partialAmount += (coefficient * sin(x + n));
		}
		partialAmount *= n;
		MPI_Send(&partialAmount, 1, MPI_DOUBLE, 0, 6, MPI_COMM_WORLD);

	}

	MPI_Finalize();

	return 0;
}

//C:\Users\mv\source\repos\mpi_integ\Debug