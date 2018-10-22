#include	<stdio.h>
#include	<mpi.h>
#include        <stddef.h>
#include        <ctype.h>
#include        <stdlib.h>
#include        <math.h>
#include        <malloc.h>
#include        <errno.h>

#define sphere_dim   10      /* chosen number of dimensions of the sphere  */
#define per_processor_tries  1500000

#define NUM_PROC   4         /* chosen number of processors. */

void srand48();
double drand48();
double pow();
double sqrt();

/*-----------------------------*/
double compute_mean(int n, int *monte_carlo_array)
{
  int i;
  double temp_mean = 0.0;

  for (i = 1; i <= n; i++)
    {
      temp_mean += ((double) monte_carlo_array[i]);
    }
  temp_mean /= ((double) n);
  return ((double) temp_mean);
}

/*-----------------------------*/
double compute_std(int n, int *monte_carlo_array, int mean)
{
  int i; 
  double temp_stdev = 0.0;

  for (i = 1; i <= n; i++)
    {
      temp_stdev += pow((double) monte_carlo_array[i] - mean, 2.0);
    }
  temp_stdev /= ((double) n);
  return((double) sqrt(temp_stdev));
}


/*-----------------------------*/
/*-----------------------------*/
void main(int argc, char *argv[])
{
  int Procs;                   /* Number of processors */
  int my_rank;                 /* Processor number */
  int i, j;
  int *sphere_darts;      /* Based on MC probability, it stores (0,1) 
                                  depending on whether dart thrown 
                                  lands inside the circle or outside it ... */
  double volume;
  double pi;
  double mean;
  double stand_dev;
  double temp_stdev;
  double sphere_volume[NUM_PROC];    /* volume computed for all procs. */
  double stand_dev_array[NUM_PROC];    /* std computed for all procs. */
  double temp_sphere_radius = 0.0;
  MPI_Status status;
  
  MPI_Init(&argc, &argv);

  MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);

  MPI_Comm_size(MPI_COMM_WORLD, &Procs);

  srand48 ((long) my_rank);
  sphere_darts = 
    (int *) malloc((per_processor_tries + 1) * sizeof(int));

/*** Step [1]: Compute a single point. ***/

  for(i=1; i <= per_processor_tries; i++)
    {
      temp_sphere_radius = 0.0;

      for(j=1; j <= sphere_dim; j++)
	{
	  temp_sphere_radius += pow(drand48(), 2.0);
	}
      if(temp_sphere_radius <= 1.0)
	{
	  sphere_darts[i] = 1;
	}
      else
	{
	  sphere_darts[i] = 0;
	}
    }

/*** Step [2]: Compute the mean. ***/
  mean = 
    (double) compute_mean(per_processor_tries, sphere_darts);

/*** Step [3]: Compute the sphere's volume. ***/
  volume = 
    (double) mean * pow(2.0, sphere_dim);

/*** Step [4]: Compute the standard deviation. ***/
  stand_dev = 
    (double) compute_std(per_processor_tries, sphere_darts, mean);


/*** Step [5]: Send/Receive global data ***/
/*** Processor (0) acts as to coordinate results/outputs 
     coming from other processors.      ***/

  if(my_rank == 0)
    {
      sphere_volume[0] = volume;
      stand_dev_array[0] = stand_dev;
      for(i=1; i < Procs; i++)
       {
         MPI_Recv(&(sphere_volume[i]),1,MPI_DOUBLE,i,0,MPI_COMM_WORLD,&status);
         MPI_Recv(&(stand_dev_array[i]),1,MPI_DOUBLE,i,0,MPI_COMM_WORLD,&status);

	 printf("In processor: %d  sphere volume: %25.16e \n",i,sphere_volume[i]);
       }
   }
   else
     {
       MPI_Send(&volume,1,MPI_DOUBLE,0,0,MPI_COMM_WORLD);
       MPI_Send(&stand_dev,1,MPI_DOUBLE,0,0,MPI_COMM_WORLD);
     }

/**---[6]-----
  Processor (0) adds up the overall results sent from other processors....
  Get the final totals to be printed.... 
 ----------**/

  if(my_rank == 0)
    {
      temp_stdev = 0.0;
      volume = 0.0;
      for(i=0; i < Procs; i++)
	{
	  volume += sphere_volume[i];
	  temp_stdev += stand_dev_array[i];
	}
      volume /= ((double) Procs);
      stand_dev = (double) temp_stdev /((double) Procs);
      printf("======================Output======================\n");
      printf("The overal final values averaged over all processors: \n");
      printf("Volume of sphere of dimension %d  is  %25.16e \n",sphere_dim, volume);
      printf("Standard dev. is %25.16e \n", stand_dev);
      printf("Total number of points is %d \n", 4*per_processor_tries);
      if(sphere_dim == 2)
	printf("In dimension 2, Volume of the Sphere = Pi = %25.16e \n", volume);
      printf("==================================================\n");
    }
  MPI_Finalize();
}

