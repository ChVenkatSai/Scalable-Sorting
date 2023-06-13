


#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/time.h>
#include "mpi.h"

#define SIZE 10



static int compare(const void *i, const void *j)
{
  if ((*(int *)i) > (*(int *)j))
    return (1);
  if ((*(int *)i) < (*(int *)j))
    return (-1);
  return (0);
}



int main (int argc, char *argv[])
{


  int 	     processors,MyRank, Root = 0;
  int 	     i,j,k, elements, elements_Bloc,
				  NoElementsToSort;
  int 	     count, temp;
  int 	     *Input, *Datatoinput;
  int 	     *Splitter, *AllSplitter;
  int 	     *Buckets, *BucketBuffer, *LocalBucket;
  int 	     *OutputBuffer, *Output;
  MPI_Status  status; 

  
  MPI_Init(&argc, &argv);
  MPI_Comm_size(MPI_COMM_WORLD, &processors);
  MPI_Comm_rank(MPI_COMM_WORLD, &MyRank);

  

  
  if (MyRank == Root){

    elements = atoi(argv[1]);
    Input = (int *) malloc (elements*sizeof(int));
	

 
    srand48((unsigned int)elements);
	 for(i=0; i< elements; i++) {
       Input[i] = rand()%1000 +1;
       
    }
  }

 
    struct timeval tv1;
    gettimeofday(&tv1, NULL);
      
      

  MPI_Bcast (&elements, 1, MPI_INT, 0, MPI_COMM_WORLD);
  

  elements_Bloc = elements / processors;
  Datatoinput = (int *) malloc (elements_Bloc * sizeof (int));

  MPI_Scatter(Input, elements_Bloc, MPI_INT, Datatoinput, 
				  elements_Bloc, MPI_INT, Root, MPI_COMM_WORLD);

  //for (int i = 0; i < n ; i++) {
    //    cout << local_list[i] << " ";
    //}
    //cout <<"My rank" << my_rank << endl;
  qsort ((char *) Datatoinput, elements_Bloc, sizeof(int), compare);

  Splitter = (int *) malloc (sizeof (int) * (processors-1));
  for (i=0; i< (processors-1); i++){
        Splitter[i] = Datatoinput[elements/(processors*processors) * (i+1)];
  } 

  
  AllSplitter = (int *) malloc (sizeof (int) * processors * (processors-1));
  MPI_Gather (Splitter, processors-1, MPI_INT, AllSplitter, processors-1, 
				  MPI_INT, Root, MPI_COMM_WORLD);

  
  if (MyRank == Root){
    qsort ((char *) AllSplitter, processors*(processors-1), sizeof(int), compare);

    for (i=0; i<processors-1; i++)
      Splitter[i] = AllSplitter[(processors-1)*(i+1)];
  }
  
  MPI_Bcast (Splitter, processors-1, MPI_INT, 0, MPI_COMM_WORLD);

  
  Buckets = (int *) malloc (sizeof (int) * (elements + processors));
  
  j = 0;
  k = 1;

  for (i=0; i<elements_Bloc; i++){
    if(j < (processors-1)){
       if (Datatoinput[i] < Splitter[j]) 
			 Buckets[((elements_Bloc + 1) * j) + k++] = Datatoinput[i]; 
       else{
	       Buckets[(elements_Bloc + 1) * j] = k-1;
		    k=1;
			 j++;
		    i--;
       }
    }
    else 
       Buckets[((elements_Bloc + 1) * j) + k++] = Datatoinput[i];
  }
  Buckets[(elements_Bloc + 1) * j] = k - 1;
      
  

  BucketBuffer = (int *) malloc (sizeof (int) * (elements + processors));

  MPI_Alltoall (Buckets, elements_Bloc + 1, MPI_INT, BucketBuffer, 
					 elements_Bloc + 1, MPI_INT, MPI_COMM_WORLD);

 
  LocalBucket = (int *) malloc (sizeof (int) * 2 * elements / processors);

  count = 1;

  for (j=0; j<processors; j++) {
  k = 1;
    for (i=0; i<BucketBuffer[(elements/processors + 1) * j]; i++) 
      LocalBucket[count++] = BucketBuffer[(elements/processors + 1) * j + k++];
  }
  LocalBucket[0] = count-1;
  // Recursively sort the less and greater lists
  //quicksort(less, my_rank, comm_size);
  //quicksort(greater, my_rank, comm_size);
  // Copy the sorted elements back into the local list
  //copy(less.begin(), less.end(), local_list.begin());
  //copy(greater.begin(), greater.end(), local_list.begin() + less.size());
  //cout << "I'm done" << my_rank << endl;
  

  NoElementsToSort = LocalBucket[0];
  qsort ((char *) &LocalBucket[1], NoElementsToSort, sizeof(int), compare); 

  
  if(MyRank == Root) {
  		OutputBuffer = (int *) malloc (sizeof(int) * 2 * elements);
  		Output = (int *) malloc (sizeof (int) * elements);
  }

  MPI_Gather (LocalBucket, 2*elements_Bloc, MPI_INT, OutputBuffer, 
				  2*elements_Bloc, MPI_INT, Root, MPI_COMM_WORLD);

  //cout << "Im done quickly" << my_rank<<endl;
     //MPI_Gather(&local_list[0], n, MPI_INT, &sorted_list[0], n, MPI_INT, 0, MPI_COMM_WORLD);
	if (MyRank == Root){
		count = 0;
		for(j=0; j<processors; j++){
          k = 1;
      	 for(i=0; i<OutputBuffer[(2 * elements/processors) * j]; i++) 
				 Output[count++] = OutputBuffer[(2*elements/processors) * j + k++];
    	}
        struct timeval tv2;
        gettimeofday(&tv2, NULL);
        long elapsed = (tv2.tv_sec - tv1.tv_sec) * 1000000 + tv2.tv_usec - tv1.tv_usec;
        printf("Time elapsed: %ld microseconds\n", elapsed);

        int elements2 = atoi(argv[1]);
        int* Input2 = (int*)malloc(elements2 * sizeof(int));
        srand48((unsigned int)elements2);
        for (i = 0; i < elements2; i++) {
            Input2[i] = rand() % 1000 + 1;
            
        }
        struct timeval t2,t1;
        gettimeofday(&t1, NULL);
        qsort((char*)Input2, elements2, sizeof(int), compare);
        gettimeofday(&t2, NULL);
        long elapsed2 = (t2.tv_sec - t1.tv_sec) * 1000000 + t2.tv_usec - t1.tv_usec;
        printf("Time elapsed: %ld microseconds\n", elapsed2);
        double speedup = (double)elapsed2 / (double)elapsed;
        double efficiency = speedup / (double)processors;
        printf("Efficiency: %f\n", efficiency);
      
		 
    	
	
   }

  	

   MPI_Finalize();
}


