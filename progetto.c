#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <mpi.h>
#include <assert.h>


#define SOFTENING 1e-9f

typedef struct {
  float x, y, z, vx, vy, vz;
} Body;

void create_bodies_file(FILE *fp, int n, Body *body);
void create_sendcnts_displs(int *sendcnts,int *displs,int n_bodies,int world_size);
void body_force(Body *body, float dt, int n_bodies, int start, int end);

int main(int argc, char** argv){

  	const float dt = 0.01f;
  

	MPI_Status status;
	MPI_Datatype body_type;	

	int rank;					//rank del processo 
	int world_size;					//numero di processi 
	int count;					//numero dei campi float da settare 
  	int n_bodies;					// numero particelle
	int bytes;					// grandezza della struttura
	int *sendcnts;					// numero elementi per processo
	int *displs;					// puntatore alla posizione iniziale della partizione
	int it;						//indice del for 
	int n_iters;					// numero d iterazioni della simulaizone 

	float *buffer;					// buffer per la memoria 

	double time_start;				//inizio tempo di esecuzione 
	double time_end;				// fine tempo di esecuzione
	double avgTime;					// media tempo di esecuzione 

	FILE *out_file;					//puntatore al file da srivere
	FILE *in_file;					//puntatore al file da leggere 

	MPI_Init(&argc, &argv);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank); 					
  	MPI_Comm_size(MPI_COMM_WORLD, &world_size);

	count = 6;
	int array_of_blocklengths[6] = {1, 1, 1, 1, 1, 1};

	// dichiarazione dei tipi dei datatype 

	MPI_Datatype array_of_types[6] = {MPI_FLOAT, MPI_FLOAT, MPI_FLOAT, MPI_FLOAT, MPI_FLOAT, MPI_FLOAT};
	
	// inizializzazione dell'array che contiene lo spostamento per ogni blocco espresso in byte  
	MPI_Aint array_of_displacements[6];
  	array_of_displacements[0] = sizeof(float) * 0;
  	array_of_displacements[1] = sizeof(float) * 1;
  	array_of_displacements[2] = sizeof(float) * 2;
  	array_of_displacements[3] = sizeof(float) * 3;
  	array_of_displacements[4] = sizeof(float) * 4;
  	array_of_displacements[5] = sizeof(float) * 5;
	
	//creazione del datatype 
	 MPI_Type_create_struct(count, array_of_blocklengths, array_of_displacements, array_of_types, &body_type);
  
	MPI_Type_commit(&body_type);

	//controllo degli argomenti in input 
	if (argc == 2){
		 n_bodies = atoi(argv[1]);
		  n_iters = 5;
	} else if (arc == 3 ){
			n_bodies = atoi(argv[1]);
			n_iters = atoi(argv[2]);
  	}else {
         n_bodies = 10000;
				 n_iters = 5;
	}

	bytes = n_bodies * sizeof(Body);
	buffer = (float*)malloc(bytes);
  	assert(buffer != NULL);
  	Body *body = (Body*)buffer;

	if(argc == 2 || argc == 3){
		out_file = fopen("file.txt","w");
	if(out_file == NULL){
		printf("Errore file\n");
      		exit(-1);
    	}
	create_bodies_file(out_file, n_bodies, body);
	fclose(out_file);
	} else if (argc > 3) {
		in_file = fopen("file.txt","r");
	if(in_file == NULL){
      		printf("Errore file\n");
     		exit(-1);
    	}
    	fread(body, sizeof(Body), n_bodies, in_file);
    	fclose(in_file);
  	}

 
  	bytes = sizeof(int) * world_size;
  	sendcnts = (int*)malloc(sizeof(int)*world_size);
  	assert(sendcnts != NULL);
  	displs = (int*)malloc(sizeof(int)*world_size);
  	assert(displs != NULL);
  	create_sendcnts_displs(sendcnts, displs, n_bodies, world_size);

 
	time_start = MPI_Wtime();

	for(it = 0; it < n_iters; it++){

	    body_force(body, dt, n_bodies, displs[rank], (displs[rank] + sendcnts[rank]));

	for(int i = sendcnts[rank]; i < displs[rank] + sendcnts[rank]; i++){
      		//Integrate position
      		body[i].x += body[i].vx*dt;
      		body[i].y += body[i].vy*dt;
      		body[i].z += body[i].vz*dt;
    	}
    
	MPI_Allgatherv(MPI_IN_PLACE, sendcnts[rank], body_type, body, sendcnts, displs, body_type, MPI_COMM_WORLD);


	}

	MPI_Barrier(MPI_COMM_WORLD);

	time_end = MPI_Wtime();

	avgTime = (time_end - time_start)/n_iters;

	if(rank == 0){
		printf("\nTempo %.4f\n", time_end-time_start);
  		printf("\n Avg Time = %.4f\n",avgTime);
  	}

	MPI_Type_free(&body_type);
  	free(sendcnts);
  	free(displs);
  	free(buffer);

	MPI_Finalize();

	return EXIT_SUCCESS;
}


void create_bodies_file(FILE *fp, int n, Body *body){
  for(int i = 0; i < n; i++){
    body[i].x = 2.0f * (rand() / (float)RAND_MAX) - 1.0f;
    body[i].y = 2.0f * (rand() / (float)RAND_MAX) - 1.0f;
    body[i].z = 2.0f * (rand() / (float)RAND_MAX) - 1.0f;
    body[i].vx = 2.0f * (rand() / (float)RAND_MAX) - 1.0f;
    body[i].vy = 2.0f * (rand() / (float)RAND_MAX) - 1.0f;
    body[i].vz = 2.0f * (rand() / (float)RAND_MAX) - 1.0f;
  }
  fwrite(body, sizeof(Body) * n, 1, fp);
}

void create_sendcnts_displs(int *sendcnts,int *displs,int n_bodies,int world_size){
  int resto = n_bodies % world_size;
  int last = 0;

  for(int i=0; i<world_size; i++){
    sendcnts[i] = n_bodies / world_size;
    if(resto > 0){
      resto--;
      sendcnts[i] = sendcnts[i] + 1;
    }
    displs[i] = last;
    last = last + sendcnts[i];
  }
}

void body_force(Body *body, float dt, int n_bodies, int start, int end){
    for(int i = start; i < end; i++){
      float Fx = 0.0f; float Fy = 0.0f; float Fz = 0.0f;

      for(int j = 0; j < n_bodies; j++){
        float dx = body[j].x - body[i].x;
        float dy = body[j].y - body[i].y;
        float dz = body[j].z - body[i].z;
        float distSqr = dx*dx + dy*dy + dz*dz + SOFTENING;
        float invDist = 1.0f / sqrtf(distSqr);
        float invDist3 = invDist * invDist * invDist;

        Fx += dx * invDist3; Fy += dy * invDist3; Fz += dz * invDist3;
      }
      body[i].vx += dt*Fx;
      body[i].vy += dt*Fy;
      body[i].vz += dt*Fz;
    }

}
