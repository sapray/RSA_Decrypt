/* nonb.c */
#include <stdio.h>
#include "mpi.h"
#include <gmp.h>

int main (int argc, char *argv[])
{
  // init
  double start_time, end_time, run_time;
  int p, my_rank, processRight, processLeft, buf, i;
  MPI_Request reqs;
  MPI_Status stats;

  int flag, count = 0;

  mpz_t n, rootn, segmentCount, segmentSize, temp;
  mpz_init(n); mpz_init(rootn); mpz_init(segmentCount); mpz_init(segmentSize); mpz_init(temp);

  mpz_set_str(n, argv[1], 10);
  mpz_sqrt(rootn, n);
  
  MPI_Init(&argc, &argv);
  MPI_Comm_size(MPI_COMM_WORLD, &p);
  MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);

  // calculate chunk size
  mpz_set_ui(segmentCount, p*p);
  mpz_cdiv_q(segmentSize, rootn, segmentCount);

  // set processRight and processLeft values
  processLeft = my_rank-1;   
  processRight = my_rank+1;
  if (my_rank == 0)
    processLeft = p - 1;
  if (my_rank == (p - 1))  
    processRight = 0;

//print Meta Info
  if (my_rank == 0){
    gmp_printf("KEY = %Zd \nsegmentSize: %Zd\nsegmentCount: %Zd\n", n, segmentSize, segmentCount);
  }

  MPI_Irecv(&buf, 1, MPI_INT, processLeft, 1, MPI_COMM_WORLD, &reqs);
  MPI_Test(&reqs, &flag, &stats);

  start_time = MPI_Wtime();
  while (count < p && !flag) {
    mpz_t start, end, mod, quot;
    mpz_init(start); mpz_init(end); mpz_init(mod); mpz_init(quot);

    mpz_set_ui(start, count * p + my_rank);
    mpz_mul(start, segmentSize, start);
    mpz_add(end, start, segmentSize);

    MPI_Test(&reqs, &flag, &stats);
// check if received
    if (flag) {
      MPI_Send(&my_rank, 1, MPI_INT, processRight, 1, MPI_COMM_WORLD);
      break;
    }

    while (mpz_cmp(start, end) <0 && !flag) {

      MPI_Test(&reqs, &flag, &stats);
      if (flag) {
        MPI_Send(&my_rank, 1, MPI_INT, processRight, 1, MPI_COMM_WORLD);
        break;
      }

      mpz_nextprime(start, start);
      mpz_mod(mod, n, start);

      if (mpz_sgn(mod) == 0) {

        mpz_cdiv_q(quot, n, start);
        if (mpz_probab_prime_p(quot, 10) >= 1){
          flag = 1;
          printf("process: %d, found the answer!!\n", my_rank);
          gmp_printf("****\n%Zd * %Zd = %Zd\n****\n", quot, start, n);
          MPI_Send(&my_rank, 1, MPI_INT, processRight, 1, MPI_COMM_WORLD);
        }
      }
    }
    count++;
  }

  while (!flag && p > 1) {
    MPI_Test(&reqs, &flag, &stats);
    // flag check
    if (flag) {
      MPI_Isend(&my_rank, 1, MPI_INT, processRight, 1, MPI_COMM_WORLD, &reqs);
      break;
    }
  }
  
  end_time = MPI_Wtime();
  run_time = end_time - start_time;

  if(my_rank != 0){
    MPI_Send(&run_time, 1, MPI_DOUBLE, 0, 3, MPI_COMM_WORLD);
  }
  else {
    FILE *fp;
    char filename[256];
    sprintf(filename,"time_%s",argv[1]);
    fp = fopen(filename,"w");

    fprintf(fp, "%d : %f\n", my_rank , run_time);
    double recv_time;
    for(i = 1 ; i < p; i++){
      MPI_Recv(&recv_time, 1, MPI_DOUBLE, i, 3, MPI_COMM_WORLD, &stats);
      fprintf(fp, "%d : %f\n", i , recv_time);
    }
    fclose(fp);
  }

  MPI_Finalize();
}
