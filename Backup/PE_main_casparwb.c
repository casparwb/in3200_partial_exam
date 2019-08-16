# include <stdlib.h>
# include <stdio.h>
# include <omp.h>
# include "PE_functions_casparwb.h"
# include <time.h>
# include <math.h>

clock_t clock(void);

int main(int narg, char **argv){
  time_t start, stop;
  double *val, t, epsilon = atof(argv[3]), d = atof(argv[2]);
  int *col_idx, *row_ptr, *dangling_idx, z, n = 3;
  int nodes, edges; // Nodes: total number of webpages (N),
                    // Edges: total number of links (lines in file)
  FILE *infile;

  double W_sum;     // variable for storing the sum of PageRank-values for the dangling wepages
  double *temp_ptr; // temporary pointer


  /*Reading the input file and extracting the total number webpages and links*/
  infile = fopen(argv[1], "r");
  char buffer[100];               // Array for storing lines that we have to skip

  if (infile == NULL){
    printf("Unable to open file \n");
    exit(0);
  }

  fgets(buffer, sizeof(buffer), infile);               // Skip the first line
  fgets(buffer, sizeof(buffer), infile);               // Skip the second line
  fscanf(infile, "%*s %*s %d %*s %d", &nodes, &edges); // Collect the nodes and edges
  fgets(buffer, sizeof(buffer), infile);               // Skip the third line
  fgets(buffer, sizeof(buffer), infile);               // Skip the fourth line

  /*In order to use the CRS arrays in other functions, we need to declare them in the main-function.
  Since we yet don't know what size they will have to be, we temporary allocate enough memory for
  their maximum possible size, and then reallocate the proper amount of memory later.*/

  // val          = (double*)malloc(edges*sizeof(double));        // array for storing non-zero values
  // col_idx      = (int*)   malloc(edges*sizeof(int));           // array for storing the column index of each nonzero value in val
  // row_ptr      = (int*)   malloc((nodes+1)*sizeof(int));       // array for storing the indices at which new rows start in val
  // dangling_idx = (int*)   malloc(nodes*sizeof(int));           // array for storing dangling webpage indices, if any

  double *x_old = (double*)malloc(nodes*sizeof(double));       // array for x_(k-1)
  double *x_new = (double*)malloc(nodes*sizeof(double));       // array for x_k
  double *spmm_result = (double*)malloc(nodes*sizeof(double)); // array for storing the result of the sparse matrix-vector multiplication

  double one_over_N = 1.0/nodes;


  t = omp_get_wtime();
  read_graph_file(infile, val, col_idx, row_ptr, dangling_idx, nodes, edges, z); // Initilize arrays with data from input file
  t = omp_get_wtime() - t;

  printf("Readfile time: %lf\n", t);
  //if (z == 0) free(dangling_idx);
  //
  // /*---Initialize x_old and x_new---*/
  // for (int i = 0; i < nodes; i++){
  //   x_old[i] = one_over_N;
  //   x_new[i] = 1;
  // }
  //
  // double diff, new_diff, max_diff = 1;
  // int counter = 0;
  //
  // t = omp_get_wtime();
  // while (max_diff > epsilon) {
  //   PageRank_iterations(val, col_idx, row_ptr, dangling_idx,
  //                       nodes, d, W_sum, x_old, x_new,
  //                       spmm_result, one_over_N, z);
  //
  //   /*Finding the largest difference between x_(k,i) and x_((k-1), i)*/
  //   diff = 0;
  //   for (int i = 0; i < nodes; i++){
  //     new_diff = fabs(x_new[i] - x_old[i]); //fabs for floats
  //     if (new_diff > diff) {
  //       diff = new_diff;
  //     } // end of if-statement
  //   } // end of for-loop
  //
  //   max_diff = diff;
  //
  //   /*Swapping pointers*/
  //   temp_ptr = x_new;
  //   x_new = x_old;
  //   x_old = temp_ptr;
  //
  //   counter ++;
  //
  //   // if (counter == 100) break;
  // } // end of while-loop
  // printf("No. of iterations: %d\n", counter);
  // t = omp_get_wtime() - t;
  // printf("PageRank time: %lf\n", t);
  //
  // double *top_pages = (double*)malloc(n*sizeof(double));
  // int *top_pages_idx = (int*)malloc(n*sizeof(int));

  // top_n_pages(x_old, top_pages, top_pages_idx, n, nodes); // Find the top n pagerank values
  //
  // for (int i = 0; i < n; i++) printf("%lf | %d\n", top_pages[i], top_pages_idx[i]);

  free(val); free(col_idx); free(row_ptr);
  free(x_old);
  free(x_new);
  free(spmm_result);

  return 0;
}

/*
Time without parallel: 2411s ~ 40 min
Time with parallel: 592.03s ~ 10 min

Speedup of almost 4 = number of cores
*/
