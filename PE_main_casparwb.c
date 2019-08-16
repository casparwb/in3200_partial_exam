# include <stdlib.h>
# include <stdio.h>
# include <omp.h>
# include "PE_functions_casparwb.h"
# include <math.h>



int main(int narg, char **argv){
  FILE *infile;
  char buffer[100];

  double t, epsilon, d, *temp_ptr;
  int z, *z_ptr, n;
  double W_sum;     // variable for storing the sum of PageRank-values for the dangling weppages
  int nodes, edges; // Nodes: total number of webpages (N),
                    // Edges: total number of links (lines in file)

  d = atof(argv[2]);        // damping factor
  epsilon = atof(argv[3]);  // stopping threshold
  n = atoi(argv[4]);        // top n wepbages to print

  z = 0;                    // counter for number of dangling webpages
  z_ptr = &z;

  /*Arrays for storing the data in CRS-format*/
  double *val;
  int *col_idx, *row_ptr, *dangling_idx;

  double *top_pages = (double*)malloc(n*sizeof(double)); // array for storing the top pagerank scores
  int *top_pages_idx = (int*)malloc(n*sizeof(int));      // array for storing the indices of the top-scoring webpages

  /*Reading the input file and extracting the total number webpages and links*/
  infile = fopen(argv[1], "r");

  if (infile == NULL){
    printf("Unable to open file \n");
    exit(0);
  }

  fgets(buffer, sizeof(buffer), infile);               // Skip the first line
  fgets(buffer, sizeof(buffer), infile);               // Skip the second line

  fscanf(infile, "%*s %*s %d %*s %d", &nodes, &edges); // Collect the nodes and edges

  fgets(buffer, sizeof(buffer), infile);               // Skip the third line
  fgets(buffer, sizeof(buffer), infile);               // Skip the fourth line


  double *x_old = (double*)malloc(nodes*sizeof(double));       // array for x_(k-1)
  double *x_new = (double*)malloc(nodes*sizeof(double));       // array for x_k
  double *spmm_result = (double*)malloc(nodes*sizeof(double)); // array for storing the result of the sparse matrix-vector multiplication

  double one_over_N = 1.0/nodes;

  t = omp_get_wtime();
  read_graph_file(infile, &val, &col_idx, &row_ptr,
                  &dangling_idx, z_ptr, nodes, edges); // Initilize arrays with data from input file
  t = omp_get_wtime() - t;

  printf("\nReadfile time: %lfs\n", t);


  /*---Initialize x_old and x_new---*/
  for (int i = 0; i < nodes; i++){
    x_old[i] = one_over_N;
    x_new[i] = 1;
  }

  double diff, new_diff, max_diff = 1;
  int counter = 0;

  t = omp_get_wtime();
  while (max_diff > epsilon) {
    PageRank_iterations(val, col_idx, row_ptr, dangling_idx,
                        nodes, d, W_sum, x_old, x_new,
                        spmm_result, one_over_N, z);

    /*Finding the largest difference between x_(k,i) and x_((k-1), i)*/
    diff = 0;
    for (int i = 0; i < nodes; i++){
      new_diff = fabs(x_new[i] - x_old[i]);

      if (new_diff > diff) {
        diff = new_diff;
      } // end of if-statement
    } // end of for-loop

    max_diff = diff;

    /*Swapping pointers*/
    temp_ptr = x_new;
    x_new = x_old;
    x_old = temp_ptr;

    counter ++;

  } // end of while-loop

  printf("\nNumber of PageRank iterations: %d\n", counter);

  t = omp_get_wtime() - t;
  printf("\nPageRank calculation time: %lfs\n", t);


  t = omp_get_wtime();
  top_n_pages(x_old, top_pages, top_pages_idx, n, nodes); // Find the top n pagerank values
  t = omp_get_wtime() - t;

  printf("\nTop pages time: %lfs\n", t);

  printf("\nPosition: | PageRank Score: | Page index:\n ------------------------------\n");

  for (int i = 0; i < n; i++) printf("        %d | %lf        | %d\n", i+1, top_pages[i], top_pages_idx[i]);

  free(val); free(col_idx); free(row_ptr);
  free(x_old);
  free(x_new);
  free(spmm_result);
  free(top_pages); free(top_pages_idx);
  z != 0 ? free(dangling_idx) : 1;

  return 0;
}

/*
$ ./a.out web-NotreDame.txt 0.85 1e-10 100

Readfile time: 1461.799546s //= 24.45min
Number of PageRank iterations: 89
PageRank calculation time: 0.364812s
Top pages time: 0.077172s

 PageRank Score: | Page index:
 ------------------------------
 0.005627        | 1963
 0.005404        | 0
 0.003326        | 124802
 0.002857        | 12129
 0.002749        | 191267
 0.002732        | 32830
 0.002590        | 3451
 0.002460        | 83606
 0.002376        | 1973
 0.002341        | 142732
 0.002234        | 24823
 0.002152        | 143218
 0.001827        | 31331
 0.001802        | 149039
 0.001477        | 140170
 0.001403        | 12838
 0.001391        | 81878
 0.001169        | 226950
 0.001110        | 73859
 0.001101        | 292009
 0.000992        | 63364
 0.000990        | 24944
 0.000962        | 88448
 0.000947        | 88118
 0.000878        | 10331
 0.000841        | 143082
 0.000839        | 235904
 0.000835        | 193592
 0.000833        | 32833
 0.000830        | 198328
 0.000826        | 46468
 0.000818        | 14
 0.000793        | 32848
 0.000790        | 10335
 0.000768        | 10336
 0.000766        | 32846
 0.000745        | 87925
 0.000737        | 3478
 0.000685        | 73873
 0.000684        | 10332
 0.000680        | 307392
 0.000672        | 199031
 0.000669        | 151241
 0.000666        | 71616
 0.000665        | 155590
 0.000650        | 130721
 0.000617        | 212843
 0.000600        | 5123
 0.000598        | 15586
 0.000594        | 27684
 0.000581        | 238
 0.000550        | 32826
 0.000542        | 73875
 0.000537        | 307409
 0.000502        | 1
 0.000499        | 267537
 0.000499        | 243449
 0.000499        | 65194
 0.000499        | 217431
 0.000498        | 124559
 0.000498        | 245889
 0.000498        | 213868
 0.000498        | 215326
 0.000498        | 175927
 0.000496        | 73874
 0.000496        | 191257
 0.000495        | 110815
 0.000495        | 102186
 0.000495        | 285409
 0.000492        | 307408
 0.000490        | 8346
 0.000481        | 212812
 0.000480        | 285918
 0.000480        | 264860
 0.000480        | 248543
 0.000480        | 124619
 0.000480        | 177376
 0.000480        | 267647
 0.000480        | 243559
 0.000480        | 246014
 0.000480        | 213994
 0.000476        | 166681
 0.000476        | 159944
 0.000476        | 292967
 0.000469        | 16
 0.000461        | 191263
 0.000459        | 9860
 0.000431        | 29701
 0.000429        | 1828
 0.000415        | 9262
 0.000414        | 63706
 0.000410        | 146957
 0.000408        | 65233
 0.000407        | 124604
 0.000407        | 245938
 0.000407        | 213917
 0.000407        | 215375
 0.000407        | 175976
 0.000407        | 177361
 0.000405        | 124788*/
