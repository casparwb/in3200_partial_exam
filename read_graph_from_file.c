#include <stdlib.h>
#include <stdio.h>
#include <omp.h>
#include <time.h>

clock_t clock(void);

void readfile(FILE *infile, double *val, int *col_idx,
              int *row_ptr, int *dangling_idx, int nodes, int edges) {

  int temp1, temp2;
  int i, j, k;

  char buffer[100]; // Array for storing lines that we have to skip

  int *FromNode = (int*)malloc(edges*sizeof(int)); //array for storing the FromNodeId-values
  int *ToNode = (int*)malloc(edges*sizeof(int)); // array for storing the ToNodeId-values
  int *num_out_links = (int*)malloc(nodes*sizeof(int)); // no. of outbound links from each webpage

  // Collecting the FromeNodeId and ToNodeId data
  i = 0;
  while (fgets(buffer, sizeof(buffer), infile) != NULL){
    sscanf(buffer, "%d %d", &temp1, &temp2);
    FromNode[i] = temp1;
    ToNode[i] = temp2;
    i++;
  } // end of while-loop
  fclose(infile);


  //Collecting total number of FromNode links for each website, and total number of nonzero values
  int nz = 0, z = 0; // no. of nonzero elements and number of zero elements
  /*The z-counter is for sizing the danging-webpages index array later*/
  int temp_counter; // Counter k


  #pragma omp parallel for private(i,j) reduction(+:nz)
  for (i = 0; i < nodes; i++){
    num_out_links[i] = 0; // Initialize L;

    for (j = 0; j < edges; j++){
      if ((FromNode[j] == i) && (ToNode[j] != FromNode[j])) {
      num_out_links[i] += 1;
      nz++;
      } // end of if-statement
    } // end of inner loop
  if (num_out_links[i] == 0) z++;
  } // end of outer loop

  if (z != 0){ // If there are any dangling webpages
    dangling_idx = (int*) realloc(dangling_idx, z*sizeof(int));
    temp_counter = 0;
    for (i = 0; i < nodes; i++){
      if (num_out_links[i] == 0){
        dangling_idx[temp_counter] = i;
        temp_counter++;
      } // end of inner if-statement
    }   // end of inner loop
  }     // end of outer if-statement
  else free(dangling_idx);


  /* Building the hyperlink matrix in CRS-format */
  val     = (double*) realloc(val, nz*sizeof(double));
  col_idx = (int*)    realloc(col_idx, nz*sizeof(int));
  row_ptr = (int*)    realloc(row_ptr, (nodes+1)*sizeof(int));


  temp_counter = 0;
  #pragma omp parallel for private(i, j) reduction(+:temp_counter)
  for (i = 0; i < nodes; i++){
    row_ptr[i] = temp_counter;
    //printf("%d\n", row_ptr[i]);
    for (j = 0; j < edges; j++){

      if ((ToNode[j] == i) && (ToNode[j] != FromNode[j])){
        val[temp_counter] = 1.0/num_out_links[FromNode[j]];
        col_idx[temp_counter] = FromNode[j];
        //printf("%lf\n", val[temp_counter]);
        //printf("%d\n", col_idx[temp_counter]);
        temp_counter++;
      } // end of if-statement
    } // end of inner loop
  } // end of outer loop
  row_ptr[nodes] = nz;
  free(FromNode);
  free(ToNode);
  free(num_out_links);

}

double dangling_PageRank(double *x_old, int *dangling_idx){
  /*Function for calculating the sum of the pagerank-score of dangling-webpages*/
  int n = sizeof(dangling_idx)/sizeof(int);
  // printf("%d\n", n);
  int sum = 0;

  for (int i = 0; i < n; i++) {
    sum += x_old[dangling_idx[i]];
    // printf("%lf\n", x_old[dangling_idx[i]]);
  }
  return sum;
}

void spmm(int nodes, double *val, int *col_idx, int *row_ptr,
          double *x_old, double *spmm_result, double d){
  /*Performing the Sparse Matrix vector multiplictation A*x^(k-1)
  spmm_result is the vector for storing the result of the multiplictation*/

  int i, j;
  for (i = 0; i < nodes; i++){
    spmm_result[i] = 0; // Initilize

    for (j = row_ptr[i]; j <= row_ptr[i+1] - 1; j++){
        spmm_result[i] += val[j] * x_old[col_idx[j]];

    } //end of inner loop

    spmm_result[i] *= d;

  }   // end of outer loop
}     // end of function

void PageRank_iterations(double *val, int *col_idx,
                         int *row_ptr, int *dangling_idx,
                         int nodes, double d, double W_sum,
                         double *x_old, double *x_new,
                         double *spmm_result, double one_over_N) {


  W_sum = dangling_PageRank(x_old, dangling_idx);
  spmm(nodes, val, col_idx, row_ptr, x_old, spmm_result, d);
  //W_sum og spmm beregnes bare én gang per iterasjon av PageRank_iterations

  for (int i = 0; i < nodes; i++){
    x_new[i] = ((1 - d + d*W_sum)*one_over_N) + spmm_result[i];
  }

  // for (int i = 0; i < nodes; i++) printf("%d  |  %lf\n", i, x_old[i]);
  // printf("\n New iteration \n");
}


int main(int narg, char **argv){
  time_t start, stop;
  double *val, t;
  int *col_idx, *row_ptr, *dangling_idx;
  int nodes, edges; // Nodes: total number of webpages (N),
                   // Edges: total number of links (lines in file)

  double epsilon, d = atof(argv[2]);
  int n; // top n webpages
  /*Reading the input file and extracting the total number webpages and links*/
  FILE *infile;
  infile = fopen(argv[1], "r");
  char buffer[100];               // Array for storing lines that we have to skip

  if (infile == NULL){
    printf("Unable to open file \n");
    exit(0);
  }

  fgets(buffer, sizeof(buffer), infile); // Skip the first line
  fgets(buffer, sizeof(buffer), infile); // Skip the second line
  fscanf(infile, "%*s %*s %d %*s %d", &nodes, &edges); // Collect the nodes and edges
  fgets(buffer, sizeof(buffer), infile); // Skip the third line
  fgets(buffer, sizeof(buffer), infile); // Skip the fourth line

  /*In order to use the CRS arrays in other functions, we need to declare them in the main-function.
  Since we yet don't know what size they will have to be, we temporary allocate enough memory for
  their maximum possible size, and then reallocate the proper amount of memory later.*/

  val          = (double*)malloc(edges*sizeof(double));    // array for storing non-zero values
  col_idx      = (int*)   malloc(edges*sizeof(int));      // array for storing the column index of each nonzero value in val
  row_ptr      = (int*)   malloc(nodes*sizeof(int));      // array for storing the indices at which new rows start in val
  dangling_idx = (int*)   malloc(nodes*sizeof(int)); // array for storing dangling webpage indices, if any

  /*Initilize arrays with data from input file*/

  readfile(infile, val, col_idx, row_ptr, dangling_idx, nodes, edges);

  double W_sum;
  double *temp_ptr;

  double *x_old = (double*)malloc(nodes*sizeof(double));
  double *x_new = (double*)malloc(nodes*sizeof(double));
  double *spmm_result = (double*)malloc(nodes*sizeof(double));

  double one_over_N = 1.0/nodes;


  for (int i = 0; i < nodes; i++){
    x_old[i] = one_over_N;
    x_new[i] = 1;
  }


  PageRank_iterations(val, col_idx, row_ptr, dangling_idx,
                      nodes, d, W_sum, x_old, x_new,
                      spmm_result, one_over_N);
  for (int i = 0; i < nodes; i++){
    printf("%lf | %lf\n", x_old[i], x_new[i]);
  }

  temp_ptr = x_old;
  x_old = x_new;
  x_new = temp_ptr;
  printf("NEW ITERTION\n");
  PageRank_iterations(val, col_idx, row_ptr, dangling_idx,
                      nodes, d, W_sum, x_old, x_new,
                      spmm_result, one_over_N);
  for (int i = 0; i < nodes; i++){
    printf("%lf | %lf\n", x_old[i], x_new[i]);
    }
  return 0;
}

/*Problemer med dangling_PageRank? dangling_idx har lengde 2 når den skal ha 1
*/


/*Expected result from example_file_small with d = 0.85:

x_0 = [1/8, 1/8, 1/8, 1/8, 1/8, 1/8, 1/8, 1/8]
spmm_result(x_0) = d*[0.0417, 0.1667, 0.0625, 0.1250, 0.1458, 0.1458, 0.1042, 0.2083]

x_1 = [0.0542, 0.1604, 0.0719, 0.1250, 0.1427, 0.1427, 0.1073, 0.1958]
x_2 = [0.0491, 0.1077, 0.0418, 0.1551, 0.1151, 0.1778, 0.1424, 0.2109]

*/

// Kan vi kollapse noen loops?
// Gjøre nodes and edges global?
