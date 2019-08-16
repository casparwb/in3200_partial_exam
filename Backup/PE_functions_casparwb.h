
void read_graph_file(FILE *infile, double *val, int *col_idx,
              int *row_ptr, int *dangling_idx, int nodes, int edges, int z) {

  int temp1, temp2, temp_counter;
  int i, j, k;
  int nz = 0; // no. of nonzero elements
  char buffer[100];  // Array for storing lines that we have to skip

  int *FromNode      = (int*)malloc(edges*sizeof(int)); // array for storing the FromNodeId-values
  int *ToNode        = (int*)malloc(edges*sizeof(int)); // array for storing the ToNodeId-values
  int *num_out_links = (int*)malloc(nodes*sizeof(int)); // array for storing the no. of outbound links from each webpage

  z = 0;
  /**********Collecting the FromeNodeId and ToNodeId data**********/
  i = 0;
  while (fgets(buffer, sizeof(buffer), infile) != NULL){
    sscanf(buffer, "%d %d", &temp1, &temp2);
    FromNode[i] = temp1;
    ToNode[i] = temp2;
    i++;
  } // end of while-loop
  fclose(infile);


  /**********Collecting total number of FromNode links for each website, and total number of nonzero values**********/
  #pragma omp parallel for private(i) reduction(+:nz)
  for (i = 0; i < nodes; i++){
    num_out_links[i] = 0; // Initialize;

    for (j = 0; j < edges; j++){
      if ((FromNode[j] == i) && (ToNode[j] != FromNode[j])) {   // if website i is in FromNode[j], and it is not a self-link
      num_out_links[i] += 1;
      nz++;
      } // end of if-statement
    } // end of inner loop
    if (num_out_links[i] == 0) z++;
  } // end of outer loop

  if (z != 0){  // If there are any dangling webpages, collect their indices
    dangling_idx = (int*) realloc(dangling_idx, (z+1)*sizeof(int));
    temp_counter = 0;

    for (i = 0; i < nodes; i++){

      if (num_out_links[i] == 0){
        dangling_idx[temp_counter] = i;
        temp_counter++;

      } // end of inner if-statement
    }   // end of inner loop
  }     // end of outer if-statement


  /**********Building the hyperlink matrix in CRS-format**********/
  val = (double*)realloc(val, nz*sizeof(double));  // Reallocate the size of the array if necessary
  col_idx = (int*)realloc(col_idx, nz*sizeof(int));

  // temp_counter = 0;
  // for (i = 0; i < nodes; i++){
  //   row_ptr[i] = temp_counter;   // Set the row_ptr index
  //
  //   for (j = 0; j < edges; j++){
  //
  //     if ((ToNode[j] == i) && (ToNode[j] != FromNode[j])){
  //       val[temp_counter] = 1.0/num_out_links[FromNode[j]]; // locate the nonzero values and put them in val
  //       col_idx[temp_counter] = FromNode[j];                // set the appropriate column index for each value
  //
  //       temp_counter++;
  //
  //     } // end of if-statement
  //
  //   } // end of inner loop
  // } // end of outer loop


  row_ptr[nodes] = nz;  // last row_ptr index should be total number of nonzero elements

  free(FromNode);
  free(ToNode);
  free(num_out_links);

} // end of read_graph_file function


void PageRank_iterations(double *val, int *col_idx,
                         int *row_ptr, int *dangling_idx,
                         int nodes, double d, double W_sum,
                         double *x_old, double *x_new,
                         double *spmm_result, double one_over_N, int z) {




    W_sum = 0;
    if (z = 0) dangling_idx[0] = 0; // to stop segmentation fault
    else{
      #pragma omp parallel for reduction(+:W_sum)
      for (int i = 0; i < z; i++) {
        W_sum += x_old[dangling_idx[i]];
      }
    }

    for (int i = 0; i < nodes; i++){
      spmm_result[i] = 0; // Initilize

      for (int j = row_ptr[i]; j < row_ptr[i+1]; j++){
          spmm_result[i] += val[j] * x_old[col_idx[j]];
      } //end of inner loop

      spmm_result[i] *= d;
    }

    for (int i = 0; i < nodes; i++){
      x_new[i] = ((1 - d + d*W_sum)*one_over_N) + spmm_result[i];
    }
} // end of PageRank-function

void top_n_pages(double *x, double *top_pages, int *top_pages_idx, int n, int nodes){
  /*Finding the top n PageRank scores and their corresponding website indices*/
  double temp_max, max_idx, max = 1;

  for (int j = 0; j < n; j++){
    temp_max = 0;

    for (int i = 0; i < nodes; i++){
      if ((x[i] > temp_max) && (x[i] < max)){
        temp_max = x[i];
        max_idx = i;
      } // end of if-statement
    } // end of inner for loop

    max = temp_max;
    top_pages_idx[j] = max_idx;
    top_pages[j] = max;

  } // end of outer for loop

} // end of function
