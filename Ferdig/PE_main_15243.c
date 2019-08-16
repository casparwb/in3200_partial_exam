# include <stdlib.h>
# include <stdio.h>
# include <omp.h>
# include "PE_functions_15243.c"




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

  printf("\nProvided arguments:\nDamping factor = %lf, epsilon = %lfe-5, n = %d\n", d, epsilon*100000, n);

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


  t = omp_get_wtime();

  read_graph_file(infile, &val, &col_idx, &row_ptr,    // Initilize arrays with data from input file
                  &dangling_idx, z_ptr, nodes, edges);

  t = omp_get_wtime() - t;
  printf("\nReadfile time: %lfs\n", t);



  t = omp_get_wtime();

  PageRank_iterations(val, col_idx, x_old, x_new, spmm_result, row_ptr, dangling_idx, // Calculate the pagerank values
                      nodes, d, z, epsilon);

  t = omp_get_wtime() - t;
  printf("\nPageRank calculation time: %lfs = %lfms\n", t, t*1000);


  t = omp_get_wtime();

  top_n_pages(x_old, top_pages, top_pages_idx, n, nodes); // Find the top n pagerank values

  t = omp_get_wtime() - t;


  printf("\nTop pages time: %lfs = %lfms\n", t, t*1000);

  printf("\nPosition: | PageRank Score: | Page index:\n ------------------------------\n");

  for (int i = 0; i < n; i++)
    printf("        %d | %lf        | %d\n", i+1, top_pages[i], top_pages_idx[i]);

  free(val); free(col_idx); free(row_ptr);
  free(x_old);
  free(x_new);
  free(spmm_result);
  free(top_pages); free(top_pages_idx);
  z != 0 ? free(dangling_idx) : 1;

  return 0;
}

/*
./a.out 100nodes_graph.txt 0.85 1e-10 50

Readfile time: 0.000390s

PageRank calculation time: 0.000836s

Top pages time: 0.000026s

Position: | PageRank Score: | Page index:
 ------------------------------
        1 | 0.045687        | 64
        2 | 0.035357        | 67
        3 | 0.032701        | 26
        4 | 0.028093        | 52
        5 | 0.022310        | 8
        6 | 0.021358        | 46
        7 | 0.020001        | 60
        8 | 0.019767        | 62
        9 | 0.017172        | 5
        10 | 0.016830        | 68
        11 | 0.016378        | 7
        12 | 0.015489        | 94
        13 | 0.015413        | 3
        14 | 0.015362        | 24
        15 | 0.014422        | 92
        16 | 0.014257        | 13
        17 | 0.014234        | 93
        18 | 0.013971        | 30
        19 | 0.013849        | 81
        20 | 0.013423        | 34
        21 | 0.013312        | 11
        22 | 0.013208        | 61
        23 | 0.013104        | 56
        24 | 0.012987        | 18
        25 | 0.012911        | 14
        26 | 0.012861        | 29
        27 | 0.012782        | 1
        28 | 0.012704        | 84
        29 | 0.012432        | 19
        30 | 0.011740        | 71
        31 | 0.011729        | 74
        32 | 0.011439        | 70
        33 | 0.011244        | 80
        34 | 0.010693        | 25
        35 | 0.010667        | 82
        36 | 0.010634        | 27
        37 | 0.010530        | 40
        38 | 0.010239        | 28
        39 | 0.010143        | 57
        40 | 0.009807        | 97
        41 | 0.009471        | 63
        42 | 0.008844        | 22
        43 | 0.008775        | 4
        44 | 0.008497        | 59
        45 | 0.008341        | 33
        46 | 0.008333        | 2
        47 | 0.008298        | 38
        48 | 0.008159        | 23
        49 | 0.007697        | 45
        50 | 0.007626        | 32*/
