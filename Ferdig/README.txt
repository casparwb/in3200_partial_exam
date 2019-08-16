The program is compiled using the GCC compiler, using the OpenMp command-line option -fopenmp. 

Running the produced output file requires the input of four arguments in the following sequence:
-------------------------------------------------------
1. A text file structured in the appropriate manner.
2. The damping factor d.
3. The convergence threshold value epsilon.
4. The n value to display the top n webpages at end of runtime.
-------------------------------------------------------

An example of how to compile and run the program:
-------------------------------------------------------
$ gcc -fopenmp PE_main_15243.c

$ ./a.out 100nodes_graph.txt 0.85 1e-5 10
--------------------------------------------------------

Which produces the output:
--------------------------------------------------------
Provided arguments:
Damping factor = 0.850000, epsilon = 1.000000e-5, n = 10

Readfile time: 0.000495s

Total number of PageRank iterations: 15

PageRank calculation time: 0.000444s = 0.443700ms

Top pages time: 0.000005s = 0.005000ms

Position: | PageRank Score: | Page index:
 ------------------------------
        1 | 0.049372        | 64
        2 | 0.038200        | 67
        3 | 0.035346        | 26
        4 | 0.030349        | 52
        5 | 0.024108        | 8
        6 | 0.023089        | 46
        7 | 0.021599        | 60
        8 | 0.021360        | 62
        9 | 0.018544        | 5
        10 | 0.018177        | 68
---------------------------------------------------------