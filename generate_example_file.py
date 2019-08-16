import numpy as np
from random import randrange

def newfile(filename, nodes):
    edges = 2*nodes
    outfile = open(filename, "w+")

    outfile.write("# Sometext1\n")
    outfile.write("# Sometext2\n")
    outfile.write("# Nodes: {} Edges: {}\n" .format(nodes, edges))
    outfile.write("# FromNodeId ToNodeId\n")


    for i in range(edges):
        outfile.write("{} {}\n" .format(randrange(0, nodes-1), randrange(0, nodes-1)))


newfile("generated_file_example.txt", 1000)
