#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
This script generates networkx graph from SIF file.
Networkx graphs are required for shortest path analysis portion of this tool.
"""
# Author: Ahmet Sinan Yavuz <asinanyavuz@sabanciuniv.edu>
# BSD-3 License
import os
import pickle
import argparse
import networkx as nx
from time import gmtime, strftime

parser = argparse.ArgumentParser(description='SIF to networkx graph conversion tool')
parser.add_argument('-bv', '--bversion', help='Version of SIF file',
                    default=strftime("%a, %d %b %Y %H:%M:%S +0000", gmtime()))
parser.add_argument('-o', '--output', help='Output directory', default="../data")
parser.add_argument('-v', '--verbose', help='Set verbosity level (report per N edges, default: 100000)', default=100000)

requiredNamed = parser.add_argument_group('required named arguments')
requiredNamed.add_argument('-i', '--input', help='SIF file', required=True)
args = parser.parse_args()

links_file = args.input
version = args.bversion
ppi = nx.Graph()
added = 0

with open(links_file) as ifh:
    for line in ifh:
        elems = line.split()
        gene1 = elems[0]
        gene_other = elems[2:]
        for gene2 in gene_other:
            ppi.add_edge(gene1.upper(), gene2.upper(), interaction_type=elems[1])
            added += 1
        if added % int(args.verbose) == 0:
            print("Added %d edges." % added)

print("Finished.")
print("Added %d edges." % added)
print("Total nodes: %d, total edges: %d." %(len(ppi.nodes()), len(ppi.edges())))
pickle.dump((version, ppi), open(args.output+"/SIF_PPI.dat", "wb"), pickle.HIGHEST_PROTOCOL)
print("Saved to \"%s\". Please use this path in \"-ppi\" parameter of analysis scripts." % os.path.abspath(args.output+"/SIF_PPI.dat"))
