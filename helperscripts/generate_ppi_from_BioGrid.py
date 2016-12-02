#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
This script generates networkx graph from BioGrid links file.
Networkx graphs are required for shortest path analysis portion of this tool.
"""
# Author: Ahmet Sinan Yavuz <asinanyavuz@sabanciuniv.edu>
# BSD-3 License
import os
import pickle
import argparse
import networkx as nx
from time import gmtime, strftime

parser = argparse.ArgumentParser(description='BioGrid to networkx graph conversion tool')
parser.add_argument('-bv', '--bversion', help='Version of BioGrid links file',
                    default=strftime("%a, %d %b %Y %H:%M:%S +0000", gmtime()))
parser.add_argument('-o', '--output', help='Output directory', default="../data")
parser.add_argument('-v', '--verbose', help='Set verbosity level (report per N edges, default: 100000)', default=100000)

requiredNamed = parser.add_argument_group('required named arguments')
requiredNamed.add_argument('-i', '--input', help='BioGrid links file', required=True)
args = parser.parse_args()

links_file = args.input
version = args.bversion
ppi = nx.Graph()
added = 0

with open(links_file) as ifh:
    next(ifh)  # Skip the header line
    for line in ifh:
        elems = line.split()
        gene1 = elems[7]
        gene2 = elems[8]
        ppi.add_edge(gene1.upper(), gene2.upper(), source=elems[14])
        added += 1
        if added % int(args.verbose) == 0:
            print("Added %d edges." % added)

print("Finished.")
print("Added %d edges." % added)
print("Total nodes: %d, total edges: %d." %(len(ppi.nodes()), len(ppi.edges())))
pickle.dump((version, ppi), open(args.output+"/BioGrid.dat", "wb"), pickle.HIGHEST_PROTOCOL)
print("Saved to \"%s\". Please use this path in \"-ppi\" parameter of analysis scripts." % os.path.abspath(args.output+"/BioGrid.dat"))
