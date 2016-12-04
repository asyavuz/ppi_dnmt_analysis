#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
PPI-based gene-list enrichment tool. (modified from tool described by Durasi & Sezerman [in preparation])
Ideally, please run this file with networkx installed PyPy2, as it reduces runtime noticeably.
However, as PyPy2 isn't compatible with scipy, hence, p-value calculation and enrichment is presented in another script.
"""
# Author: Ahmet Sinan Yavuz <asinanyavuz@sabanciuniv.edu>
# BSD-3 License

import os
import sys
import math
import pickle
import warnings
import argparse
import networkx as nx
import multiprocessing
from builtins import range
from time import gmtime, strftime

orig_cmd = " ".join(sys.argv)
NUM_PROCS = multiprocessing.cpu_count()

parser = argparse.ArgumentParser(description='PPI-based gene-list enrichment tool - shortest pathway generation')
parser.add_argument('-m', '--mirna', help='List of important miRNAs identified by the experiment', default="")
parser.add_argument('-lg', '--litgene', help='Known gene list associated with the condition',
                    default="./data/LiteratureGenes.txt")
parser.add_argument('-lm', '--litmirna', help='Known miRNA list associated with the condition',
                    default="./data/LiteraturemiRNAs.txt")
parser.add_argument('-p', '--ppi', help='PPI network file (only pickled networkx graphs)',
                    default="./data/BioGrid.dat")
parser.add_argument('-s', '--start', help="Start nodes (allowed: expt, known, both, or custom, default: expt)",
                    default="expt")
parser.add_argument('-e', '--end', help="End nodes (allowed: expt, known, both, or custom, default: known)",
                    default="known")
parser.add_argument('-cs', '--cstart', help="Custom start node file", default="")
parser.add_argument('-ce', '--cend', help="Custom end node file", default="")
parser.add_argument('-sc', '--score', help="Score nodes in list (allowed: expt, known, both, default: both)",
                    default="both")
parser.add_argument('-t', '--tarbase', help='TARBASE database file', default="./data/tarbase_data.csv")
parser.add_argument('-tf', '--tarbasefilter', help='TARBASE filter string (separated by ;, default: '
                                                   '\'species=Homo sapiens\')',
                    default="species=Homo sapiens")
parser.add_argument('-exp', '--expressed', help="List of expressed genes", default="")
parser.add_argument('-o', '--output', help='Output folder name (default: Output)', default='Output')
parser.add_argument('-nc', '--cores', help='Number of processes (default: %d)' % NUM_PROCS, default=NUM_PROCS)
parser.add_argument('-v', '--verbose', help='Set verbosity level (report per N paths, default: 10000)', default=10000)

requiredNamed = parser.add_argument_group('required named arguments')
requiredNamed.add_argument('-g', '--gene', help='List of important genes identified by the experiment', required=True)
args = parser.parse_args()

# Input files
ppi_network = args.ppi
tarbase_file = args.tarbase
known_genes_file = args.litgene
known_mirna_file = args.litmirna

identified_genes_file = args.gene
identified_mirna_file = args.mirna

expressed_gene_file = args.expressed

# Output folder
outfolder = args.output

# TARBASE filter, only equality filters are allowed.
tmp = args.tarbasefilter
tarbase_filter = tmp.split(";")

# Start and end nodes "expt" for experimental genes/miRNAs, "known" for known genes/miRNAs,
# "both" for both "expt" and "known", or "custom" for custom file. Custom file may contain mixture of genes and miRNAs.
start_node = args.start
end_node = args.end
custom_start = args.cstart
custom_end = args.cend

if start_node == "custom" and not custom_start:
    raise ValueError("Cannot set start node mode as custom without specifying the custom node list file with --cstart "
                     "parameter!")

if end_node == "custom" and not custom_end:
    raise ValueError("Cannot set end node mode as custom without specifying the custom node list file with --cend "
                     "parameter!")

report_per = int(args.verbose)

# ################################################ #
#       DO NOT TOUCH BELOW HERE UNLESS NEEDED      #
# ################################################ #

if args.cores:
    NUM_PROCS = int(args.cores)

if not os.path.exists(outfolder):
    os.makedirs(outfolder)
stats_fh = open(outfolder + "/Run_Information.txt", "a")
stats_fh.write("--------------------------------------------\n")
stats_fh.write("Run Time: %s\n" % strftime("%a, %d %b %Y %H:%M:%S +0000", gmtime()))
stats_fh.write("Run Command: %s\n" % orig_cmd)

print("Loading PPI network...")
version, ppi = pickle.load(open(ppi_network, "rb"))
ppigenes = ppi.nodes()
print("-Loaded.")
stats_fh.write("PPI File: %s (%d Nodes, %d Edges, Version: %s)\n" % (ppi_network, len(ppi.nodes()), len(ppi.edges()),
                                                                     version))
if (identified_mirna_file or known_mirna_file) and tarbase_file:
    mirna_targets = {}
    print("Loading TARBASE data (Filters: %s)..." % (" AND ".join(tarbase_filter)))
    with open(tarbase_file) as f:
        header = {}
        total = 0
        unfiltered_total = 0
        for i, line in enumerate(f):
            if i == 0:
                elems = line.split("\t")

                # Store header element indexes to use it in filtering
                for j, elem in enumerate(elems):
                    header[elem.strip()] = j
            else:
                elems = line.split("\t")
                elems[-1] = elems[-1].strip()
                unfiltered_total += 1

                filter_line = False
                if tarbase_filter:
                    for filt in tarbase_filter:
                        filt_elems = filt.split("=")
                        if len(filt_elems) < 2 or filt_elems[0] not in header.keys():
                            warnings.warn("Unrecognised filter: \"%s\". Passing..." % filt, RuntimeWarning)
                            continue

                        if elems[header[filt_elems[0]]] != filt_elems[1]:
                            filter_line = True
                if not filter_line:
                    if elems[2].upper() in mirna_targets:
                        mirna_targets[elems[2].upper()].append(elems[1].upper())
                    else:
                        mirna_targets[elems[2].upper()] = [elems[1].upper()]

                    ppi.add_edge(elems[2].upper(), elems[1].upper(), source="TARBASE")
                    total += 1

    print("-Loaded. Parsed %d miRNA-Target interactions and added these into PPI." % total)

    stats_fh.write("miRNA Target Database: %s (Unfiltered: %d, Filtered: %d Interactions)\n" % (tarbase_file,
                                                                                                unfiltered_total,
                                                                                                total))
    if tarbase_filter:
        stats_fh.write("miRNA Target Database Filter: %s\n" % tarbase_filter)
    stats_fh.write("PPI Size After Added miRNA Interactions: %d Nodes, %d Edges)\n" % (len(ppi.nodes()),
                                                                                       len(ppi.edges())))

known_genes = []
known_mirna = []
experimental_genes = []
experimental_mirna = []
ppigenes = ppi.nodes()
print("Loading known genes (File: %s)..." % known_genes_file)
unrecognised = 0
total = 0
with open(known_genes_file) as f:
    for i, line in enumerate(f):
        total += 1
        try:
            # If gene is in PPI line below will work without any problem, if not it will raise ValueError
            b = ppigenes.index(line.strip().upper())
            known_genes.append(line.strip().upper())
        except ValueError:
            unrecognised += 1

known_genes = list(set(known_genes))
print("-Loaded. Recognised %d unique genes as in PPI out of %d entries." % (len(known_genes), total))
stats_fh.write("Known Genes: %s (%d Entries, %d recognised as in PPI)\n" % (known_genes_file, total, len(known_genes)))

if known_mirna_file:
    print("Loading known miRNAs (File: %s)..." % known_mirna_file)
    unrecognised = 0
    total = 0
    with open(known_mirna_file) as f:
        for i, line in enumerate(f):
            total += 1
            try:
                # If gene is in PPI line below will work without any problem, if not it will raise ValueError
                b = list(mirna_targets.keys()).index(line.strip().upper())
                known_mirna.append(line.strip().upper())
            except ValueError:
                unrecognised += 1
    known_mirna = list(set(known_mirna))
    print("-Loaded. Recognised %d unique miRNAs as in PPI out of %d entries." % (len(known_mirna), total))
    stats_fh.write("Known miRNAs: %s (%d miRNAs, %d recognised as in PPI)\n" % (known_mirna_file, total,
                                                                                len(known_mirna)))

print("Loading experimental genes (File: %s)..." % identified_genes_file)
unrecognised = 0
total = 0
with open(identified_genes_file) as f:
    for i, line in enumerate(f):
        total += 1
        try:
            # If gene is in PPI line below will work without any problem, if not it will raise ValueError
            b = ppigenes.index(line.strip().upper())
            experimental_genes.append(line.strip().upper())
        except ValueError:
            unrecognised += 1
experimental_genes = list(set(experimental_genes))
print("-Loaded. Recognised %d unique genes as in PPI out of %d entries." % (len(experimental_genes), total))
stats_fh.write("Experimental Genes: %s (%d entries, %d recognised as in PPI)\n" % (identified_genes_file,
                                                                                   total,
                                                                                   len(experimental_genes)))

if identified_mirna_file:
    print("Loading experimental miRNAs (File: %s)..." % identified_mirna_file)
    unrecognised = 0
    total = 0
    with open(identified_mirna_file) as f:
        for i, line in enumerate(f):
            total += 1
            try:
                # If gene is in PPI line below will work without any problem, if not it will raise ValueError
                b = list(mirna_targets.keys()).index(line.strip().upper())
                experimental_mirna.append(line.strip().upper())
            except ValueError:
                unrecognised += 1

    experimental_mirna = list(set(experimental_mirna))
    print("-Loaded. Recognised %d unique miRNAs as in PPI out of %d entries." % (len(experimental_mirna), total))
    stats_fh.write("Experimental Genes: %s (%d entries, %d recognised as in PPI)\n" % (identified_mirna_file, total,
                                                                                       len(experimental_mirna)))

expressed_genes = []
ppigenes = ppi.nodes()
if expressed_gene_file:
    print("Loading expressed genes (File: %s)..." % expressed_gene_file)
    unrecognised = 0
    total = 0
    with open(expressed_gene_file) as f:
        for i, line in enumerate(f):
            total += 1
            try:
                # If gene is in PPI line below will work without any problem, if not it will raise ValueError
                b = ppigenes.index(line.strip().upper())
                expressed_genes.append(line.strip().upper())
            except ValueError:
                unrecognised += 1

    expressed_genes = list(set(expressed_genes))
    print("-Loaded. Recognised %d unique expressed genes as in PPI out of %d entries." % (len(expressed_genes),
                                                                                          total))
    stats_fh.write("Expressed Genes: %s (%d entries, %d recognised as in PPI)\n" % (expressed_gene_file, total,
                                                                                    len(experimental_genes)))

entrez_to_hgnc = {}
hgnc_to_entrez = {}
print("Loading gene ID conversions...")
with open("./data/Entrez_to_HGNC.txt") as f:
    for i, line in enumerate(f):
        if i == 0:
            continue
        elems = line.split()
        entrez_to_hgnc[elems[0]] = elems[1].strip()
        hgnc_to_entrez[elems[1].strip()] = elems[0]
print("Loaded mapping data for %d IDs." % len(entrez_to_hgnc))

custom_end_elems = []
custom_start_elems = []
bg_genes = []

if custom_start:
    print("Loading custom start node labels (File: %s)..." % custom_start)
    unrecognised = 0
    total = 0
    with open(custom_start) as f:
        for i, line in enumerate(f):
            total += 1
            try:
                b = ppigenes.index(line.strip().upper())
                custom_start_elems.append(line.strip().upper())
            except ValueError:
                unrecognised += 1

    custom_start_elems = list(set(custom_start_elems))
    print("-Loaded. Recognised %d unique genes as in PPI out of %d entries." % (len(custom_start_elems), total))
    stats_fh.write("Custom Start Nodes File: %s (%d entries,  %d recognised as in PPI)\n" % (custom_start, total,
                                                                                             len(custom_start_elems)))

if custom_end:
    print("Loading custom end node labels (File: %s)..." % custom_end)
    unrecognised = 0
    total = 0
    with open(custom_end) as f:
        for i, line in enumerate(f):
            total += 1
            try:
                b = ppigenes.index(line.strip().upper())
                custom_end_elems.append(line.strip().upper())
            except ValueError:
                unrecognised += 1

    custom_end_elems = list(set(custom_end_elems))
    print("-Loaded. Recognised %d unique genes as in PPI out of %d entries." % (len(custom_end_elems), total))
    stats_fh.write("Custom End Nodes File: %s (%d entries, %d recognised as in PPI)\n" % (custom_end, total,
                                                                                          len(custom_end_elems)))


start_nodes = []
end_nodes = []

# In order to get rid of potential repeated genes, convert to set then back to list
all_nodes = list(set(experimental_genes + experimental_mirna + known_genes + known_mirna))
all_mirnas = experimental_mirna + known_mirna
if start_node == "expt":
    start_nodes = list(set(experimental_genes + experimental_mirna))
elif start_node == "known":
    start_nodes = list(set(known_genes + known_mirna))
elif start_node == "both":
    start_nodes = all_nodes
elif start_node == "custom":
    start_nodes = custom_start_elems
else:
    raise ValueError("Unrecognised start node selection: %s!" % start_node)

if end_node == "expt":
    end_nodes = list(set(experimental_genes + experimental_mirna))
elif end_node == "known":
    end_nodes = list(set(known_genes + known_mirna))
elif end_node == "both":
    end_nodes = all_nodes
elif end_node == "custom":
    end_nodes = custom_end_elems
else:
    raise ValueError("Unrecognised end node selection: %s!" % end_node)

score = []
if args.score == "expt":
    score = list(set(experimental_genes + experimental_mirna))
elif args.score == "known":
    score = list(set(known_genes + known_mirna))
elif args.score == "both":
    score = all_nodes

if expressed_genes:
    orig_start = start_nodes
    orig_end = end_nodes
    orig_score = score
    start_nodes = list(set(start_nodes) & set(expressed_genes))
    end_nodes = list(set(end_nodes) & set(expressed_genes))
    score = list(set(score) & set(expressed_genes))
    if custom_end:
        orig_cend = custom_end_elems
        custom_end_elems = list(set(custom_end_elems) & set(expressed_genes))
    if custom_start:
        orig_cstart = custom_start_elems
        custom_start_elems = list(set(custom_start_elems) & set(expressed_genes))

def chunks(l, n):
    """Yield successive n-sized chunks from l."""
    for i in range(0, len(l), n):
        yield l[i:i + n]

total_no = len(start_nodes) * len(end_nodes)

# Do not calculate same short path twice (graph is undirected!)
common_nodes = list(set(start_nodes) & set(end_nodes))

# Filter PPI to expressed genes
if expressed_genes:
    print("Filtering PPI to %d genes..." % len(expressed_genes))
    ppi = ppi.subgraph(expressed_genes)
    ppigenes = ppi.nodes()
    ppiedges = ppi.edges()
    print("PPI Expression Filter (%d genes): %d nodes, %d edges" % (len(expressed_genes), len(ppigenes),
                                                                    len(ppiedges)))
    stats_fh.write("PPI Expression Filter (%d genes): %d nodes, %d edges\n" % (len(expressed_genes), len(ppigenes),
                                                                               len(ppiedges)))
    stats_fh.write("PPI Expression Filter - Start: %d nodes (%d), End: %d nodes (%d), "
                   "Score: %d nodes (%d)\n" % (len(start_nodes), len(orig_start), len(end_nodes), len(orig_end),
                                               len(score), len(orig_score)))
    if custom_start_elems:
        stats_fh.write("PPI Expression Filter - Custom Start: %d nodes (%d)\n" % (len(custom_start_elems),
                                                                                  len(orig_cstart)))
    if custom_end_elems:
        stats_fh.write("PPI Expression Filter - Custom End: %d nodes (%d)\n" % (len(custom_end_elems),
                                                                                  len(orig_cend)))
all_pairs = []
for start_node in start_nodes:
    for end_node in end_nodes:
        if start_node != end_node:
            all_pairs.append((start_node, end_node))

chunked_pairs = list(chunks(all_pairs, int(math.floor(float(len(all_pairs))/(NUM_PROCS-1)))))


def calc_shortest_path(chunk):
    all_paths = []
    for iter_no, pair in enumerate(chunk):
        if iter_no % report_per == 0:
            print("Completed %d." % iter_no)
        start = pair[0]
        end = pair[1]

        try:
            # Instead of taking a random short path, pick most relevant shortest path among
            # competing shortest paths
            common = {}
            shortest_path_len = len(nx.shortest_path(ppi, source=start, target=end))
            accepted_path = False
            for currpath in nx.shortest_simple_paths(ppi, source=start, target=end):
                if len(currpath) > shortest_path_len and accepted_path:
                    break

                curr_common = len(set(currpath[1:-1]) & set(score))
                curr_path_mirna = set(currpath[1:-1]) & set(all_mirnas)

                if len(curr_path_mirna) > 0:
                    curr_accepted_path = False
                else:
                    curr_accepted_path = True

                if curr_accepted_path:
                    accepted_path = True
                    mark_path = []
                    for item in currpath:
                        marked_item = item
                        if item in experimental_mirna + experimental_genes:
                            marked_item += "*"
                        if item in known_genes + known_mirna:
                            marked_item += "**"
                        mark_path.append(marked_item)
                    if curr_common in common:
                        common[curr_common].append((currpath, mark_path))
                    else:
                        common[curr_common] = [(currpath, mark_path)]
            if accepted_path:
                all_paths.extend(common[max(list(common.keys()))])
            else:
                raise nx.exception.NetworkXNoPath

        except nx.exception.NetworkXNoPath:
            all_paths.append(-1)

    return all_paths

print("Starting calculating shortest path between start-end node pairs "
      "(Total: %d, Max Chunk Size: %d)..." % (len(all_pairs), int(math.floor(float(len(all_pairs))/(NUM_PROCS-1)))))
p = multiprocessing.Pool(processes=NUM_PROCS)
short_paths_raw = p.map(calc_shortest_path, chunked_pairs)

short_paths = [item for sublist in short_paths_raw for item in sublist]

no_path = short_paths.count(-1)
with_path = len(short_paths)-no_path
print("-Completed with %d paths left for p-value calculation." % with_path)
print("--No path between %d out of %d (%%%.2f) start-end node can be found." % (no_path, len(short_paths),
                                                                                100 * (float(no_path) /
                                                                                       len(short_paths))))
stats_fh.write("Total calculated paths: %d\n" % with_path)
stats_fh.write("-No available paths between %d start-end node pairs\n" % no_path)

print("Saving shortest paths..")
pickle.dump(short_paths, open(outfolder+"/Shortest_Paths.dat", "wb"), pickle.HIGHEST_PROTOCOL)
print("-Completed.")
stats_fh.write("--------------------------------------------\n")
stats_fh.close()