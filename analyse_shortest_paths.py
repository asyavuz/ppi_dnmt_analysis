#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
PPI-based gene-list enrichment tool.
P-value calculation and KEGG enrichment part.
"""
# Author: Ahmet Sinan Yavuz <asinanyavuz@sabanciuniv.edu>
# BSD-3 License

import os
import sys
import math
import pickle
import warnings
import argparse
import scipy.stats
import numpy as np
import networkx as nx
import urllib.request
import multiprocessing
from builtins import range
from time import gmtime, strftime

orig_cmd = " ".join(sys.argv)
NUM_PROCS = multiprocessing.cpu_count()

parser = argparse.ArgumentParser(description='PPI-based gene-list enrichment tool - path evaluation and enrichment')
parser.add_argument('-m', '--mirna', help='List of important miRNAs identified by the experiment', default="")
parser.add_argument('-lg', '--litgene', help='Known gene list associated with the condition',
                    default="./data/LiteratureGenes.txt")
parser.add_argument('-lm', '--litmirna', help='Known miRNA list associated with the condition',
                    default="./data/LiteraturemiRNAs.txt")
parser.add_argument('-p', '--ppi', help='PPI network file (only pickled networkx graphs)',
                    default="./data/BioGrid.dat")
parser.add_argument('-sc', '--score', help="Score nodes in list (allowed: expt, known, both, default: both)",
                    default="both")
parser.add_argument('-min', '--minnode', help="Minimum node in a path required for pathway enrichment (including "
                                              "start and end nodes, default: 10)",
                    default=10)
parser.add_argument('-cp', '--correctpath', help="Correction method for p-value calculation in path evaluation "
                                                 "(allowed: fdr, bonferroni, default: fdr)",
                    default="fdr")
parser.add_argument('-ce', '--correctenr', help="Correction method for p-value calculation in pathway enrichment "
                                                "(allowed: fdr, bonferroni, default: fdr)",
                    default="fdr")
parser.add_argument('-pp', '--pvalpath', help="Corrected p-value cut off for path evaluation (default: 0.01)",
                    default=0.01)
parser.add_argument('-pe', '--pvalenr', help="Corrected p-value cut off for pathway enrichment (default: 0.01)",
                    default=0.01)
parser.add_argument('-t', '--tarbase', help='TARBASE database file', default="./data/tarbase_data.csv")
parser.add_argument('-tf', '--tarbasefilter', help='TARBASE filter string (separated by ;, default: '
                                                   '\'species=Homo sapiens\')',
                    default="species=Homo sapiens")
parser.add_argument('-bg', '--background', help='Background gene list for pathway enrichment '
                                                '(default: blank [for Human - KEGG])', default="")
parser.add_argument('-org', '--organism', help="KEGG organism code (default: hsa)", default="hsa")
parser.add_argument('-on', '--online', help='Skip cache and get latest KEGG pathways online', action='store_true')
parser.add_argument('-nc', '--cores', help='Number of processes (default: %d)' % NUM_PROCS, default=NUM_PROCS)
requiredNamed = parser.add_argument_group('required named arguments')
requiredNamed.add_argument('-sp', '--paths', help='Pickled shortest paths', required=True)
requiredNamed.add_argument('-g', '--gene', help='List of important genes identified by the experiment', required=True)
args = parser.parse_args()

# Input files
ppi_network = args.ppi
tarbase_file = args.tarbase
known_genes_file = args.litgene
known_mirna_file = args.litmirna
identified_genes_file = args.gene
identified_mirna_file = args.mirna
shortest_path_file = args.paths

# TARBASE filter, only equality filters are allowed.
tmp = args.tarbasefilter
tarbase_filter = tmp.split(";")

# P-value filters
correct = args.correctpath
p_value_filter = float(args.pvalpath)

enrichment_correct = args.correctenr
p_value_filter_enrichment = float(args.pvalenr)

# Background gene list
custom_bg_file = args.background

# Advanced controls
online = args.online
kegg_cache = "./data/KEGG_%s.db" % str(args.organism)

# ################################################ #
#       DO NOT TOUCH BELOW HERE UNLESS NEEDED      #
# ################################################ #

if args.cores:
    NUM_PROCS = int(args.cores)


def bonferroni_correction(pvals, alpha):
    # Bonferroni correction

    cutoff = alpha / float(len(pvals))
    reject = pvals <= cutoff
    pvals_corrected = np.minimum(np.ones(len(pvals)), pvals * float(len(pvals)))

    return reject, pvals_corrected


def fdr_correction(pvals, alpha):
    # Benjamini/Hochberg p-value correction for false discovery rate
    # Code is directly obtained from statsmodels package.
    # Included here for decoupling dependency on statsmodels.

    pvals_sortind = np.argsort(pvals)
    pvals_sorted = np.take(pvals, pvals_sortind)
    nobs = len(pvals_sorted)
    ecdffactor = np.arange(1, nobs + 1) / float(nobs)

    reject = pvals_sorted <= ecdffactor * alpha
    if reject.any():
        rejectmax = max(np.nonzero(reject)[0])
        reject[:rejectmax] = True

    pvals_corrected_raw = pvals_sorted / ecdffactor
    pvals_corrected = np.minimum.accumulate(pvals_corrected_raw[::-1])[::-1]
    pvals_corrected[pvals_corrected > 1] = 1

    # reorder p-values and rejection mask to original order of pvals
    pvals_corrected_ = np.empty_like(pvals_corrected)
    pvals_corrected_[pvals_sortind] = pvals_corrected

    reject_ = np.empty_like(reject)
    reject_[pvals_sortind] = reject

    return reject_, pvals_corrected_


def chunks(l, n):
    """Yield successive n-sized chunks from l."""
    for i in range(0, len(l), n):
        yield l[i:i + n]


if correct not in ["fdr", "bonferroni"]:
    raise ValueError("Unsupported p-value correction method: %s." % correct)

outfolder = os.path.dirname(os.path.abspath(shortest_path_file))

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

if (known_mirna_file or identified_mirna_file) and tarbase_file:
    mirna_targets = {}
    print("Loading TARBASE data (Filters: %s)..." % (" AND ".join(tarbase_filter)))
    with open(tarbase_file) as f:
        header = {}
        total = 0
        unfiltered_total = 0
        for i, line in enumerate(f):
            if i == 0:
                elems = line.split("\t")
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

print("Loading known genes...")
unrecognised = 0
total = 0
with open(known_genes_file) as f:
    for i, line in enumerate(f):
        total += 1
        try:
            b = ppigenes.index(line.strip().upper())
            known_genes.append(line.strip().upper())
        except ValueError:
            unrecognised += 1

known_genes = list(set(known_genes))
print("-Loaded. Recognised %d unique genes as in PPI out of %d entries." % (total - unrecognised, total))
stats_fh.write("Known Genes: %s (%d Entries, %d Unique Genes, %d recognised as in PPI)\n" % (known_genes_file, total,
                                                                                             len(known_genes),
                                                                                             total - unrecognised))
if known_mirna_file:
    print("Loading known miRNAs...")
    unrecognised = 0
    total = 0
    with open(known_mirna_file) as f:
        for i, line in enumerate(f):
            total += 1
            try:
                b = ppigenes.index(line.strip().upper())
                known_mirna.append(line.strip().upper())
            except ValueError:
                unrecognised += 1
    known_mirna = list(set(known_mirna))
    print("-Loaded. Recognised %d unique miRNAs as in PPI out of %d entries." % (total - unrecognised, total))
    stats_fh.write(
        "Known miRNAs: %s (%d miRNAs, %d Unique miRNAs, %d recognised as in PPI)\n" % (known_mirna_file, total,
                                                                                       len(known_mirna),
                                                                                       total - unrecognised))
print("Loading experimental genes...")
unrecognised = 0
total = 0
with open(identified_genes_file) as f:
    for i, line in enumerate(f):
        total += 1
        try:
            b = ppigenes.index(line.strip().upper())
            experimental_genes.append(line.strip().upper())
        except ValueError:
            unrecognised += 1

print("-Loaded. Recognised %d genes out of %d entries." % (unrecognised, total))
stats_fh.write("Experimental Genes: %s (%d Genes, %d recognised as in PPI)\n" % (identified_genes_file, total,
                                                                                 total - unrecognised))
if identified_mirna_file:
    print("Loading experimental miRNAs...")
    unrecognised = 0
    total = 0
    with open(identified_mirna_file) as f:
        for i, line in enumerate(f):
            total += 1
            try:
                b = ppigenes.index(line.strip().upper())
                experimental_mirna.append(line.strip().upper())
            except ValueError:
                unrecognised += 1

    print("-Loaded. Recognised %d miRNAs out of %d entries." % (unrecognised, total))
    stats_fh.write("Experimental miRNAs: %s (%d miRNAs, %d recognised as in PPI)\n" % (identified_mirna_file, total,
                                                                                       total - unrecognised))

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

bg_genes = []

score = []
if args.score == "expt":
    score = experimental_genes + experimental_mirna
elif args.score == "known":
    score = known_genes + known_mirna
elif args.score == "both":
    score = list(set(experimental_genes + experimental_mirna + known_genes + known_mirna))

# Start processing shortest paths.
short_paths = pickle.load(open(shortest_path_file, "rb"))

no_path = short_paths.count(-1)
with_path = len(short_paths) - no_path
print("-Completed with %d paths left for p-value calculation." % with_path)
print("--No path between %d out of %d (%%%.2f) start-end node can be found." % (no_path, len(short_paths),
                                                                                100 * (float(no_path) /
                                                                                       len(short_paths))))
stats_fh.write("Total calculated paths: %d\n" % with_path)
stats_fh.write("-No available paths between %d start-end node pairs\n" % no_path)

no_of_background = len(ppi.nodes())
no_of_known_in_ppi = len(set(ppi.nodes()) & set(score))


def calculate_pval(chunk):
    returnset = []
    for item in chunk:
        curr_path = item[0]
        no_of_hotpath = len(set(curr_path) & set(score))
        p_val = scipy.stats.hypergeom.sf(no_of_hotpath, no_of_background, no_of_known_in_ppi, len(curr_path) - 2)

        returnset.append((item[0], item[1], p_val))
    return returnset


short_paths = [item for item in short_paths if isinstance(item, tuple)]
print("Calculating p-values...")

pval_chunks = list(chunks(short_paths, int(math.floor(len(short_paths) / float((NUM_PROCS - 1))))))

p2 = multiprocessing.Pool(processes=NUM_PROCS)
short_paths_raw = p2.map(calculate_pval, pval_chunks)

short_paths = [item for sublist in short_paths_raw for item in sublist]

corrected_p_values = None
p_value_mask = None
print("-Correcting p-values with %s..." % correct)
pvals = np.asarray([item[2] for item in short_paths], dtype=np.float64)

if correct == "bonferroni":
    p_value_mask, corrected_p_values = bonferroni_correction(pvals, p_value_filter)
    p_value_mask2, corrected_p_values2 = bonferroni_correction(pvals, p_value_filter)
elif correct == "fdr":
    p_value_mask, corrected_p_values = fdr_correction(pvals, p_value_filter)
print("-Done.")

short_paths = [(item[0], item[1], item[2], corrected_p_values[k]) for k, item in enumerate(short_paths)]

sorted_ind = np.argsort(corrected_p_values)
print("Writing out paths file...")
written = 0
with open(outfolder + "/Paths.txt", "w") as o:
    o.write("\t".join(["#", "Start Point", "End Point", "Effective Path Size", "Hot Node in Path", "PPI Size",
                       "Hot Node in PPI", "P-value", "%s P-value" % correct, "Marked Path", "Raw Path"]) + "\n")
    for path_no, ind in enumerate(sorted_ind):
        path = short_paths[ind][0]
        no_of_known_in_path = len(set(path[1:-1]) & set(score))
        marked_path = short_paths[ind][1]
        pval = short_paths[ind][2]
        cor_pval = short_paths[ind][3]
        if cor_pval <= p_value_filter:
            o.write("\t".join([str(path_no), path[0], path[-1], str(len(path) - 2), str(no_of_known_in_path),
                               str(len(ppi.nodes())), str(no_of_known_in_ppi),
                               "%g" % pval, "%g" % cor_pval,
                               ",".join(marked_path), ",".join(path)]) + "\n")
            written += 1
print("-Done. Written %d paths that satisfy %g p-value threshold." % (written, p_value_filter))
stats_fh.write("-Paths satisfying %g p-value threshold: %d\n" % (p_value_filter, written))
stats_fh.write("\n")
stats_fh.write("Enrichment with KEGG")
print("Starting enrichment analyses for %d paths passing p-value threshold..." % written)
enrichment_source = {}

if kegg_cache:
    try:
        enrichment_source = pickle.load(open(kegg_cache, "rb"))
    except FileNotFoundError:
        warnings.warn("No KEGG cache found. Running via KEGG REST API...",
                      RuntimeWarning)
        online = True

if online:
    print("Loading KEGG pathways from KEGG REST API...")
    kegg_pathways = {}

    # First let's get current KEGG information
    resp = urllib.request.urlopen("http://rest.kegg.jp/info/pathway")
    kegg_info = ""
    for line in resp:
        kegg_info += line.decode("utf-8")
    kegg_pathways["DBINFO"] = kegg_info

    # Now retrieve list of pathways
    resp = urllib.request.urlopen("http://rest.kegg.jp/list/pathway/%s" % str(args.organism))
    for line in resp:
        line = line.decode("utf-8").strip()
        elems = line.split("\t")
        path_id = elems[0].split(":")[1]
        path_name = elems[1].split(" - ")[0]
        kegg_pathways[path_id] = {"Name": path_name, "Genes": []}
    print("-Found %d pathways for organism %s..." % (len(kegg_pathways), str(args.organism)))

    all_kegg_genes = []
    # Retrieve genes for each of the pathway
    for path_id in kegg_pathways:
        if path_id == "DBINFO":
            continue
        print("-Fetching pathway: %s" % path_id)
        resp = urllib.request.urlopen("http://rest.kegg.jp/link/genes/%s" % path_id)
        for line in resp:
            line = line.decode("utf-8").strip()
            elems = line.split("\t")
            gene_id = elems[-1].strip().split(":")[1]
            all_kegg_genes.append(gene_id)
            try:
                kegg_pathways[path_id]["Genes"].append(entrez_to_hgnc[gene_id])
            except KeyError:
                warnings.warn("Cannot find gene ID conversion for gene ID %s. Instead ID is used." % gene_id)
                kegg_pathways[path_id]["Genes"].append(gene_id)
    if kegg_cache:
        pickle.dump(kegg_pathways, open(kegg_cache, "wb"), pickle.HIGHEST_PROTOCOL)
        with open("./data/KEGG_%s.txt" % str(args.organism), "w") as o:
            for curgene in list(set(all_kegg_genes)):
                try:
                    o.write(entrez_to_hgnc[curgene] + "\n")
                except KeyError:
                    warnings.warn("Cannot find gene ID conversion for gene ID %s. Instead ID is used." % curgene)
                    o.write(curgene + "\n")

    enrichment_source = kegg_pathways
    print("-KEGG pathways loaded.")

stats_fh.write("KEGG Info:\n")
stats_fh.write(enrichment_source["DBINFO"] + "\n\n")

if args.background:
    stats_fh.write("Selected Enrichment Background: %s\n" % args.background)
    print("Loading custom background genes...")
    unrecognised = 0
    total = 0
    with open(args.bg) as f:
        for i, line in enumerate(f):
            total += 1
            bg_genes.append(line.strip())

    print("-Loaded %d genes as background set." % total)
    stats_fh.write("Custom Enrichment Background File: %s (%d entries)\n" % (custom_bg_file, total))
else:
    stats_fh.write("Selected Enrichment Background: KEGG_%s\n" % str(args.organism))
    print("Loading human KEGG background set...")
    total = 0
    with open("./data/KEGG_%s.txt" % str(args.organism)) as f:
        for line in f:
            bg_genes.append(line.strip())
            total += 1
    print("-Loaded %d genes as background set." % total)
    stats_fh.write("Enrichment Background File: %s (%d entries)\n" % ("./data/Human_KEGG.txt", total))


def enrich(chunk):
    curr_enrichments = []
    filtered = 0
    for item in chunk:
        curr_ind = item[0]
        curr_path_no = item[1]

        curr_path = short_paths[curr_ind][0]

        # Skip path if it does not contain minimum required amount of nodes
        if len(curr_path) < int(args.minnode):
            filtered += 1
            continue

        curr_marked_path = short_paths[curr_ind][1]
        pathways = []
        enrichment_raw_pvals = []

        # Run enrichment test for each pathway in the enrichment source
        for path_id in enrichment_source:
            if path_id == "DBINFO":
                continue

            common_genes = len(set(curr_path) & set(enrichment_source[path_id]["Genes"]))
            path_genes = len(enrichment_source[path_id]["Genes"])

            # Skip pathway if there's no common gene
            if common_genes == 0:
                continue

            pathways.append(path_id)
            path_pval = scipy.stats.hypergeom.sf(common_genes, len(bg_genes), path_genes, len(curr_path))
            enrichment_raw_pvals.append(path_pval)

        # Correct p-values
        path_mask = None
        enrichment_corrected_p_vals = None
        if enrichment_correct == "bonferroni":
            path_mask, enrichment_corrected_p_vals = bonferroni_correction(np.asarray(enrichment_raw_pvals,
                                                                                      dtype=np.float64),
                                                                           p_value_filter_enrichment)
        elif enrichment_correct == "fdr":
            path_mask, enrichment_corrected_p_vals = fdr_correction(np.asarray(enrichment_raw_pvals,
                                                                               dtype=np.float64),
                                                                    p_value_filter_enrichment)

        # Add filtered enriched terms to general enrichment list
        sorted_ind = list(np.argsort(enrichment_corrected_p_vals))
        for pathway_no in sorted_ind:
            if path_mask[pathway_no]:
                common_genes = list(set(curr_path) & set(enrichment_source[pathways[pathway_no]]["Genes"]))
                genes_entrez = "+".join([hgnc_to_entrez[gene] for gene in common_genes])
                kegg_link = "http://www.genome.jp/kegg-bin/show_pathway?%s+%s" % (pathways[pathway_no],
                                                                                  genes_entrez)
                curr_enrichments.append((pathways[pathway_no],
                                         enrichment_source[pathways[pathway_no]]["Name"],
                                         len(enrichment_source[pathways[pathway_no]]["Genes"]),
                                         len(curr_path),
                                         len(set(curr_path) & set(
                                             enrichment_source[pathways[pathway_no]]["Genes"])),
                                         "%g" % enrichment_raw_pvals[pathway_no],
                                         "%g" % enrichment_corrected_p_vals[pathway_no], curr_path_no,
                                         ",".join(curr_marked_path),
                                         ",".join(curr_path), kegg_link))
    return curr_enrichments, filtered


enrich_pairs = []

for path_no, ind in enumerate(sorted_ind):
    if not p_value_mask[ind]:
        continue
    enrich_pairs.append((ind, path_no))

enrich_chunks = list(chunks(enrich_pairs, int(math.floor(len(enrich_pairs) / float((NUM_PROCS - 1))))))

p3 = multiprocessing.Pool(processes=NUM_PROCS)
enrichments_raw = p3.map(enrich, enrich_chunks)

enrichments = []
total_filtered = 0
for res in enrichments_raw:
    total_filtered += res[1]
    for sublist in res[0]:
        enrichments.append(sublist)

print(len(enrichments))
# enrichments = [item for sublist in enrichments_raw for item in sublist]

enrichments.sort(key=lambda tup: tup[7])
print("-Done.")

print("Writing enrichment results...")
written = 0
with open(outfolder + "/Enriched_Terms.txt", "w") as o:
    o.write("Pathway ID\tPathway Name\tTerm\tQuery\tCommon\tP-value\t%s P-value\tPath #"
            "\tMarked Path\tRaw Path\tMarked KEGG Link\n" % enrichment_correct)
    for enrichment in enrichments:
        if enrichment[7] <= p_value_filter_enrichment:
            o.write("\t".join(list(map(str, enrichment))) + "\n")
            written += 1
print("-Done.")
stats_fh.write("Paths filtered due to minimum node threshold: %d (Minimum Node Threshold: %d)\n" % (total_filtered,
                                                                                                    int(args.minnode)))
stats_fh.write("Pathways enriched for all filtered paths filtered with enrichment "
               "p-value threshold of %g: %d \n" % (p_value_filter_enrichment, written))

print("Analysis completed. Results can be found in %s directory." % outfolder)
stats_fh.write("--------------------------------------------\n")
stats_fh.close()
