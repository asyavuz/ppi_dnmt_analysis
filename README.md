#PPI Network Analysis Tool
===========================

This repository hosts a PPI-based gene list enrichment tool. 

* * *

Description:
---------------
This tool calculates shortest paths between desired two node lists and assesses their statistical significance using hypergeometric tests. Statistical significance is assessed through presence of validated disease-associated genes in calculated paths.  Shortest paths were then enriched with KEGG pathways.

* * *

Installation
-------------
You can get this tool with 
```text
git clone https://github.com/asyavuz/ppi_dnmt_analysis.git
``` 
command.

This tool is tested only on Mac OSX (El Capitan). If you observe any problems, please report it at the [*Issues*](https://github.com/asyavuz/ppi_dnmt_analysis/issues) section.

Please refer individual websites of required libraries/tools for their installation instructions. 

### Requirements ###
* [networkx](https://networkx.github.io)

* * *

Usage
-----
Initially, you need to create the PPI network. You can use [BioGrid](https://thebiogrid.org) links file or SIF files.

In order to generate PPI from a BioGrid links file, you may use the following command:
~~~~
python helperscripts/generate_ppi_from_BioGrid.py -i [BIOGRID_LINKS_FILE_LOCATION]
~~~~

or for SIF files you can use:

~~~~
python helperscripts/generate_ppi_from_sif.py -i [SIF_FILE_LOCATION]
~~~~

For additional parameters of these scripts please type:

~~~~
python [SCRIPT_NAME] --help
~~~~

* * * 

After generation of PPI network, you need to calculate shortest paths. You may run this part using [PyPy](http://pypy.org) as it reduces run time significantly.

Please run following command, and find out necessary parameters for your analysis:
~~~~
python generate_paths.py --help
~~~~

Lastly, you may enrich calculated shortest paths with analysis script. You can view necessary parameters of this script by
~~~~
python analyse_shortest_paths.py --help
~~~~

Please open an issue in this page if you have any questions.

* * *

License
-------
This tool uses the BSD-3 License (see [LICENCE](https://github.com/asyavuz/ppi_dnmt_analysis/blob/master/LICENSE)). 
* * *

Reference
---------
If you would like to use tool in your publications, please consider citing:
>   Yavuz, A.S. (2017). Predictive analysis of epigenetic variability (Unpublished doctoral dissertation). Sabanci University, Istanbul, Turkey.
