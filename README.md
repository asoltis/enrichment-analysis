Repository for Gene ontology and other pathway enrichment analysis scripts

## GO directory
- ontology_scripts sub-directory
    - parseOBO_and_map.py script downloads and parses GO OBO file (ontology graph) and GAF (gene annotations) from geneontology.org
    - New version uses urllib to obtain files and obonet package to create graph from OBO file
    - legacy/ folder contains old versions of code
- GO_enrichment sub-directory
    - GO_enrichment.py is script for running GO enrichment analysis of foreground against background gene lists
    - Cytoscape graph creation/style files
        - create_cytoscape_files_from_GO_output.py - create graph files from GO output files
        - GOstyle.xml - Cytoscape style file for GO graph data

## pathway_enrichments directory
- includes "pathway_enrichment.py" for executing pathway enrichment analysis of foreground versus background
  gene lists against gene sets (valid sets listed in accompanying file). Currently reads MSigDB *.gmt files
  to load/create gene sets for analysis. 
- "valid_gene_sets.txt" file lists (current) valid sets to choose from

## stattests directory
- helper statistical functions for enrichment analysis routines

