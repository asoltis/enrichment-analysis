# Kinase-substrate enrichment analysis code

This repository contains code for performing three kinase-substrate enrichment analysis procedures.

As implemented in [Soltis et al. Cell Rep Med 2022](https://pmc.ncbi.nlm.nih.gov/articles/PMC9729884/pdf/main.pdf):
- Methods described in "Phosphoprotein and kinase enrichment analysis" section
- Results from procedures are presented in Figure 5 (panel A in particular) and Figure S7A

"KS_regression" method
-----------------------
Runs kinase-substrate regression procedure with either standard OLS or Ridge regression. Phosphosite abundances for a given sample 
(median-centered across all samples) are regressed against a binary matrix of kinase-phosphosite links, producing effect size (beta) and
p-value estimates for each kinase's activity. These values or derived scores can then be used to estimate enrichments between sample 
groups (e.g. by subtypes). 

Scripts:
- ks_regression.py: code to run kinase regressions 
- DE_ks_scores_by_exp_subtypes.py: sample script for performing differential analysis of kinase scores between groups of samples

"KSEA" method
----------------
Runs Gene Set Enrichment Analysis (GSEA) statistical tool on pre-ranked set of phosphosites (e.g. by signed -log10(p-value))
against kinase-substrate links database (in GMT format; e.g. data from PhosphoSitePlus). Output is GSEA permutation-based statistical
enrichments for kinases based on ranked phosphosite data. 

Scripts:
- run_ksea.sh: main driver script that runs GSEA command line tool on input data
- prep_MS_for_GSEA.py: prep script run by above script to prepare differential phosphosite data for input to GSEA (ranked sites file
  with scores calculated by selected metric)

"HGT_enrichment" method
------------------------
Performs enrichment analysis using hypergeometric testing of pre-defined foreground (e.g. differential in a subtype) phosphosites against 
a background of all measured sites, considering kinase-substrate links as the database. 

Scripts:
- PTM_enrichment.py: code to perform actual enrichment analysis
- prep_MS_for_KS_enrichment.py: input prep script for running analysis with above code

Additional scripts:
--------------------
combine_ksr_ksea_ptmHGT.py
- combines enrichment p-values from above three methods with Fisher's method
- Reports combined p-value and flags sign conflicts (e.g. if effect size direction differs between methods)


