import sys, os, glob
import numpy as np
from scipy.stats import combine_pvalues as CPVAL
from stattests import fdr_BH

args = sys.argv[1:]
if len(args) < 3: 
    print 'Usage: python <prog> <K-S regression DE file> <K-S regression kinase peptide counts file> <KSEA output directory> <HGT enrichment output directory> <output file name>'
    sys.exit()

# set arguments
inKSR = args[0]
inKSR_pfn = args[1]
inKSEA = args[2]
inHGT = args[3]
ofn = args[4]
if not inKSEA.endswith('/'): inKSEA += '/'
if not inHGT.endswith('/'):  inHGT += '/'

# options
nsubs_thresh = 3 # minimum number of substrates associated with kinase
ksea_min = 1e-5 # minimum p-value for KSEA; should be 1/nperm

# GSEA inputs
inKSEA_pos = glob.glob(inKSEA + 'gsea_report*pos*.xls')[0]
inKSEA_neg = glob.glob(inKSEA + 'gsea_report*neg*.xls')[0]

# HGT inputs
inHGT_up = glob.glob(inHGT + '*up*FDR*txt')[0]
inHGT_down=glob.glob(inHGT + '*down*FDR*txt')[0]

# Load number of substrates per kinase - count as total modifications for now
Nsubstrates = {}
for line in open(inKSR_pfn).readlines()[1:]:
    l = line.strip().split('\t')
    kin, sites = l[0], l[2]
    NS = 0
    for site in sites.split(';'):
        for mod in site.split('_')[1:]:
            NS += 1
    Nsubstrates[kin] = NS

# Get results - KS regression
RESULTS = {}                                     
EFFSIZES = {}
for line in open(inKSR).readlines()[1:]:
    l = line.strip().split('\t')               
    k, eff, pv = l[1], float(l[2]), float(l[3])                  
    if Nsubstrates[k] < nsubs_thresh: continue                           
    RESULTS[k] = [pv, 1, 1]
    EFFSIZES[k] = [eff, 0, 0]

# Get KSEA results
for line in open(inKSEA_pos).readlines()[1:]:
    l = line.strip().split('\t')         
    k, eff, pv = l[0], float(l[5]), float(l[6])   
    if pv == 0: pv = ksea_min 
    if k not in RESULTS: continue
    RESULTS[k][1] = pv
    EFFSIZES[k][1] = eff

for line in open(inKSEA_neg).readlines()[1:]:
    l = line.strip().split('\t')         
    k, eff, pv = l[0], float(l[5]), float(l[6])   
    if pv == 0: pv = ksea_min    
    if k not in RESULTS: continue
    RESULTS[k][1] = pv
    EFFSIZES[k][1] = eff

# Get HGT results
HGT_dict = {}
for line in open(inHGT_up).readlines()[1:]:
    l = line.strip().split('\t')
    k, enr, pv = l[0], float(l[5]), float(l[6])
    if k not in RESULTS: continue
    HGT_dict[k] = [[enr, pv]]
for line in open(inHGT_down).readlines()[1:]:
    l = line.strip().split('\t')
    k, enr, pv = l[0], float(l[5]), float(l[6])
    if k not in RESULTS: continue
    if k not in HGT_dict: HGT_dict[k] = []
    HGT_dict[k].append([-enr, pv]) # take negative enrichment for down-regulation
for k in HGT_dict:
    hgt = HGT_dict[k]
    if len(hgt) == 1:
        RESULTS[k][2] = hgt[0][1]
        EFFSIZES[k][2] = hgt[0][0]
    else:
        hgts = sorted(hgt, key = lambda x: x[1])
        RESULTS[k][2] = hgts[0][1]
        EFFSIZES[k][2] = hgts[0][0]

# combine results
kinases = RESULTS.keys()
for k in kinases:
    cpv = CPVAL(RESULTS[k], method = 'fisher')[1]
    RESULTS[k].append(cpv)

# FDR correction and sort results
FDR = fdr_BH(np.array([RESULTS[k][-1] for k in kinases]))
results = []
for (k, fdr) in zip(kinases, FDR):
    results.append(RESULTS[k] + [fdr])
si = np.argsort(np.array([r[-2] for r in results]))
kinases = np.array(kinases)[si]
results = np.array(results)[si]

# write output
of = open(ofn, 'w')
of.writelines('\t'.join(['Kinase','Reg_beta','KSEA_NES','HGT_enrich','Reg_pval','KSEA_pval','HGT_pval','Fisher_pval','FDR','signConflict'])+'\n')
for i, (k,r) in enumerate(zip(kinases, results)):
    ES = EFFSIZES[k]
    er, ek, eh = ES
    pr, pk, ph = results[i,0], results[i,1], results[i,2]
    fpv = results[i,3]
    fqv = results[i,4]
    signConf = 'TRUE'
    if np.sign(er) == np.sign(ek) == np.sign(eh): signConf='FALSE'

    ol = '%s\t%0.3f\t%0.3f\t%0.3f\t%0.3g\t%0.3g\t%0.3g\t%0.3g\t%0.3g\t%s\n'%(k, er, ek, eh, pr, pk, ph, fpv, fqv, signConf)
    of.writelines(ol)
of.close()


