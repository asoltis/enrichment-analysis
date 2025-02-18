import sys, os
import numpy as np
from scipy import stats
from sklearn.linear_model import LinearRegression
from stattests import fishers_exact_test as FET
import csv
import statsmodels.api as sm

def fdr_BH(pvals):
    '''
    Compute Benjamini-Hochberg FDR q-values on sorted array of p-values.
    '''
    total = pvals.size
    si = np.argsort(pvals) 
    pvals_sort = pvals[si]
    
    fdrs0 = total*pvals_sort/range(1,total+1)
    fdrs = []
    # Preserve monotonicity
    for i in range(0,total):
        fdrs.append(min(fdrs0[i:]))
    
    # preserve original sort order 
    pairs = [(i,j) for i, j in enumerate(si)]
    fdrs_resort = np.zeros(pvals.shape)
    for pair in pairs: fdrs_resort[pair[1]] = fdrs[pair[0]]
    return fdrs_resort

## KS scores file
sidr = './'
scores_fn = sdir + 'ks_regression_score.txt'

## Expression sub-types
exp_st_fn = sdir + 'sample_expression_subtypes.txt'

# Sexes, ancestries
sex_fn = sdir + 'sample_sexes.txt'
anc_fn = sdir + 'sample_ancestries.txt'

##########
# OUTPUT #
##########
outdir = '/'.join(scores_fn.split('/')[0:-1]) + '/DE_scores_by_exp_subtypes/'
os.system('mkdir -p %s'%(outdir))

# load subtype info
clusters = {}
clust_nums = set()
for line in open(exp_st_fn).readlines():
    l = line.strip().split('\t') 
    typD = {'PI':0,'PP':1,'TRU':2}
    s, typ = l[0], typD[l[1]]
    clusters[s] = typ
    clust_nums.add(typ)
tot_clusters = len(clust_nums)
clust_samps = sorted(clusters.keys())
tot_clust_samps = len(clusters)

# Make samples by cluster binary matrix
clust_bin_mat = np.zeros((tot_clust_samps, tot_clusters))
for i, s in enumerate(clust_samps):
    cn = clusters[s] 
    clust_bin_mat[i, cn] = 1
clust_bin_mat = np.array(clust_bin_mat)

# covariate data
samp_data = {}
for line in open(sex_fn).readlines():
    l = line.strip().split('\t')
    s, sex = l
    samp_data[s] = [sex]

# Load kinase score data
kinases, ksamps, kmat = [], [], []
for i, line in enumerate(open(scores_fn).readlines()):
    l = line.strip().split('\t')
    if i == 0:
        ksamps = l[1:]
    else:
        kin = l[0]
        v = map(float, l[1:])
        kinases.append(kin)
        kmat.append(v)
kmat = np.array(kmat)
ksamps = np.array(ksamps)
kinases = np.array(kinases)

# sort protein info to match cluster data; also, remove samples with no cluster assignment
psi = [0] * len(clust_samps)
for i, s in enumerate(ksamps):
    if s not in clust_samps: continue
    else:
        ind = np.argwhere(np.array(clust_samps) == s)[0][0]
        psi[ind] = i
ksamps = ksamps[psi]
kmat = kmat[:, psi].T

# Get covariate data
COV_DATA = []
for s in ksamps:
    #sex, anc = samp_data[s]
    sex = samp_data[s][0]
    SEX, ANC = None, None
    
    if sex == 'female': SEX = 1
    else: SEX = 0

    #if anc == 'EUR': ANC = 0
    #else: ANC = 1
 
    COV_DATA.append([SEX])
COV_DATA = np.array(COV_DATA)

############
# ANALYSIS #
############

# Test for sig proteintides
sig_kinases_pvals = np.ones((len(kinases), tot_clusters))
sig_kinases_betas = np.ones((len(kinases), tot_clusters))
non_cases = []
for cn in range(0, tot_clusters):
    clust_inds = clust_bin_mat[:, cn] #.reshape(-1,1)
    
    # Compute significance for each protein in cluster vs. others
    for pi in range(0, len(kinases)):
        protein = kinases[pi]
        Yp = kmat[:, pi].reshape(-1,1)
        
        # standardize
        #Ys = Yg - np.mean(Yg)
        #if np.std(Yg) > 0: Ys /= np.std(Yg)
        #Xs = clust_inds - np.mean(clust_inds)
        #if np.std(clust_inds) > 0: Xs /= np.std(clust_inds)
   
        Xs = np.zeros((COV_DATA.shape[0], 1+COV_DATA.shape[1]))
        Xs[:,0] = clust_inds
        for xi in range(0, COV_DATA.shape[1]):
            Xs[:, xi+1] = COV_DATA[:, xi]
        Ys = Yp
    
        nclust = len(Ys[(Xs[:,0]==1)])
        nnon = len(Ys[(Xs[:,0]==0)])

        # Handle missing data
        if sum(np.isnan((Ys[(Xs[:,0]==1)]))) >= (0.8 * nclust): # skip if 80% or more of cluster samples are missing a measurement
            sig_kinases_pvals[pi, cn] = 1.0
            sig_kinases_betas[pi, cn] = 0.0
            continue
        elif sum(np.isnan((Ys[(Xs[:,0]==0)]))) == nnon: # skip if all non-cluster cases have no measurement; save info
            non_cases.append([pi, cn])
            sig_kinases_pvals[pi, cn] = 1.0
            sig_kinases_betas[pi, cn] = 0.0
            continue
        
        # remove nan data for passing kinases
        notnani = ~np.isnan(Ys)
        Ysf = Ys[notnani].reshape(-1,1)
        Xsf0 = Xs[np.argwhere(notnani==True)[:,0],:]

        # add intercept
        Xsf = np.ones((Xsf0.shape[0], 1+Xsf0.shape[1]))
        for xi in range(0, Xsf0.shape[1]):
            Xsf[:, xi+1] = Xsf0[:, xi]

        # Run regression
        clf = LinearRegression(fit_intercept = False)
        clf.fit(Xsf, Ysf)
       
        SSE = np.sum((clf.predict(Xsf) - Ysf) ** 2, axis = 0) / float(Xsf.shape[0] - Xsf.shape[1])
        SE = np.array([np.sqrt(np.diagonal(SSE[i] * np.linalg.inv(np.dot(Xsf.T, Xsf)))) for i in range(SSE.shape[0])])
        tstat = clf.coef_ / SE
        pvals = 2 * (1-stats.t.cdf(np.abs(tstat), Ysf.shape[0] - Xsf.shape[1])) 
        pv = pvals[0, 1]
        beta = clf.coef_[0, 1]
        
        if np.isnan(pv): pv = 1.0
        #if beta < 0: pv = 1.0 # want scores to be higher in cluster
        sig_kinases_pvals[pi, cn] = pv
        sig_kinases_betas[pi, cn] = beta

print non_cases

# get FDR p-values for each cluster
sig_kinases_fdr = np.ones(sig_kinases_pvals.shape)
for cn in range(0, tot_clusters):
    sig_kinases_fdr[:, cn] = fdr_BH(sig_kinases_pvals[:, cn])

##########
# OUTPUT #
##########
# write significance results
fdr_thresh = 0.1
for cn in range(0, tot_clusters):
    # set up output
    if cn == 0: cname = 'PI'
    elif cn == 1: cname = 'PP'
    elif cn == 2: cname = 'TRU'
    print cn, cname
    ofn = outdir + 'subtype_%s_DE_kinase_scores_FDR_lt_%s.txt'%(cname, str(fdr_thresh))
    of = open(ofn, 'w')
    of.writelines('Subtype\tKinase\tbeta\tpval\tFDR\n')
    
    ofn2 = outdir + 'subtype_%s_DE_kinase_scores_all.txt'%(cname)
    of2 = open(ofn2, 'w')
    of2.writelines('Subtype\tKinase\tbeta\tpval\tFDR\n')

    olines, olines2 = [] ,[]
    for pi in range(0, len(kinases)):
        fdr = sig_kinases_fdr[pi, cn]
        protein = kinases[pi]
        beta = sig_kinases_betas[pi, cn]
        pv = sig_kinases_pvals[pi, cn]
        
        if fdr < 0.1:
            olines.append([cname, protein, beta, pv, fdr])
        olines2.append([cname, protein, beta, pv, fdr])

    # sort output lines and write
    olines.sort(key = lambda x: x[-2])
    for ol in olines:
        of.writelines('%s\t%s\t%0.3f\t%0.3g\t%0.3g\n' % (tuple(ol)))
    of.close()

    olines2.sort(key = lambda x: x[-2])
    for ol in olines2:
        of2.writelines('%s\t%s\t%0.3f\t%0.3g\t%0.3g\n' % (tuple(ol)))
    of2.close()

# write all observed kinases file
of = open(outdir + 'all_observed_kinases.txt', 'w')
pset = sorted(list(set([x for x in kinases])))
for p in pset: of.writelines('%s\n'%(p))
of.close()

