from __future__ import division, print_function
import sys, os, csv
if sys.version_info.major == 2:
    import cPickle 
elif sys.version_info.major == 3:
    import _pickle as cPickle
import numpy as np
from scipy import stats
from sklearn.preprocessing import StandardScaler
from sklearn.linear_model import LinearRegression, Ridge, Lasso, ElasticNet
from sklearn.model_selection import ShuffleSplit, cross_val_score
from sklearn.impute import KNNImputer
import warnings
import multiprocessing as mp
import statsmodels.api as SM

def warn(*args, **kwargs):
    pass

def run_regression(X, Y, regType, set_alpha = None, nrep = 20):
    '''
    Run various regression models on X and Y data
    '''
    
    # OLS
    if regType == 'OLS':
        clf = LinearRegression(fit_intercept = False)
        clf.fit(X, Y)
        SSE = np.sum((clf.predict(X) - Y) ** 2, axis = 0) / float(X.shape[0] - X.shape[1])
        SE = np.array([np.sqrt(np.diagonal(SSE[i] * np.linalg.inv(np.dot(X.T, X)))) for i in range(SSE.shape[0])])
        tstat = clf.coef_ / SE
        pvals = 2 * (1-stats.t.cdf(np.abs(tstat), Y.shape[0] - X.shape[1]))
        beta = clf.coef_[0]
        pvals = pvals[0]
        Rsquare = clf.score(X, Y)
  
    # Ridge
    elif regType == 'Ridge':

        alphas = ([pow(10, x) for x in [x*1e-2 for x in range(-400, 200, 5)]])
        train_res = []
        if set_alpha == None:
            for alpha in alphas:
                clf_cv = Ridge(alpha = alpha, fit_intercept = False)
                cv = ShuffleSplit(n_splits = nrep, test_size = 0.3, random_state = 0)
                scores = cross_val_score(clf_cv, X, Y, cv = cv, scoring = 'neg_mean_squared_error')
                train_res.append([alpha, np.mean(-scores)])
            train_res.sort(key = lambda x: x[1], reverse=False)
        
        # opt fit
        if set_alpha == None: 
            alpha = train_res[0][0]
        else: alpha = set_alpha
        clf = Ridge(alpha = alpha, fit_intercept = False)
        clf_red = Ridge(alpha = alpha, fit_intercept = False)
        
        clf.fit(X, Y)
        beta = clf.coef_[0]
        Rsquare = clf.score(X, Y)
        Ypred = clf.predict(X)
        RSSf = np.sum((Ypred - Y) ** 2)
        np_full = X.shape[1]
        df_full = X.shape[0] - np_full

        # compute F-test p-values for each coefficient
        pvals = []
        clf_for_pv = Ridge(alpha = alpha, fit_intercept = False)
        for i, b in enumerate(beta):
            inds = [j for j in range(0, len(beta)) if j != i]
            gx = X[:, inds]
            
            clf_for_pv.fit(gx, Y)
            Ypred_pv = clf_for_pv.predict(gx)
            RSSr = np.sum((Ypred_pv - Y) ** 2)
            npr = gx.shape[1]
            dpr = gx.shape[0] - npr
            
            fs = ((RSSr - RSSf) / (np_full - npr)) / (RSSf / df_full)
            if fs < 0: fpv = 1.0
            else: fpv = 1 - stats.f.cdf(fs, np_full - npr, df_full)
            pvals.append(fpv)
    
    
    return beta, pvals


def main():
    warnings.warn = warn

    ## OPTIONS ##
    tot_norm = False # NOT IMPLEMENTED
    regType = 'Ridge'
    set_alpha = 0.1 # set to None if not using
    nrep = 20
    do_impute = True
    kneighbor = 10
    measured_fraction = 0.5 # fraction of samples with measurement to keep; set to None if not using 
    normalization = 'medcenter' # set as 'rZscore' for robust Z-score or 'medcenter' for median centering
    min_peptides = 3 # Minimum number of peptides linked to kinase to include in regression; set to None if not using
    includeNKIN = False # NOT IMPLEMENTED

    ## INPUT ## 
    # Proteomics files
    phos_table = 'CPTAC_LUAD/Gillette_etal_2020_data/supp_files/1-s2.0-S0092867420307443-mmc3_TableS3B.csv'

    # PhosphoSitePlus data
    psp_fn = 'PhosphoSitePlus/parsed_files/substrate_sites_to_kinase_dict.human_to_human.20190304.pkl'
    ps_dict = cPickle.load(open(psp_fn, 'rb'))

    # Include NetworKIN predictions
    if includeNKIN: 
        nkin_fn = ''
        for line in open(nkin_fn).readlines()[1:]:
            l = line.strip().split('\t')
            subs, kin = l[0], l[2]
            prot, mod = subs.split('_')
            if prot not in ps_dict: ps_dict[prot] = {}
            if mod not in ps_dict[prot]: ps_dict[prot][mod] = []
            if kin not in ps_dict[prot][mod]: ps_dict[prot][mod].append(kin)

    ## OUTPUT ##
    adate = 'date'
    if includeNKIN:
        dbtype = 'PSP_NetworKIN_KS_interactions'
    else:
        dbtype = 'PSP_KS_interactions'
    if tot_norm:
        vals_type = 'norm_to_global_MS/phos_MS'
    else:
        vals_type = 'raw_phospho/phos_MS'
    if measured_fraction == None:
        meas_string = 'all_peps' 
    else:
        meas_string = 'measured_fraction_%s'%(str(measured_fraction))
    meas_string += '_%s'%(normalization)
    if do_impute:
        outdir = '../results/%s/%s/%s_KNNImpute_k_%d_%s/' % (adate, dbtype, vals_type, kneighbor, meas_string)
    else:
        outdir = '../results/%s/%s/%s_native_values_%s/' % (adate, dbtype, vals_type, meas_string)
    if set_alpha != None:
        outdir += '%s_regression_setAlpha_%s'%(regType, str(set_alpha))
    else:
        outdir += '%s_regression_nrep_%d'%(regType, nrep)
    if min_peptides == None:
        outdir += '_all_peptides/'
    else:
        outdir += '_ge_%d_peptides/'%(min_peptides)
    os.system('mkdir -p %s'%(outdir))
    print(outdir)
    
    ## LOAD DATA ##
    # Load phospho data - MS
    peps, psamps, pmat = [], [], []
    stypes, nmfc, expST, sexes, ancestry = [], [], [], [], []
    sindex = 23 # starting index for patient value measurements
    tsinduse = [] # list for storing column indices for tumor samples
    with open(phos_table) as fn:
        reader = csv.reader(fn, delimiter = ',', quotechar = '"')
        for li, l in enumerate(reader):
            if li < 2: continue
            if l[0] == '': continue

            # sample info rows
            if li < 71:
                if l[0] == 'Sample.ID':
                    psamps = l[sindex:]
                elif l[0] == 'Type':
                    stypes = l[sindex:]
                    for si in range(0, len(stypes)):
                        if stypes[si] == 'Tumor':
                            tsinduse.append(si)
                elif l[0] == 'Gender':
                    sexes = l[sindex:]
                elif l[0] == 'Ethnicity':
                    ancestry = l[sindex:]
                elif l[0] == 'mRNA.Expression.Subtype.TCGA':
                    expST = l[sindex:]
                elif l[0] == 'NMF.consensus':
                    nmfc = l[sindex:]
                    
            # Proteomics rows start at row 71 (0-index)
            else:
                GS, mods = l[2], [x for x in l[9].split(' ') if x != '']
                mods = '_'.join(sorted([x[0:-1] for x in mods]))
                pep = GS + '_' + mods  
                peps.append(pep)

                vals0 = l[sindex:]
                vals = []
                for v in vals0:
                    if v == 'NA': v = np.nan
                    else: v = float(v)
                    vals.append(v)
                pmat.append(vals)
              
    pmat = np.array(pmat)[:, tsinduse]
    psamps = np.array(psamps)[tsinduse]
    peps = np.array(peps)
    print(pmat.shape)

    expST = np.array(expST)[tsinduse] 
    expST[expST == 'Proximal-inflammatory'] = 'PI'
    expST[expST == 'Proximal-proliferative'] = 'PP'
    expST[expST == 'Terminal Respiratory Unit'] = 'TRU'
    nmfc = np.array(nmfc)[tsinduse]
    sexes = np.array(sexes)[tsinduse]
    ancestry = np.array(ancestry)[tsinduse]

    # Reduce to valid cohort samples and align
    psinds = []
    for si, s in enumerate(psamps):
        ST = expST[si]
        if ST in ['PI','PP','TRU']: psinds.append(si)
    pmat = pmat[:, psinds]
    psamps = psamps[psinds]
    expST = expST[psinds]
    nmfc = nmfc[psinds]
    sexes = sexes[psinds]
    ancestry = ancestry[psinds]
    print(pmat.shape) 
   
    # Adjust by number of measurements
    if measured_fraction != None:
        pnansumsi = (np.sum(~np.isnan(pmat), axis = 1) >= (measured_fraction * pmat.shape[1]))
        pmat = pmat[pnansumsi, :]
        peps = peps[pnansumsi]
    print(pmat.shape)
    
    # Imputation
    if do_impute:
        imputer = KNNImputer(n_neighbors = kneighbor, weights = 'uniform', metric = 'nan_euclidean')
        pmat = imputer.fit_transform(pmat)
 
    # write out sample to expression subtype, sex, and ancestry files
    ofst = open(outdir + 'sample_expression_subtypes.txt', 'w')
    ofnmf = open(outdir + 'sample_NMF_clusters.txt', 'w')
    ofsex = open(outdir + 'sample_sexes.txt', 'w')
    ofanc = open(outdir + 'sample_ancestries.txt', 'w')
    for (s, st, nmf, sex, anc) in zip(psamps, expST, nmfc, sexes, ancestry):
        ofst.writelines('%s\t%s\n'%(s, st))
        ofnmf.writelines('%s\t%s\n'%(s, nmf))
        ofsex.writelines('%s\t%s\n'%(s, sex))
        ofanc.writelines('%s\t%s\n'%(s, anc))
    ofst.close()
    ofnmf.close()
    ofsex.close()
    ofanc.close()
    
    ## EXECUTE ##
    OUTRES = {} # dictionary to store patient results
    KIN2PEP = {} # dictionary for kinase to peptide matches
    # Loop over samples and find evidence for kinase sub-strates; then perform regressions
    for si, samp in enumerate(psamps):
        pnonnani = ~np.isnan(pmat[:, si]) # get non NaN indexes for MS data
 
        mat_samp = pmat[pnonnani, si].reshape(-1,1)
        pep_samp = peps[pnonnani]
        mat_all = pmat[pnonnani, :] 
        
        # find kinases
        kinases, pep_inds = {}, set()
        for pi, p in enumerate(pep_samp):
            prot, mods = p.split('_')[0], p.split('_')[1:]
            if prot in ps_dict:
                for mod in mods:
                    if mod in ps_dict[prot]:
                        KIN = ps_dict[prot][mod]
                        pep_inds.add(pi)
                        for kin in KIN:
                            if kin not in kinases: kinases[kin] = []
                            if kin not in KIN2PEP: KIN2PEP[kin] = set()
                            kinases[kin].append(pi)
                            KIN2PEP[kin].add(p)
        
        # Restrict to minimum number of peptides linked to kinases
        if min_peptides != None:
            kinases_new = {}
            for kin in kinases:
                if len(kinases[kin]) < min_peptides: continue
                kinases_new[kin] = kinases[kin]
            kinases = kinases_new
           
        # Handle kinases with same behavior in group
        Klist0 = sorted(kinases.keys())
        kgroups = []
        for ki in range(0, len(Klist0)):
            kg = set()
            for kj in range(0, len(Klist0)):
                if ki == kj: continue
                if kinases[Klist0[ki]] == kinases[Klist0[kj]]:
                    kg.add(ki)
                    kg.add(kj)
            kg = sorted(list(kg))
            if len(kg) > 0 and kg not in kgroups:
                kgroups.append(kg)
        
        Klist = [x for x in Klist0]
        for kg in kgroups:
            kcomb = '_'.join(Klist[x] for x in kg)
            Klist[kg[0]] = kcomb
            for xi in kg[1:]: Klist[xi] = 'REMOVER'
        Klist = [x for x in Klist if x != 'REMOVER']

        # create X and Y data for regression
        pep_inds = sorted(list(pep_inds))
        Xmat = np.zeros((len(pep_inds), 1+len(Klist))) # add one for intercept
        Ymat = np.zeros((len(pep_inds), 1))
        pi_map = {}
        for ii, pi in enumerate(pep_inds):
            vmed = np.nanmedian(mat_all[pi, :])
            vmad = np.nanmedian(abs(mat_all[pi, :] - vmed))
            if vmad == 0: vmad = 1 # can occur if n = 1
            if normalization == 'medcenter':
                v = mat_samp[pi, 0] - vmed # median centering
            elif normalization == 'rZscore':
                v = (mat_samp[pi, 0] - vmed) / vmad # Robust z-score 
            Ymat[ii, 0] = v
            pi_map[pi] = ii
        Xmat[:, 0] = 1+Xmat[:,0]
        for ki, kin in enumerate(Klist):
            kin = kin.split('_')[0] # for cases of merged kinases
            pi = kinases[kin]
            pii = [pi_map[x] for x in pi]
            Xmat[pii, ki+1] = 1
        
        # run regression
        BETA, PV = run_regression(Xmat, Ymat, regType, set_alpha = 0.1, nrep = nrep)
        #mod = SM.OLS(Ymat, Xmat).fit()
        #print mod.summary()
        
        # Collect results
        for i in range(1, len(BETA)): # start at 1 to ignore intercept
            KIN = Klist[i-1]
            for kin in KIN.split('_'):
                if kin not in OUTRES: OUTRES[kin] = {}
                OUTRES[kin][samp] = [BETA[i], PV[i]]

        print(samp, si)
        #if si > 5:
        #    break

    # collect data
    KOUT = sorted(OUTRES.keys())
    OUTBETA = np.zeros((len(OUTRES), len(psamps)))
    OUTPV = np.zeros((len(OUTRES), len(psamps)))
    OUTSCORE = np.zeros((len(OUTRES), len(psamps)))
    for ki, kin in enumerate(KOUT):
        for si, samp in enumerate(psamps):
            if samp not in OUTRES[kin]:
                b, p = [0.0, 1.0]
            else: b,p = OUTRES[kin][samp]
            
            OUTBETA[ki, si] = b
            OUTPV[ki, si] = p
            
            if b >= 0:
                score = -np.log10(p) 
            elif b < 0:
                score = (-1) * -np.log10(p)
            
            OUTSCORE[ki, si] = score

    # write out
    ofbeta = open(outdir + 'ks_regression_beta.txt','w')
    ofpv = open(outdir + 'ks_regression_pv.txt','w')
    ofscore = open(outdir + 'ks_regression_score.txt','w')
    ofk2p = open(outdir + 'ks_kinases_to_peptides.txt','w')
    ofbeta.writelines('\t'.join(['KINASE'] + [ps for ps in psamps])+'\n')
    ofpv.writelines('\t'.join(['KINASE'] + [ps for ps in psamps])+'\n')
    ofscore.writelines('\t'.join(['KINASE'] + [ps for ps in psamps])+'\n')
    ofk2p.writelines('\t'.join(['KINASE','NPEPTIDES','PEPTIDES'])+'\n')
    for ki, kin in enumerate(KOUT):
        vb = OUTBETA[ki, :]
        olb = [kin] + ['%0.4f'%(x) for x in vb]
        ofbeta.writelines('\t'.join(olb)+'\n')

        vp = OUTPV[ki, :]
        olp = [kin] + ['%0.4g'%(x) for x in vp]
        ofpv.writelines('\t'.join(olp)+'\n')

        vs = OUTSCORE[ki, :]
        ols = [kin] + ['%0.4f'%(x) for x in vs]
        ofscore.writelines('\t'.join(ols)+'\n')

    for kin in sorted(KIN2PEP.keys()):
        ntarg = len(KIN2PEP[kin])
        targets = ';'.join(sorted(list(KIN2PEP[kin])))
        ofk2p.writelines('%s\t%d\t%s\n'%(kin, ntarg, targets))

    ofbeta.close()
    ofpv.close()
    ofscore.close()
    ofk2p.close()

if __name__ == '__main__': main()


