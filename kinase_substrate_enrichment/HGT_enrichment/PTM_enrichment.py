from __future__ import division
import sys,os
import cPickle
import numpy as np
import argparse
from stattests import hypgeo_test as HGT

def main():

    ###############################
    # ARGUMENT AND OPTION PARSING #
    ###############################
    
    usage='python %s <foreground_sites.txt> <background_sites.txt> [options]'%(sys.argv[0].split('/')[-1])
    description='''Post-translational modification (PTM) enrichment analysis.'''

    parser = argparse.ArgumentParser(usage=usage,description=description,formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('foreground',metavar='foreground.txt',type=str,
                        help='Input foreground protein sites file, one per line (Required).')
    parser.add_argument('background',metavar='background.txt',type=str,
                        help='Background protein sites file, one per line (Required).')
    parser.add_argument('-o','--outdir',type=str,default='./',dest='outdir',
                        help='Provide output directory for analysis.')
    parser.add_argument('--prefix',type=str,default='ptm_enrichment',dest='prefix',
                        help='Provide prefix for analysis.')
    parser.add_argument('--db-data',type=str,dest='db_data',
                        default='PTMsigDB/processed/db_files/uniprot_human_v1.9.0/PSP-KS',
                        help='Provide path and prefix to database, e.g.: \
                              <db_data>.pep2path.pkl and <db_data>.path2pep.pkl') 
    parser.add_argument('-t','--threshold',type=float,dest='threshold',default=0.05,
                        help='Provide Benjamini-Hochberg FDR adjusted p-value threshold for output reporting.')
    
    args = parser.parse_args()
    vargs = vars(args)

    # Parse required
    fg_fn = vargs['foreground']
    if not os.path.isfile(fg_fn): parser.error('Foreground sites file does not exist!')
    bg_fn = vargs['background']
    if not os.path.isfile(bg_fn): parser.error('Background sites file does not exist!')
    
    # Options
    outdir = vargs['outdir']
    os.system('mkdir -p %s'%(outdir))
    if not outdir.endswith('/'): outdir += '/'
    prefix = vargs['prefix'] + '_'
    db_data = vargs['db_data']
    if not os.path.isfile(db_data+'.pep2path.pkl'): parser.error('Gene to path info dictionary not found!')
    if not os.path.isfile(db_data+'.path2pep.pkl'): parser.error('Path to gene dictionary not found!') 
    threshold = vargs['threshold']
    
    #############
    # EXECUTION #
    #############

    # Load dictionaries
    print 'Loading dictionaries...'
    path2pep = cPickle.load(open(db_data+'.path2pep.pkl'))
    pep2path = cPickle.load(open(db_data+'.pep2path.pkl'))
    print 'Dictionaries loaded.\n'

    # Read foreground and background files
    background = [x.strip() for x in open(bg_fn).readlines()]
    foreground = []
    for l in open(fg_fn).readlines():
        g = l.strip()
        if g not in background:
            print ' %s not in background gene set. Skipping.'%(g)
        else: foreground.append(g)
    # convert to sets to eliminate any duplicates
    foreground = set(foreground)
    background = set(background) 
   
    # check format of input sites
    site_len = set([len(x.split('_')) for x in background]) 
    if len(site_len) > 1:
        print 'Variable format input data. Exiting.'
        sys.exit()
    site_len = list(site_len)[0]
    if site_len == 2:
        print 'Modification information and direction not provided. Assuming sites are phospho ("p") and ignoring direction.'
    if site_len == 3:
        print 'Direction not provided. Ignoring direction information.'

    # Count terms from foreground and background
    print 'Computing PTM database counts...'
    tot_terms = set()
    nfg_not_found,nbg_not_found = 0,0
    syn_convert = 0
    fg_counts,bg_counts = {},{}
    fg_sites = {}
    for g in background:
        g_used = g
        FOUND = False
        terms = [] 
        
        if g_used in pep2path: terms += pep2path[g_used]
        if len(terms) > 0: FOUND = True

        if not FOUND:
            nbg_not_found += 1
            if g in foreground: nfg_not_found += 1
            continue

        for t in terms:
            tot_terms.add(t)
            if t not in bg_counts: bg_counts[t] = 0
            bg_counts[t] += 1
            if g in foreground:
                if t not in fg_counts: fg_counts[t] = 0
                if t not in fg_sites: fg_sites[t] = []
                fg_counts[t] += 1
                fg_sites[t].append(g)

    print 'Term counting finished.'
    print '  %d values converted by synonyms.'%(syn_convert)
    print '  %d PTM database sets associated with sites.'%(len(tot_terms))
    print '  %d of %d foreground sites not in database.'%(nfg_not_found,len(foreground))
    print '  %d of %d background sites not in database.'%(nbg_not_found,len(background))
    
    # Compute enrichment statistics for terms
    terms = []
    for t in sorted(list(tot_terms)):
        if t not in fg_counts: continue
        terms.append(t)
    print '  %d PTM databases sets associated to at least one foreground PTM site used for analysis.'%(len(terms))
    print 'Computing hypergeometric tests...'
    olines = []
    for t in terms:
        K = bg_counts[t]
        x = fg_counts[t]
        M = len(background) #- nbg_not_found
        N = len(foreground) #- nfg_not_found

        [d1,pv,d2] = HGT(x,M,K,N) # take right tail p-value
        enrich = (x / N) / (K / M)
        ol = [t,x,K,N,M,enrich,pv]
        olines.append(ol)

    print 'Hypergeometric tests complete.'
    print 'Correcting p-values with B-H procedure...'
    olines = sorted(olines,key=lambda x:x[-1])
    pvals = []
    for ol in olines: pvals.append(ol[-1])
    fdrs = fdr_BH(np.array(pvals))
    print 'FDR p-values computed.'

    # Write output 
    ofn = outdir + prefix + 'FDR_lt_%s.txt'%(str(threshold))
    of = open(ofn,'w')
    of.writelines('\t'.join(['PTM_set','num_fg','num_bg','tot_fg','tot_bg','enrichment','pvalue','FDR','sites'])+'\n')
    for i,ol in enumerate(olines):
        if fdrs[i] < threshold:
            ol.append(fdrs[i])
            ol.append(','.join(sorted(fg_sites[ol[0]])))
            oll = '%s\t%d\t%d\t%d\t%d\t%0.3f\t%0.3g\t%0.3g\t%s\n'%(tuple(ol))
            of.writelines(oll)
    of.close()    

##################
# MATH FUNCTIONS #
##################

def fdr_BH(pvals):
    '''
    Compute Benjamini-Hochberg FDR q-values on sorted array of p-values.
    '''
    total = pvals.size
    fdrs0 = total*pvals/range(1,total+1)
    fdrs = []
    # Preserve monotonicity
    for i in range(0,total):
        fdrs.append(min(fdrs0[i:]))
    return fdrs


if __name__ == '__main__': main()

