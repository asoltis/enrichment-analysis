from __future__ import division
import sys,os
import _pickle
import numpy as np
import argparse
sys.path.insert(0, '') # Path to stattests code
from stattests import hypgeo_test as HGT

'''
Updates:

- 2019-01-02
    - Included --terms option to allow user to supply list of GO terms to which analysis should be restricted.
'''

def main():

    ###############################
    # ARGUMENT AND OPTION PARSING #
    ###############################

    usage='python %s <foreground_genes.txt> <background_genes.txt> [options]'%(sys.argv[0].split('/')[-1])
    description='''
    Performs Gene Ontology term enrichment on a set of foreground genes (e.g. differentially expressed genes) against a background set
    (e.g. all genes, all expressed genes, etc.). Foreground should be subset of background.
    '''
    parser = argparse.ArgumentParser(usage=usage,description=description,formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('foreground',metavar='foreground.txt',type=str,
                        help='Input foreground genes file, one gene per line (Required).')
    parser.add_argument('background',metavar='background.txt',type=str,
                        help='Background genes file, one per line (Required).')

    parser.add_argument('-o','--outdir',type=str,default='./',dest='outdir',
                        help='Provide output directory for analysis.')
    parser.add_argument('--prefix',type=str,default='GO_enrichment',dest='prefix',
                        help='Provide prefix for analysis.')
    parser.add_argument('--GO-data',type=str,dest='GO_data',
                        default='GeneOntology/human/20230727/goa_human_20230727',
                        help='Provide path and prefix to gene ontology mapping data. Code will look for files for analysis: \
                              1) <GO_data>_GO_term_info_dictionary.pkl \
                              2) <GO_data>_GO2Gene_dictionary.pkl \
                              3) <GO_data>_Gene2GO_dictionary.pkl')
    parser.add_argument('-s','--synonyms',type=str,dest='synonyms',default=None,
                        help='Provide a dictionary of synonyms for gene names.')
    parser.add_argument('-t','--threshold',type=float,dest='threshold',default=0.05,
                        help='Provide Benjamini-Hochberg FDR adjusted p-value threshold for output reporting.')
    parser.add_argument('--ontology',type=str,default=None,dest='ontology',
                        help='Restrict analysis to only one of biological_process, molecular_function, or cellular_component.')
    parser.add_argument('--restrict',type=str,default=None,dest='restrict',
                        help='Provide filters on number of genes associated with GO term as <min>,<max> (comma-separated list). \
                             This tells program to only consider GO terms with > <min> genes and less than <max> genes.')
    parser.add_argument('--terms',type=str,default = None, dest = 'terms_fn',
                        help='Provide a list of GO term IDs on which analysis should be run. Results from these specific terms will \
                        be reported only.')

    args = parser.parse_args()
    vargs = vars(args)

    # Parse required
    fg_fn = vargs['foreground']
    if not os.path.isfile(fg_fn): parser.error('Foreground genes file does not exist!')
    bg_fn = vargs['background']
    if not os.path.isfile(bg_fn): parser.error('Background genes file does not exist!')

    # Options
    outdir = vargs['outdir']
    os.system('mkdir -p %s'%(outdir))
    if not outdir.endswith('/'): outdir += '/'
    prefix = vargs['prefix'] + '_'
    GO_data = vargs['GO_data']
    if not os.path.isfile(GO_data+'_GO_term_info_dictionary.pkl'): parser.error('GO term info dictionary not found!')
    if not os.path.isfile(GO_data+'_GO2Gene_dictionary.pkl'): parser.error('GO2Gene dictionary not found!')
    if not os.path.isfile(GO_data+'_Gene2GO_dictionary.pkl'): parser.error('Gene2GO dictionary not found!')
    synonyms = vargs['synonyms']
    if synonyms != None:
        if not os.path.isfile(synonyms): parser.error('Synonyms file does not exist')
        synonyms = _pickle.load(open(synonyms, 'rb'))
    threshold = vargs['threshold']
    restrict = vargs['restrict']
    if restrict != None: restrict = [int(x) for x in restrict.split(',')]
    terms_fn = vargs['terms_fn']
    if terms_fn != None:
        terms_to_use = []
        for line in open(terms_fn).readlines(): terms_to_use.append(line.strip())
    ontology = vargs['ontology']
    if ontology not in ['biological_process','molecular_function','cellular_component']:
        parser.error('Provided ontology is not valid.')
    #nprocessors = min(mp.cpu_count(),vargs['nprocessors'])

    #############
    # EXECUTION #
    #############

    # Load dictionaries
    print('Loading dictionaries...')
    GO_info = _pickle.load(open(GO_data+'_GO_term_info_dictionary.pkl', 'rb'))
    GO2Gene = _pickle.load(open(GO_data+'_GO2Gene_dictionary.pkl', 'rb'))
    Gene2GO = _pickle.load(open(GO_data+'_Gene2GO_dictionary.pkl', 'rb'))
    print('Dictionaries loaded.')

    # Read foreground and background files
    background = [x.strip() for x in open(bg_fn).readlines()]
    foreground = []
    for l in open(fg_fn).readlines():
        g = l.strip()
        if g not in background:
            print(' %s not in background gene set. Skipping.'%(g))
        else: foreground.append(g)
    # convert to sets to eliminate any duplicates
    foreground = set(foreground)
    background = set(background)

    # Count terms from foreground and background
    print('Computing GO term counts...')
    tot_terms = set()
    nfg_not_found,nbg_not_found = 0,0
    syn_convert = 0
    fg_counts,bg_counts = {},{}
    fg_genes = {}
    for g in background:
        g_used = g
        try:
            terms = Gene2GO[g]
        except KeyError:
            if synonyms != None:
                if g in synonyms:
                    g_used = synonyms[g]
                    if g_used in Gene2GO:
                        syn_convert += 1
                        terms = Gene2GO[g_used]
                    else:
                        nbg_not_found += 1
                        if g in foreground: nfg_not_found += 1
                        continue
                        #print(' gene %s synonyms %s not in Gene2GO dictionary. Skipping'%(g,g_used)
                else:
                    nbg_not_found += 1
                    if g in foreground: nfg_not_found += 1
                    continue
                    #print(' gene %s not in dictionaries or supplied synonyms dictionary. Skipping'%(g)
            else:
                #print(' gene %s not found in dictionaries. Skipping.'%(g)
                nbg_not_found += 1
                if g in foreground: nfg_not_found += 1
                continue

        for t in terms:
            if ontology != None:
                if GO_info[t]['namespace'] != ontology: continue
            if restrict != None:
                nGenes = len(GO2Gene[t])
                if (nGenes < restrict[0]) or (nGenes > restrict[1]): continue
            if terms_fn != None:
                if t not in terms_to_use: continue
            tot_terms.add(t)
            if t not in bg_counts: bg_counts[t] = 0
            bg_counts[t] += 1
            if g in foreground:
                if t not in fg_counts: fg_counts[t] = 0
                if t not in fg_genes: fg_genes[t] = []
                fg_counts[t] += 1
                fg_genes[t].append(g)

    print('Term counting finished.')
    print('  %d genes converted by synonym.'%(syn_convert))
    print('  %d GO terms associated with genes.'%(len(tot_terms)))
    print('  %d of %d foreground genes not found.'%(nfg_not_found,len(foreground)))
    print('  %d of %d background genes not found.'%(nbg_not_found,len(background)))

    # Compute enrichment statistics for terms
    terms = []
    for t in sorted(list(tot_terms)):
        if t not in fg_counts: continue
        terms.append(t)
    print('  %d terms associated to at least one foreground gene used for analysis.'%(len(terms)))
    print('Computing hypergeometric tests...')
    olines = []
    for t in terms:
        K = bg_counts[t]
        x = fg_counts[t]
        M = len(background) #- nbg_not_found
        N = len(foreground) #- nfg_not_found

        [d1,pv,d2] = HGT(x,M,K,N) # take right tail p-value
        enrich = (x / N) / (K / M)
        ol = [t,GO_info[t]['name'],GO_info[t]['namespace'],x,K,N,M,enrich,pv]
        olines.append(ol)

    print('Hypergeometric tests complete.')
    print('Correcting p-values with B-H procedure...')
    olines = sorted(olines,key=lambda x:x[-1])
    pvals = []
    for ol in olines: pvals.append(ol[-1])
    fdrs = fdr_BH(np.array(pvals))
    print('FDR p-values computed.')

    # Write output
    ofn = outdir + prefix + 'GO_FDR_lt_%s'%(str(threshold))
    if ontology != None:
        ofn += '_%s' %(ontology)
    if restrict != None:
        ofn += '_restrict_%d_%d'%(restrict[0],restrict[1])
    ofn += '.txt'
    of = open(ofn,'w')
    of.writelines('\t'.join(['GO_term','name','ontology','num_fg','num_bg','tot_fg','tot_bg','enrichment','pvalue','FDR','genes'])+'\n')
    for i,ol in enumerate(olines):
        if fdrs[i] < threshold:
            ol.append(fdrs[i])
            ol.append(','.join(sorted(fg_genes[ol[0]])))
            oll = '%s\t%s\t%s\t%d\t%d\t%d\t%d\t%0.3f\t%0.3g\t%0.3g\t%s\n'%(tuple(ol))
            of.writelines(oll)
    of.close()

    print('Done.')

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

