from __future__ import division
import sys, os, glob, csv
import _pickle
import numpy as np
import argparse
sys.path.insert(0, '') # Path to stattests
from stattests import hypgeo_test as HGT
from stattests import fdr_BH as FDRBH

######################
# Gene set functions #
######################

def geneSet_to_msigdb_id(set_name):
    '''
    Get MSigDB IDs for pathway sets
    '''
    setid = None
    if set_name == 'Hallmark':
        setid = 'h.all'
    elif set_name == 'GO_BP':
        setid = 'c5.go.bp'
    elif set_name == 'GO_MF':
        setid = 'c5.go.mf'
    elif set_name == 'GO_CC':
        setid = 'c5.go.cc'
    elif set_name == 'canonical_pathways':
        setid = 'c2.cp'
    elif set_name == 'Reactome':
        setid = 'c2.cp.reactome'
    elif set_name == 'KEGG':
        setid = 'c2.cp.kegg'
    elif set_name == 'TFtargets':
        setid = 'c3.tft'
    elif set_name == 'oncoSigs':
        setid = 'c6.all'

    return setid

def load_pathway_data(path_dir, set_name, valid_set_names):
    '''
    Function to load pathway set data. If pre-parsed data not found, function looks for *.gmt file
    for gene set as secondary check and creates needed dictionaries from there.
    '''

    if set_name not in valid_set_names:
        print('Valide gene sets:', valid_set_names)
        print(set_name, 'not a current valid gene set. Exiting.')
        sys.exit()

    # Check if processed files are available
    g2p_found, p2g_found = False, False
    if os.path.isdir(path_dir + '/%s' % (set_name)):
        if os.path.isfile(path_dir+'/%s/%s.gene2path.pkl'%(set_name,set_name)): g2p_found = True
        if os.path.isfile(path_dir+'/%s/%s.path2gene.pkl'%(set_name,set_name)): p2g_found = True

    if g2p_found and p2g_found:
        path2gene = _pickle.load(open(path_dir+'/%s/%s.path2gene.pkl'%(set_name,set_name),'rb'))
        gene2path = _Pickle.load(open(path_dir+'/%s/%s.gene2path.pkl'%(set_name,set_name),'rb'))

        return path2gene, gene2path

    # If not found, look for .gmt file and parse from there
    else:
        setid = geneSet_to_msigdb_id(set_name)
        gmt_file = glob.glob(path_dir + '/%s*' % (setid))

        if len(gmt_file) == 0:
            return None, None

        # If file found, parse to dictionaries
        gmt_file = gmt_file[0]
        path2gene, gene2path = {}, {}
        with open(gmt_file) as GMT:
            reader = csv.reader(GMT, delimiter = '\t')
            for ri, row in enumerate(reader):
                pname = row[0]
                genes = row[2:]
                path2gene[pname] = genes
                for gene in genes:
                    if gene not in gene2path: gene2path[gene] = []
                    gene2path[gene].append(pname)

        # Return parsed dictionaries
        return path2gene, gene2path

def load_pathway_data_from_gmt_file(gmt_file):
    '''
    Load gene set from direct gmt file input.
    '''

    # Parse gmt file
    path2gene, gene2path = {}, {}
    with open(gmt_file) as GMT:
        reader = csv.reader(GMT, delimiter = '\t')
        for ri, row in enumerate(reader):
            pname = row[0]
            genes = row[2:]
            path2gene[pname] = genes
            for gene in genes:
                if gene not in gene2path: gene2path[gene] = []
                gene2path[gene].append(pname)

    # Return parsed dictionaries
    return path2gene, gene2path

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

def main():

    ###############################
    # ARGUMENT AND OPTION PARSING #
    ###############################

    usage='python %s <foreground_genes.txt> <background_genes.txt> <gene_set> [options]'%(sys.argv[0].split('/')[-1])
    description='''
    Performs Gene Ontology term enrichment on a set of foreground genes (e.g. differentially expressed genes) against a background set
    (e.g. all genes, all expressed genes, etc.). Foreground should be subset of background.
    '''
    parser = argparse.ArgumentParser(usage=usage,description=description,formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('foreground',metavar='foreground.txt',type=str,
                        help='Input foreground genes file, one gene per line (Required).')
    parser.add_argument('background',metavar='background.txt',type=str,
                        help='Background genes file, one per line (Required).')
    parser.add_argument('gene_set',metavar='gene_set',type=str,
                        help='Name of gene set (Required). Valid gene sets are listed in accompanying "valid_gene_sets.txt" file.')
    parser.add_argument('-o','--outdir',type=str,default='./',dest='outdir',
                        help='Provide output directory for analysis.')
    parser.add_argument('--prefix',type=str,default='pathway_enrichment',dest='prefix',
                        help='Provide prefix for analysis.')
    parser.add_argument('--pathway-data',type=str,dest='path_data',
                        default='datasets/MSigDB/v2023.1.Hs/',
                        help='Provide path for pathway data. Code looks for files: \
                              <path_data>/<geneSet>.gene2path.pkl and <path_data>/<geneSet>.path2gene.pkl or looks for \
                              .gmt files in directory based on gene_set input. Ignored if direct --gmt input is used.')
    parser.add_argument('--gmt',type=str,dest='gmt_file',default = None,
                        help='Provide full path to a .gmt file for use as a gene set. This can be used instead of --pathway-data \
                        and the <gene_set> input argument will be used as the set name.')
    parser.add_argument('-s','--synonyms',type=str,dest='synonyms',default=None,
                        help='Provide a dictionary of synonyms for gene names.')
    parser.add_argument('-t','--threshold',type=float,dest='threshold',default=0.05,
                        help='Provide Benjamini-Hochberg FDR adjusted p-value threshold for output reporting.')

    args = parser.parse_args()
    vargs = vars(args)

    # Parse required
    fg_fn = vargs['foreground']
    if not os.path.isfile(fg_fn): parser.error('Foreground genes file does not exist!')
    bg_fn = vargs['background']
    if not os.path.isfile(bg_fn): parser.error('Background genes file does not exist!')
    gene_set = vargs['gene_set']

    # Options
    outdir = vargs['outdir']
    os.system('mkdir -p %s'%(outdir))
    if not outdir.endswith('/'): outdir += '/'
    prefix = vargs['prefix'] + '_'
    path_data = vargs['path_data']
    gmt_file = vargs['gmt_file']
    synonyms = vargs['synonyms']
    if synonyms != None:
        if not os.path.isfile(synonyms): parser.error('Synonyms file does not exist')
        synonyms = cPickle.load(open(synonyms))
    threshold = vargs['threshold']

    #############
    # EXECUTION #
    #############

    ## Valid gene sets definition (load from accompanying file)
    codeDir = os.path.abspath(os.path.dirname(__file__))
    valid_sets = [x.strip() for x in open(codeDir + '/valid_gene_sets.txt')]

    # Load pathway data
    print('Loading pathway data...')
    if gmt_file != None:
        if not os.path.isfile(gmt_file): parser.error('Provided gmt file path does not exist!')
        path2gene, gene2path = load_pathway_data_from_gmt_file(gmt_file)
    else:
        path2gene, gene2path = load_pathway_data(path_data, gene_set, valid_sets)
    print('Pathway data loaded.')

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
    print('Computing pathway counts...')
    tot_terms = set()
    nfg_not_found,nbg_not_found = 0,0
    syn_convert = 0
    fg_counts,bg_counts = {},{}
    fg_genes = {}
    for g in background:
        g_used = g
        try:
            terms = gene2path[g]
        except KeyError:
            if synonyms != None:
                if g in synonyms:
                    g_used = synonyms[g]
                    if g_used in gene2path:
                        syn_convert += 1
                        terms = gene2path[g_used]
                    else:
                        nbg_not_found += 1
                        if g in foreground: nfg_not_found += 1
                        continue
                else:
                    nbg_not_found += 1
                    if g in foreground: nfg_not_found += 1
                    continue
            else:
                nbg_not_found += 1
                if g in foreground: nfg_not_found += 1
                continue

        for t in terms:
            tot_terms.add(t)
            if t not in bg_counts: bg_counts[t] = 0
            bg_counts[t] += 1
            if g in foreground:
                if t not in fg_counts: fg_counts[t] = 0
                if t not in fg_genes: fg_genes[t] = []
                fg_counts[t] += 1
                fg_genes[t].append(g)

    print('Term counting finished.')
    print('  %d genes converted by synonyms.'%(syn_convert))
    print('  %d pathways associated with genes.'%(len(tot_terms)))
    print('  %d of %d foreground genes not in database.'%(nfg_not_found,len(foreground)))
    print('  %d of %d background genes not in database.'%(nbg_not_found,len(background)))

    # Compute enrichment statistics for terms
    terms = []
    for t in sorted(list(tot_terms)):
        if t not in fg_counts: continue
        terms.append(t)
    print('  %d pathways associated to at least one foreground gene used for analysis.'%(len(terms)))
    print('Computing hypergeometric tests...')
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

    print('Hypergeometric tests complete.')
    print('Correcting p-values with B-H procedure...')
    olines = sorted(olines,key=lambda x:x[-1])
    pvals = []
    for ol in olines: pvals.append(ol[-1])
    #fdrs = fdr_BH(np.array(pvals))
    fdrs = FDRBH(np.array(pvals))
    print('FDR p-values computed.')

    # Write output
    print('Writing output...')
    ofn = outdir + prefix + '%s_FDR_lt_%s.txt'%(gene_set,str(threshold))
    of = open(ofn,'w')
    of.writelines('\t'.join(['pathway','num_fg','num_bg','tot_fg','tot_bg','enrichment','pvalue','FDR','genes'])+'\n')
    for i,ol in enumerate(olines):
        if fdrs[i] < threshold:
            ol.append(fdrs[i])
            ol.append(','.join(sorted(fg_genes[ol[0]])))
            oll = '%s\t%d\t%d\t%d\t%d\t%0.3f\t%0.3g\t%0.3g\t%s\n'%(tuple(ol))
            of.writelines(oll)
    of.close()

    print('Done.')

if __name__ == '__main__': main()

