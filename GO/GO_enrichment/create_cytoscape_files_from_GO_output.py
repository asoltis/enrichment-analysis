from __future__ import division
import sys, os
import numpy as np
import argparse

def overlap(set1,set2,method):
    '''
    '''
    olap = None

    if method == 'jaccard':
        intersect = [x for x in set1 if x in set2]
        union = list(set(set1+set2))
        olap = len(intersect) / len(union)

    elif method == 'coefficient':
        intersect = [x for x in set1 if x in set2]
        denom = min(len(set1),len(set2))
        olap = len(intersect) / denom

    return olap

def main():
    ###############################
    # ARGUMENT AND OPTION PARSING #
    ###############################

    usage='python %s <GO_enrichment_output_file> [options]'%(sys.argv[0].split('/')[-1])
    description='''
    Creates ...
    '''
    parser = argparse.ArgumentParser(usage=usage,description=description,formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('GOfn',metavar='GOfn.txt',type=str,
                        help='Input GO enrichment output file (Required).')

    parser.add_argument('-o','--outdir',type=str,default='./',dest='outdir',
                        help='Provide output directory for analysis.')
    parser.add_argument('--prefix',type=str,default='GO_enrichment',dest='prefix',
                        help='Provide prefix for analysis.')
    parser.add_argument('-q','--qvalue-thresh',type=float,default=0.1,dest='qvalT',
                        help='Select q-value threshold for included GO terms.')
    parser.add_argument('--olap-stat',type=str,default='jaccard',dest='method',
                        help="Select overlap statistic; valid choices are: \
                              1) 'jaccard' - intersection over union; \
                              2) 'coefficient' - intersection over min(|A|,|B|)")
    parser.add_argument('--olap-thresh',type=float,default=0.25,dest='olapT',
                        help='Select overlap index cut-off threshold creating edges in output file.')
    parser.add_argument('--sig-metric',type=str,default='neg_log10_qval',dest='sig_metric',
                        help='Select significance metric for output. Valid choices are: neg_log10_qval, enrichment.')

    args = parser.parse_args()
    vargs = vars(args)

    # parse values
    GOfn = vargs['GOfn']

    # parse options
    outdir = vargs['outdir']
    os.system('mkdir -p %s'%(outdir))
    if not outdir.endswith('/'): outdir += '/'
    prefix = vargs['prefix'] + '_'
    qvalT = vargs['qvalT']
    method = vargs['method']
    if method not in ['jaccard','coefficient']:
        parser.error('Invalid overlap statistic selected.')
    olapT = vargs['olapT']
    sig_metric = vargs['sig_metric']
    if sig_metric not in ['neg_log10_qval','enrichment']:
        parser.error('Invalid significance metric')

    # Load in GO output file info
    GO_data = {}
    for line in open(GOfn).readlines()[1:]:
        l = line.strip().split('\t')
        term,name,count,enrich,qval,genes = l[0],l[1],int(l[3]),float(l[7]),float(l[9]),l[10]
        genes = genes.split(',')

        if qval < qvalT:

            if sig_metric == 'neg_log10_qval': sig = -np.log10(qval)
            elif sig_metric == 'enrichment': sig = enrich

            GO_data[term] = {}
            GO_data[term]['name'] = name
            GO_data[term]['count'] = count
            GO_data[term]['enrichment'] = enrich
            GO_data[term]['significance'] = sig
            GO_data[term]['genes'] = genes

        else: continue

    # Create graph from data
    graph = []
    terms = list(GO_data.keys())
    for i in range(0,len(terms)-1):
        for j in range(i+1,len(terms)):
            t1,t2 = terms[i],terms[j]
            s1,s2 = GO_data[t1]['genes'],GO_data[t2]['genes']
            olap = overlap(s1,s2,method)

            if olap >= olapT:
                graph.append((t1,t2,olap))

    # Write outputs
    gof = open(outdir+prefix+'graph.sif','w')
    eaf = open(outdir+prefix+'edgeAttributes.txt','w')
    eaf.writelines('Edge\tWeight\n')
    for e in graph:
        gof.writelines('%s\tpp\t%s\n'%(e[0],e[1]))
        eaf.writelines('%s (pp) %s\t%0.3f\n'%(e))
    gof.close()
    eaf.close()

    naf = open(outdir+prefix+'nodeAttributes.txt','w')
    naf.writelines('Node\tTermName\tSet_size\tSignificance\n')
    for node in GO_data:
        name = GO_data[node]['name']
        count = GO_data[node]['count']
        sig = GO_data[node]['significance']
        naf.writelines('%s\t%s\t%d\t%0.3f\n'%(node,name,count,sig))
    naf.close()

if __name__ == '__main__': main()

