import sys, os
import numpy as np

########
# PREP #
########

args = sys.argv[1:]
if len(args) < 2:
    print 'Usage: python <prog> <DE sites file> <output file name> [metric]'
    print "Metric options: 'beta' (default)"
    sys.exit()
inf = args[0]
ofn = args[1]

# defaults
metric = 'beta'
if len(args) > 2:
    metric = args[2]

#######
# RUN #
#######

# Read file, split multi-sites, and deal with ties/multiply quantified sites
PSITES = {}
psi, meti = None, None
if metric == 'beta':
    metval = 'beta'
elif metric == 'signlog10pval':
    metval = 'pval'
    betval = 'beta'
for i,line in enumerate(open(inf).readlines()):
    l = line.strip().split('\t')
    if i == 0:
        for li, lv in enumerate(l):
            if lv == 'ProteinSite': psi = li
            elif lv == metval: meti = li
            if metric == 'signlog10pval':
                if lv == 'beta': betai = li  
        continue

    if psi == None:
        print 'ProteinSite not found in input file. Exiting.'
        sys.exit()
    if meti == None:
        print 'Metric ID for metric %s not found in input file. Exiting.'%(metric)
        sys.exit()

    psite = l[psi]
    if metric == 'signlog10pval':
        pv = float(l[meti])
        if pv == 0: pv = 1e-50
        l10pv = -np.log10(pv)
        bsign = np.sign(float(l[betai]))
        val = bsign * l10pv
    else:
        val = float(l[meti])
    prot, sites = psite.split('_')[0], psite.split('_')[1:]
    for site in sites:
        ps = prot + '_' + site
        if ps not in PSITES: PSITES[ps] = []
        PSITES[ps].append(val)

# for multiples, take absolute maximum value
PSL = []
for ps in PSITES:
    if len(PSITES[ps])>1:
        vals = PSITES[ps]
        ui = np.argwhere(np.abs(np.array(vals))==np.max(np.abs(vals)))[0][0]
        val = vals[ui]
    else:
        val = PSITES[ps][0]

    PSL.append([ps, val])

# sort by metric in descending order
PSL.sort(key = lambda x: x[1], reverse = True)

# write output
of = open(ofn, 'w')
for ps in PSL:
    ol = '%s\t%f\n'%(tuple(ps))
    of.writelines(ol)
of.close()

