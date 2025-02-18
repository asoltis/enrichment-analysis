import sys, os
import numpy as np

########
# PREP #
########

args = sys.argv[1:]
if len(args) < 2:
    print 'Usage: python <prog> <DE MS file> <output file base> [FDR threshold] [RPPA multi-map flag]'
    print "Output files are <output_file base>_bg.txt and <output_file_base>_<direction>_fg.txt, where <direction> is 'all', 'up', or 'down'."
    print "FDR threshold used to consider significant changes vs. non."
    print "Set RPPA multi-map flag to True to consider multi-mappers in analysis; set to False to skip; default is True"
    sys.exit()
inMS = args[0]
ofbase = args[1]

# defaults
fdrT = 0.05
if len(args) > 2:
    fdrT = float(args[2])

#######
# RUN #
#######

# Read file, split multi-sites, and deal with ties/multiply quantified sites
PSITES = {}
metric = 'beta'

# MS
psi, fdri, meti = None, None, None
for i,line in enumerate(open(inMS).readlines()):
    l = line.strip().split('\t')
    if i == 0:
        for li, lv in enumerate(l):
            if lv == 'ProteinSite': psi = li
            elif lv == 'FDR': fdri = li
            elif lv == metric: meti = li
        continue

    if psi == None:
        print 'ProteinSite not found in input MS file. Exiting.'
        sys.exit()
    if fdri == None:
        print 'FDR column not found in input MS file. Exiting.'
        sys.exit()
    if meti == None:
        print 'Metric field for %s not found in MS file. Exiting.'%(metric)
        sys.exit()

    psite, val = l[psi], float(l[fdri])
    met = float(l[meti])
    prot, sites = psite.split('_')[0], psite.split('_')[1:]
    for site in sites:
        ps = prot + '_' + site
        if ps not in PSITES: PSITES[ps] = []
        PSITES[ps].append([met, val])

# Get foreground and background lists
# For up/down, take most significant for multiples
bgset = set()
fgset, fgUpset, fgDownset = set(), set(), set()
for ps in PSITES:
    psinfo = PSITES[ps]
    fdrs = np.array([x[1] for x in psinfo])
    mi = np.argwhere(fdrs == min(fdrs))[0][0]
    mmet, mfdr = psinfo[mi]

    if mfdr < fdrT:
        fgset.add(ps)
        if mmet < 0: fgDownset.add(ps)
        elif mmet > 0: fgUpset.add(ps)
    bgset.add(ps)

# write output
offgA = open(ofbase + '_all_fg.txt','w')
offgU = open(ofbase + '_up_fg.txt','w')
offgD = open(ofbase + '_down_fg.txt','w')
ofbg = open(ofbase + '_bg.txt','w')
for ps in sorted(list(fgset)): offgA.writelines('%s\n'%(ps))
for ps in sorted(list(fgUpset)): offgU.writelines('%s\n'%(ps))
for ps in sorted(list(fgDownset)): offgD.writelines('%s\n'%(ps))
for ps in sorted(list(bgset)): ofbg.writelines('%s\n'%(ps))
offgA.close()
offgU.close()
offgD.close()
ofbg.close()
