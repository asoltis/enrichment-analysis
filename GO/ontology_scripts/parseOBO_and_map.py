import sys, os, gzip
#from OBOparser import OBOparser
from urllib.request import urlopen
import networkx as nx
import obonet
import _pickle

### Inputs for analysis
# URLs
OBOurl = "http://purl.obolibrary.org/obo/go/go-basic.obo"
GAFurl = "http://geneontology.org/gene-associations/goa_human.gaf.gz"

# File names
date = '20230727' # convention is to use GO database date
OBOfn = 'OBO_files/go-basic-%s.obo' % (date)
gaf_fn = 'human/%s/goa_human.gaf' % (date)
obase = 'human/%s/goa_human_%s_' % (date, date)

# Directories
os.makedirs('OBO_files', exist_ok = True)
os.makedirs('human/%s/' % (date), exist_ok = True)

# Get OBO file from url and write file
print('Obtaining OBO file...')
oboFH = urlopen(OBOurl)
OBO = open(OBOfn, 'w')
for line in oboFH.readlines():
    OBO.writelines(line.decode())
OBO.close()
print('OBO file obtained.')

# Download GAF human file and write to file
print('Obtaining GAF file...')
gafFH = urlopen(GAFurl)
gaf_lines = None
with gzip.GzipFile(fileobj = gafFH) as GF:
    gaf_lines = GF.read()
with open(gaf_fn, 'wb') as GOF:
    GOF.write(gaf_lines)
GOF.close()
print('GAF file obtained.')

# Create dictionary for GO terms (names, ancestors)
print('Creating GO dictionary.')
GO_dict = {}
goid, goname, gons = None, None, None
for line in open(OBOfn):
    l = line.strip()
    if l == '[Term]':
        if goid != None:
            GO_dict[goid] = {}
            GO_dict[goid]['name'] = goname
            GO_dict[goid]['namespace'] = gons
    elif l.startswith('id: '):
        goid0 = l.split('id: ')[1].strip()
        if not goid0.startswith('GO:'): break
        goid = goid0
    elif l.startswith('name: '):
        goname = l.split('name: ')[1].strip()
    elif l.startswith('namespace: '):
        gons = l.split('namespace: ')[1].strip()
if goid not in GO_dict: # Handle last entry
    GO_dict[goid] = {}
    GO_dict[goid]['name'] = goname
    GO_dict[goid]['namespace'] = gons

# Create GO graph with obonet
print('Loading OBO graph...')
DAG = obonet.read_obo(OBOfn).reverse() # have to reverse edges here to children GO terms point back to ancestors
print('Graph loaded')
print('Getting ancestors...')
for term in DAG.nodes:
    anc = nx.ancestors(DAG, term)
    GO_dict[term]['ancestors'] = anc
print('Ancestors processed.')

# Now read gene associations file
# Create GO2Gene and Gene2GO dictionaries
print('Creating dictionaries for analysis...')
Gene2GO,GO2Gene = {},{}
for line in open(gaf_fn).readlines():
    if line.startswith('!'): continue # header
    l = line.strip().split('\t')
    gname, GO = l[2],l[4]

    if GO not in DAG.nodes:
        print('  Term %s not found. Likely obsolete.'%(GO))
        continue

    ancestors = GO_dict[GO]['ancestors']

    if gname not in Gene2GO: Gene2GO[gname] = set()
    Gene2GO[gname].add(GO)
    for anc in ancestors: Gene2GO[gname].add(anc)

    if GO not in GO2Gene: GO2Gene[GO] = set()
    GO2Gene[GO].add(gname)
    for anc in ancestors:
        if anc not in GO2Gene: GO2Gene[anc] = set()
        GO2Gene[anc].add(gname)

# save data
print('Writing output...')
# GO term info
_pickle.dump(GO_dict,open(obase+'GO_term_info_dictionary.pkl','wb'))
# DAG
_pickle.dump(DAG,open(obase+'DAG.networkx.pkl','wb'))
# GO2Gene
_pickle.dump(GO2Gene,open(obase+'GO2Gene_dictionary.pkl','wb'))
# Gene2GO
_pickle.dump(Gene2GO,open(obase+'Gene2GO_dictionary.pkl','wb'))

print('Done.')

