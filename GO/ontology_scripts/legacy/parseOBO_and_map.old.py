import sys,os
from OBOparser import OBOparser
import networkx as nx
import _pickle

# Inputs for analysis
url = "http://purl.obolibrary.org/obo/go/go-basic.obo"
date = '20230727'
OBOfn = 'OBO_files/go-basic-%s.obo' % (date)
ga_fn = 'human/%s/goa_human.gaf' % (date)
obase = 'human/%s/goa_human_%s_' % (date, date)
os.makedirs('OBO_files', exist_ok = True)
os.makedirs('human/%s/' % (date), exist_ok = True)

# Initialize parser
parser = OBOparser()
parser.setUrl(url)

# Read OBO, write to output file
OBOlines = parser.readFile()
OBO = open(OBOfn,'w')
for line in OBOlines: OBO.writelines(line.decode())

# Create the ontology
ontology = parser.createOntologyFromOBOFile(parser.oboFile)
#ontology = parser.createOntologyFromOBOFile(open(OBOfn).readlines()) # can use with pre-downloaded OBO or if URL get is not working
terms = ontology.terms

# Create Digraph and dictionary of term info
DAG = nx.DiGraph()
GO_dict = {}
alt2GO = {} # Dictionary for alt IDs to GO IDs
for i,term in enumerate(terms):
    GOid = term.getId()
    if GOid == '': continue
    if term.isObsolete == 'true': continue

    GO_dict[GOid] = {}
    GO_dict[GOid]['name'] = term.getName()
    GO_dict[GOid]['namespace'] = term.getNamespace()
    altId = term.getAlternativeId()
    if altId != '': alt2GO[altId] = GOid

    isA = term.getIsA()
    relation = term.getRelationship()

    # Add edges for isA
    for ancest in isA:
        DAG.add_edge(ancest,GOid)

    # Add edges for part of
    if 'part_of' in relation:
        rel = relation.split(' ')[1]
        DAG.add_edge(rel,GOid)

    if (i+1) % 5000 == 0:
        print('  processed %d terms of %d'%(i+1,len(terms)))

# Get ancestors for GO terms
for GO in GO_dict:
    anc = nx.ancestors(DAG,GO)
    GO_dict[GO]['ancestors'] = anc

# Now read gene associations file
# Create GO2Gene and Gene2GO dictionaries
Gene2GO,GO2Gene = {},{}
for line in open(ga_fn).readlines():
    if line.startswith('!'): continue # header
    l = line.strip().split('\t')
    gname,GO = l[2],l[4]

    if GO not in GO_dict:
        if GO in alt2GO:
            print('  Using altID %s for %s'%(alt2GO[GO],GO))
            GO = alt2GO[GO]
        else:
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
# GO term info
_pickle.dump(GO_dict,open(obase+'GO_term_info_dictionary.pkl','wb'))
# DAG
_pickle.dump(DAG,open(obase+'DAG.networkx.pkl','wb'))
# GO2Gene
_pickle.dump(GO2Gene,open(obase+'GO2Gene_dictionary.pkl','wb'))
# Gene2GO
_pickle.dump(Gene2GO,open(obase+'Gene2GO_dictionary.pkl','wb'))

