# Code and K-S database paths
gsea=GSEA/GSEA_Linux_4.0.3/gsea-cli.sh # path to gsea client shell script (command line tool)
prep=prep_MS_for_GSEA.py # prep script for input to GSEA tool
db= # GMT file with kinase-substrate links

# input file
inMS= # input text file with DE proteins (including stats)
ST= # ID for subtype/group in question

# ranking metric
metric=signlog10pval

# output directories
ODIR=results_rank_${metric}/$ST
mkdir -p $ODIR

# input file info
files=$ODIR/input_files.txt
> $files
echo $inMS >> $files
echo $db >> $files

# rank file name
rnk=$ODIR/${ST}_MS_phos.rnk

# run prep
python $prep $inMS $rnk $metric

# params
scheme=weighted # weighted is for p=1 weighting
norm=meandiv # meandiv is NES
nperm=100000
label=$ST
seed=1
minset=3

$gsea GSEAPreranked -rnk $rnk -gmx $db -scoring_scheme $scheme -norm $norm -nperm $nperm -rpt_label $label -rnd_seed $seed -out $ODIR -set_min $minset



