## **A simple pipeline for selection COVI19 amplicons with host specificity and RNA secondary structures prediction**
By Anibal Arce and Ignacio Ibarra

### Motivation
---------------------------------------------------------
Design and selection of primer combinations that generate stable amplicons is desirable for Molecular Biologists. Additionally, filtering amplicons by formation of secondary structures is a step that aids in the shortlisting of amplicons for downstreams detection tools.

### Solution
----------------
This script (i) selects desirable primers for downstream selection, and (ii) uses additional scores to filter them by (a) unwanted specificity with other viral genomes, and (b) selection of

```
python generate_primer_pairs.py --help
python generate_primer_pairs.py --test 100
python generate_primer_pairs.py --others
```

### Workflow steps
1. Primers are scanned against COVID19's GenBank entry (`--fastaid` to modify with custom fasta), 
and filtered by desirable %GC content rules, Tm and low duplex.
2. Primers are scanned against 7 viral genomes, to reduce amplification of untargeted viral genomes.
3. Primers are grouped into pairs and filtered by amplicon lengths.
4. Predicted RNA amplicons are assessed for RNA secondary structures using LinearFold [Huang et al. 2019](https://academic.oup.com/bioinformatics/article/35/14/i295/5529205)). Z-scores are calculated uses a regression model that corrects for the amplicon length, allowing interpretability of values.

### Requirements
- Programming knowledge in Python and Data Science skills are required to run these examples.
- Install [LinearFold](https://github.com/LinearFold/LinearFold). Include to `$PATH` or define path with `--linearfold`.

### Dependencies
- `Python 3` https://www.python.org/
- `pandas numpy seaborn scikit-learn`

#### Installation (typical time: less than 5 minutes)
```git clone primer_design_linearfold_covid19.git
cd primer_design_linearfold_covid19`
python generate_primer_pairs.py --ntest
```

Feedback, errors, or further questions

Please report in Issues, to [Anibal Arce](aaarce@uc.cl) or to [Ignacio Ibarra](ignacio.ibarra@helmholtz-muenchen.de).

### Funding
OpenCovid19 Initiative - Just One Giant Lab.
More information at https://app.jogl.io/project/188