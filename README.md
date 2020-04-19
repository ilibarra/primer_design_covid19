## **Screening of primers for COVID-19 with genome specificity and amplicons with stable RNA secondary structure**
By Anibal Arce and Ignacio Ibarra

### Motivation
---------------------------------------------------------
- Design and selection of primer combinations that specifically recognize target viruses is desirable for Molecular Biologists.
- Filtering amplicons by formation of secondary structures aids in the shortlisting of candidates for downstreams detection tools.

### Solution
----------------

This Python workflow:
1. Selects desirable primers for downstream selection in the COVID-19 genome. Unwanted specificity with other viral genomes is checked. 
2. Amplicon products obtained from primer pairs are checked for RNA secondary structure
stability (reported as free energy and Z-scores).

### Workflow steps
1. Primers are scanned against COVID-19's GenBank entry (`--fastaid` to modify with custom fasta), 
and filtered by desirable %GC content rules, Tm and low duplex.
2. Primers are scanned against 7 viral genomes selected to reduce amplification of untargeted viral genomes:
    - RSV, Ketapneumovirus, Parainfluenza 4a, Adenovirus E, Influenza B, Influenza A.
3. Primers are grouped into pairs and filtered by amplicon lengths.
4. Predicted RNA amplicons are assessed for RNA secondary structures 
using LinearFold ([Huang et al. 2019](https://academic.oup.com/bioinformatics/article/35/14/i295/5529205)). 
Z-scores are calculated using a regression model that models the free energy versus amplicon length. This allows comparing free energies 
of long and short amplicons.

#### Installation and running (typical time: less than 5 minutes)
```
git clone primer_design_linearfold_covid19.git
cd primer_design_linearfold_covid19`
```

### Requirements
- `Python 3` https://www.python.org/
- `pandas numpy seaborn scikit-learn`
- We use [LinearFold](https://github.com/LinearFold/LinearFold) to assess free energy of RNA secondary structures.
- Please include to `$PATH` or define path with `--linearfold`.

# Running examples
```
python generate_primer_pairs.py --help # see help (fasta ID, gc content filter, min-max primer/amplicon lengths, etc.)
python generate_primer_pairs.py --ntest 100 # 
python generate_primer_pairs.py --checkothers # check for other viral genomes
```

# Output
- An Excel file called `output/GCF_009858895.2_CDS_pairs.xlsx` contain details for all shortlisted pairs selected by input criteris.
- Check Z-scores for RNA-stability of primers. More negative values are correlated with stable 2D-RNA amplicons.
    - To visualize, you can submit amplicon sequences in [NUPACK website](http://www.nupack.org/partition/new). 

# Running time
- Around 10-20 minutes for full execution (one CPU, default parameters, verification against other genomes and
RNA secondary structure assessment).
- Adding more background genomes increases running time linearly.
- Reducing the primer lengths increases exponentially the running time.

### Feedback, errors, or further questions
Please report in Issues, to [Anibal Arce](aaarce@uc.cl) or to [Ignacio Ibarra](ignacio.ibarra@helmholtz-muenchen.de).

### Misc
- More background viral genomes can be added manually in `input/other_viruses` as fasta files. Please add also an entry to the file `names.tsv` with a keyword for that name.
- You can replace the FASTA_ID in input to genomes of your interest (not only COVID-19).

### Funding
OpenCovid19 Initiative - Just One Giant Lab.
More information at https://app.jogl.io/project/188