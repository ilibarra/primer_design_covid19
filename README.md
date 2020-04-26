<img src="./about/jogl_logo.png" width="180px" height="250px" align="right">

## **Screening of primers for COVID-19 with genome specificity and amplicons with stable single-stranded RNA structure**
By Anibal Arce and Ignacio Ibarra

### Motivation
---------------------------------------------------------
- Sequence-specific isothermal amplification of targets from viral samples is essential for implementing low-cost diagnostic tools based on nucleic acid detection.
- Design and selection of primer combinations for this propuse can be automated considering the optimal parameters for a particular methodology, such as NASBA isothermal amplification of RNA.
- Analysis of the secondary structure of the amplicon generated, aids in the shortlisting of candidates for downstreams detection tools such as Toeholds RNA Sensors.

### Solution
----------------

This Python workflow:
1. Selects desirable primers for experiments targeting in the COVID-19 genome. Unwanted specificity with other viral genomes is checked and reported. 
2. Amplicon products obtained from primer pairs are checked for RNA secondary structure
stability (reported as minimum free energy and Z-scores).

### Workflow steps
1. Primers are scanned against COVID-19's GenBank entry (`--fastaid` to modify with custom fasta), 
and filtered by desirable %GC content rules, Tm and local alignment between primers (duplex formation).
2. Primers are scanned against 7 viral genomes selected to reduce amplification of untargeted viral genomes:
    - RSV, Ketapneumovirus, Parainfluenza 4a, Adenovirus E, Influenza B, Influenza A.
3. Primers are grouped into pairs and filtered by amplicon lengths.
4. Predicted RNA amplicons are assessed for RNA secondary structures 
using LinearFold ([Huang et al. 2019](https://academic.oup.com/bioinformatics/article/35/14/i295/5529205)). 
Z-scores are calculated using a regression model that models the free energy versus amplicon length. This allows comparing free energies 
of long and short amplicons.

### Installation and running (typical time: less than 5 minutes)
```
git clone primer_design_linearfold_covid19.git
cd primer_design_linearfold_covid19
```

### Environment requirements
- `Python 3` https://www.python.org/
- Data Science packages for Python: `pandas numpy seaborn scikit-learn`
- [LinearFold](https://github.com/LinearFold/LinearFold) must be installed to assess free energy of RNA secondary structures. Please include to `$PATH` or define path with `--linearfold`.
- **New (April/2020)** MUSCLE : Necessary to generate MSA files between reference FASTA and included aln for virus of interest.

### Execution examples
```
python generate_primer_pairs.py --help # see help (fasta ID, gc content filter, min-max primer/amplicon lengths, etc.)
# LOAD TESTS
python generate_primer_pairs.py --ntest 100 --overwrite2 # test only with the first 100 primer pairs (run in 5 minutes).
python generate_primer_pairs.py --ntest 1000 --overwrite2 # test only with the first 1000 pairs.
# FULL RUN
python generate_primer_pairs.py --overwrite2 # full run (checking RNA secondary structures).
python generate_primer_pairs.py --checkothers --overwrite1 --overwrite2 # full run (checking against background viral genomes).
```

### Output
- An Excel file called `output/FASTA_ID_pairs.xlsx` contains details for all shortlisted primer pairs and amplicons selected.
- Check Z-scores for RNA-stability of primers and amplicons. More negative values are correlated with more stable RNA secondary structures.
    - You can visualize 2D-stable amplicon sequences using [NUPACK](http://www.nupack.org/partition/new). 

### Running time
- Around 5-15 minutes for load tests between 100 and 2000 (primer pairs and their respective amplicons).
- ~90-120 minutes for full execution (one CPU, default parameters, verification against other genomes and
RNA secondary structure assessment).
- Adding more background genomes increases running time linearly.
- Reducing the primer lengths increases exponentially the running time.

### Misc
- More background viral genomes can be added manually in `input/other_viruses` as fasta files. Please add also an entry to the file `names.tsv` with a keyword for that name.
- You can replace the FASTA_ID in input to genomes of your interest (not only COVID-19).

### Feedback, errors, or further questions
- Report in Issues.
- You can also Send an e-mail to Anibal Arce (aaarce@uc.cl) or Ignacio Ibarra (ignacio.ibarra@helmholtz-muenchen.de).

### Funding
OpenCovid19 Initiative - Just One Giant Lab.
More information at https://app.jogl.io/project/188
