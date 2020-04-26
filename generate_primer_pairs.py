'''
Created on 3/26/2020

Fast analyzer for designing primer pairs with 2D-stable structure

@author: ignacio
'''
import tempfile
from difflib import SequenceMatcher
from itertools import combinations
import numpy as np
import pandas as pd
import utilities
from utilities import *
from collections import Counter

# Main script function
def run(pmin, pmax, gcmin, gcmax, tm, amplicon_min, amplicon_max, fasta_id, input_dir, output_dir,
        linearfold_bin, check_others, **kwargs):


    print('Reading MSA data (for mismatches mapping)')
    # calculate iteratively an MSA by calling muscle
    msa_path = join(input_dir, 'lcl_mod_emb-LR757997.1 and 105 other sequences.aln')
    msa_dir = join(input_dir, 'msa')
    if not exists(msa_dir):
        mkdir(msa_dir)
    msa_bkp = join(msa_dir, 'lcl_mod_emb-LR757997.1 and 105 other sequences_n_GCF_009858895.2_CDS.fasta')
    if not exists(msa_bkp):
        sequences = FastaAnalyzer.get_fastas(join(input_dir, "%s.fasta" % fasta_id))
        for h, s in sequences:
            next_fa = join(msa_dir, h.split(" ")[1][1:-1] + ".afa")
            print(exists(next_fa), next_fa)
            if not exists(next_fa):
                fa_in = join(msa_dir, h.split(" ")[1][1:-1] + ".fa")
                FastaAnalyzer.write_fasta_from_sequences([[h, s]], fa_in)
                muscle_cmd = 'muscle -profile -in1 %s -in2 \"%s\" -out %s' % (fa_in, msa_path, next_fa)
                print(muscle_cmd)
                system(muscle_cmd)

    # read msa info: This dictionary will be referred to during the last step
    variants_by_h = {}
    for f in listdir(msa_dir):
        if not f.endswith('.afa'):
            continue
        code = f.replace(".afa", '')
        curr_position = 0
        msa = FastaAnalyzer.get_fastas(join(msa_dir, f))
        for si, nt in enumerate(msa[0][1]):
            if nt != '-':
                if not code in variants_by_h:
                    variants_by_h[code] = {}
                variants_by_h[code][curr_position] = Counter([s2[si] for h2, s2 in msa])
                curr_position += 1

    hits_by_h2 = {}
    for h, s in sequences:
        print(h, len(s))
        for h2, s2 in msa:
            if s in s2.replace("-", ''):
                print(h, h2, 'found')
                hits_by_h2[h2] = 1 if not h2 in hits_by_h2 else hits_by_h2[h2] + 1



    print('Primer pairs generator + with background genome check (step 1) and RNA secondary structure check (step 2)')
    sequences = FastaAnalyzer.get_fastas(join(input_dir, "%s.fasta" % fasta_id), as_dict=True)
    tagseq = kwargs.get('tagseq', '')
    overwrite1 = kwargs.get('overwrite1', 0)
    overwrite2 = kwargs.get('overwrite2', 0)



    # RULES 1-8 + VIRUSES COMPETITION
    print('STEP 1')
    bkp_path_df = join(output_dir, "%s.tsv.gz" % fasta_id)

    if not exists(bkp_path_df) or overwrite1:
        table = []
        for h in sequences:
            s = sequences[h]
            for primer_len in range(pmin, pmax):
                print('scanning %s for primers primer of len %i' % (fasta_id, primer_len))
                for si in range(len(s) - primer_len + 1):
                    primer = s[si: si + primer_len]

                    for seq, direction in zip([primer, SequenceMethods.get_complementary_seq(primer)],
                                               ['+', '-']):
                        # GC content between gcmin and gcmax
                        gc = SequenceMethods.get_gc_content(seq)
                        rule1 = gc >= gcmin and gc <= gcmax

                        # GC content of first 6 nt
                        gc5prime = SequenceMethods.get_gc_content(seq[:6])
                        gc3prime = SequenceMethods.get_gc_content(seq[-6:])


                        # rule2
                        # Tm = 81.5 + 16.6(log10([Na+])) + .41*(%GC) - 600/length,
                        na_concentration = 0.1
                        primer_tm = 81.5 + 16.6*(np.log10(na_concentration)) + .41*(gc) - 600/len(seq)
                        rule2 = primer_tm > tm

                        # homopolymers of length 4
                        rule3 = SequenceMethods.has_homopolymer(seq, 4)

                        rule4 = seq[-1] == 'A'

                        if not rule1 or not rule2 or not rule3 or not rule4:
                            continue
                        table.append([fasta_id, h, si, primer_len, seq, direction, gc, primer_tm,
                                      rule1, rule2, rule3, rule4, gc5prime, gc3prime])


        df = pd.DataFrame(table, columns=['fasta.id', 'fa.name', 'fasta.position', 'primer.len', 'seq', 'direction',
                                          'GC', 'Tm', 'rule.1', 'rule.2', 'rule.3', 'rule.4', 'gc.5p6nt', 'gc.3p6n'])
        df['k'] = df['fasta.id'] + "_" + df['fasta.position'].astype(str) + "_" + df['primer.len'].astype(str)


        # check for MSA alternate variants within the primers
        msa_flagged_primers = []
        flag_details = []
        for ri, r in df.iterrows():
            flag = False
            flag_details.append('')
            k = r['fa.name'].split(" ")[1][1:-1]
            for pi in range(r['fasta.position'], r['fasta.position'] + r['primer.len']):
                n_variants = len({nt for nt in variants_by_h[k][pi].keys() if nt not in {'N'}})
                if n_variants >= 2:
                    flag = True
                    flag_details[-1] += str(pi) + ":" + str(dict(variants_by_h[k][pi])) + ";"
            msa_flagged_primers.append(flag)
        df['msa.flagged.primers'] = msa_flagged_primers
        df['msa.flagged.primers.desc'] = flag_details
        df['n.msa.flagged.primers.desc'] = [len(x.split(";")) - 1 for x in df['msa.flagged.primers.desc']]


        # scan whether primers intersect with other viruses
        other_viruses_dir = join(input_dir, "other_viruses")
        column_names_by_f = DataFrameAnalyzer.get_dict(DataFrameAnalyzer.read_tsv(join(input_dir, 'other_viruses',
                                                                                       'names.tsv')),
                                                       'FILENAME', 'VIRUS')

        # check against all background viral genomes
        if check_others:
            for f in listdir(other_viruses_dir):
                if f.endswith('.tsv'):
                    continue
                print('Scanning against background viral genomes using file ...', f, column_names_by_f[f])
                best_hits = []
                fa = FastaAnalyzer.get_fastas(join(other_viruses_dir, f))
                for ri, r in df.iterrows():
                    a = r['seq']
                    if ri % 200 == 0:
                        print(ri, 'primers out of', df.shape[0], 'matched against', column_names_by_f[f], a)
                    cmp_a = SequenceMethods.get_complementary_seq(a)
                    n_best_match = 0
                    for h, b in fa:
                        match_ab = SequenceMatcher(None, a, b,
                                                   autojunk=False).find_longest_match(0, len(a), 0, len(b))
                        match_cmpa_b = SequenceMatcher(None, cmp_a, b,
                                                       autojunk=False).find_longest_match(0, len(cmp_a), 0, len(b))
                        if match_ab.size == 0 and match_cmpa_b.size == 0:
                            print('problem with sequence matches. Please check')
                            assert 1 > 2
                        n_best_match = max(match_ab.size, match_cmpa_b.size, n_best_match)

                    best_hits.append(n_best_match)
                df[f] = best_hits
            df = df.rename(columns=column_names_by_f)
            df['best.hit.others'] = df[[c for c in df if c in column_names_by_f.values()]].max(axis=1)
        else:
            df['best.hit.others'] = -1


        print('Saving selected primers:')
        DataFrameAnalyzer.to_tsv(df, join(output_dir, "%s.tsv.gz" % fasta_id))
        df.to_excel(join(output_dir, "%s.xlsx" % fasta_id))
    else:
        print('skip STEP 1 (file exists and overwrite1=False')
    df = DataFrameAnalyzer.read_tsv_gz(bkp_path_df)
    best_hit_by_k = DataFrameAnalyzer.get_dict(df, 'k', 'best.hit.others')


    # Analyze group primers by pairs and filter ones that are not good amplicon length min-max
    # and Run LinearFold to get score estimates
    print('STEP 2')
    print('\nSelection of 1-2 primer pairs (amplicon length + RNA secondary structure)...')
    bkp_path_df2 = join(output_dir, "%s_pairs.tsv.gz" % fasta_id)
    if not exists(bkp_path_df2) or overwrite2:
        iloc_by_idx = {idx: df.iloc[idx] for idx in df.index}
        faname_by_idx = {idx: df.iloc[idx]['fa.name'] for idx in df.index}
        direction_by_idx = {idx: df.iloc[idx]['direction'] for idx in df.index}

        print('Generating primer pairs (within ORFs). This is the slowest step (please wait 3-5 in all ORFs...)')
        accepted_pairs = []
        for orf, grp in df.groupby('fa.name'):
            ntest = kwargs.get('ntest')
            print('Generating primer pairs for ORF:', orf[:100], "... # primers=%i" % grp.shape[0])
            accepted_pairs += [[a, b] for a, b in combinations(grp.head(ntest if ntest is not None else grp.shape[0]).index, 2) if
                              (abs(iloc_by_idx[b]['fasta.position'] - iloc_by_idx[a]['fasta.position'] + 1) >= amplicon_min) and
                              (abs(iloc_by_idx[b]['fasta.position'] - iloc_by_idx[a]['fasta.position'] + 1) <= amplicon_max) and
                               (direction_by_idx[a] == '+') and (direction_by_idx[b] == '-')]


        df2 = pd.DataFrame([[a, iloc_by_idx[a]['k'], b, iloc_by_idx[b]['k']] for a, b in accepted_pairs], columns=['i', 'vi', 'j', 'vj'])

        for symbol in ['i', 'j']:
            df2['seq.%s' % symbol] = [iloc_by_idx[r[symbol]]['seq'] for ri, r in df2.iterrows()]
            df2['gc.%s' % symbol] =  df2['seq.%s' % symbol].apply(SequenceMethods.get_gc_content)
            df2['direction.%s' % symbol] = [direction_by_idx[r[symbol]] for ri, r in df2.iterrows()]
            df2['tag.seq.%s' % symbol] = tagseq + df2['seq.%s' % symbol]


        assert sum(df2['i'].map(faname_by_idx) != df2['j'].map(faname_by_idx)) == 0

        df2['cds'] = df2['i'].map(faname_by_idx)
        df2['amplicon.len.idx'] = df2['vj'].str.split("_").str[-2].astype(int) - \
                                  df2['vi'].str.split("_").str[-2].astype(int) + 1
        df2 = df2[df2['amplicon.len.idx'] > 0]

        # find the strongest local match between primer pairs
        longest_local_match = []
        for ri, r in df2.iterrows():
            if ri % 100 == 0:
                print("# Scanning for local primer pair hits", ri, 'out of', df2.shape[0])
            a, b = r['seq.i'], r['seq.j']
            cmp_b = SequenceMethods.get_complementary_seq(b)
            match_a_cmpb = [a, cmp_b, SequenceMatcher(None, a, cmp_b, autojunk=False).find_longest_match(0, len(a), 0, len(b))]
            longest_local_match.append(a[match_a_cmpb[-1].a: match_a_cmpb[-1].a + match_a_cmpb[-1].size] + "/" + \
                                       SequenceMethods.get_complementary_seq(cmp_b[match_a_cmpb[-1].b:  match_a_cmpb[-1].b +
                                                                                                        match_a_cmpb[-1].size]))

        df2['longest.local.match'] = longest_local_match
        df2['longest.local.match.len'] = df2['longest.local.match'].str.split("/").str[0].str.len()


        # map the best hit with other viruses
        df2['best.hit.others'] = [max(a, b) for a, b in zip(list(df2['vi'].map(best_hit_by_k)),
                                                            list(df2['vj'].map(best_hit_by_k)))]

        df2 = df2[df2['longest.local.match.len'].abs() <= 5]

        df2['amplicon.fwd'] = [sequences[r['cds']][int(r['vi'].split("_")[-2]):int(r['vj'].split("_")[-2]) +
                                                                               int(r['vj'].split("_")[-1])]
                               for ri, r in df2.iterrows()]
        df2['amplicon.rev'] = df2['amplicon.fwd'].apply(SequenceMethods.get_complementary_seq)
        df2['tag.amplicon.fwd'] = tagseq + df2['amplicon.fwd']
        df2['tag.amplicon.rev'] = tagseq + df2['amplicon.rev']

        df2['amplicon.len.str'] = df2['vj'].str.split("_").str[-2].astype(int) -\
                                  df2['vi'].str.split("_").str[-2].astype(int) + \
                                  df2['vi'].str.split("_").str[-1].astype(int) + 1

        tmppath = tempfile.mkstemp()[1]
        inpath = tmppath

        queries_linearfold = ['tag.seq.i', 'tag.seq.j',
                              'amplicon.fwd', 'amplicon.rev',
                              'tag.amplicon.fwd', 'tag.amplicon.rev']

        linearfold = '%s -V' % linearfold_bin
        inpath = tmppath + ".in"

        # SLOWEST STEP. RUN ONCE PER SEQUENCE SET
        for qi, label in enumerate(queries_linearfold):
            scores_col = 'dG.LFold.%s' % label
            if scores_col in df2:
                continue
            print(qi, label)
            DataFrameAnalyzer.write_list(df2[label].replace("T", "U"), inpath)
            print('running LinearFold with columns %s (# queries=%i)' % (label, df2.shape[0]))
            out_linearfold = os.popen('cat %s | %s' % (inpath, linearfold)).read()
            scores = [float(r.split(" ")[1].replace("(", '').replace(")", ''))
                      for r in out_linearfold.split("\n") if '.' in r]
            df2['dG.LFold.%s' % label] = scores

        # calculate Z-scores based on mean by length
        # do this twice (i) only for tag+primers and (ii) for tag+primer+amplicon
        # mean is estimated by linear model
        # more negative Z-scores = less reliable primers/amplicons due to unexpected stability
        for sub_queries in queries_linearfold[:2], queries_linearfold[2:]:
            z_df = []
            for qi, q in enumerate(sub_queries):
                print(q, q in df2)
                if not "dG.LFold.%s" % q in df2:
                    continue
                sel = df2[[q, "dG.LFold.%s" % q]]
                sel.columns = ['len', 'dG']
                sel['len'] = sel['len'].str.len()
                print(sel.head())
                z_df.append(sel)
            z_df = pd.concat(z_df).reset_index(drop=True)
            from sklearn.linear_model import LinearRegression
            # calculate a coefficient to estimate expected dG by len

            print('generate Z-score using linear model')
            lm = LinearRegression().fit(np.array(z_df[['len']]), np.array(z_df['dG']))
            sigma = float(z_df.groupby('len').std().mean())

            for qi, q in enumerate(sub_queries):
                df2['z.score.dG.%s' % q] = (df2['dG.LFold.%s' % q] - lm.predict(np.array(df2[q].str.len()).reshape(-1, 1))) / sigma



        # check for MSA alternate variants within the primers
        msa_flagged_primers = []
        flag_details = []
        for ri, r in df2.iterrows():
            flag = False
            flag_details.append('')
            k = r['cds'].split(" ")[1][1:-1]
            start, end, primer_len = int(r["vi"].split("_")[-2]), int(r["vj"].split("_")[-2]), int(r["vi"].split("_")[-1])
            for pi in range(start, end + primer_len):
                n_variants = len({nt for nt in variants_by_h[k][pi].keys() if nt not in {'N'}})
                if n_variants >= 2:
                    flag = True
                    flag_details[-1] += str(pi) + ":" + str(dict(variants_by_h[k][pi])) + ";"
            msa_flagged_primers.append(flag)
        df2['msa.flagged.primers'] = msa_flagged_primers
        df2['msa.flagged.primers.desc'] = flag_details
        df2['n.msa.flagged.primers.desc'] = [len(x.split(";")) - 1 for x in df2['msa.flagged.primers.desc']]


        DataFrameAnalyzer.to_tsv_gz(df2, join(output_dir, "%s_pairs.tsv.gz" % fasta_id))
        df2.to_excel(join(output_dir, "%s_pairs.xlsx" % fasta_id), index=None)
    else:
        print('skip STEP 2 (file exists and overwrite2=False')

if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument("--pmin", type=int, default=20, help='minimum primer length (def. 20)')
    parser.add_argument("--pmax", type=int, default=24, help='minimum primer length (def. 24)')

    parser.add_argument("--gcmin", type=float, default=40,
                        help='minimum GC content for primers (def. 40 Celsius).')
    parser.add_argument("--gcmax", type=float, default=60,
                        help='maximum GC content for primers (def. 60 Celsius).')

    parser.add_argument("--tmmin", type=float, default=41,
                        help='minimum melting temperature (def. 41 Celsius)')

    parser.add_argument("--ampliconmin", type=int, default=120, help='minimum amplicon length (def. 120)')
    parser.add_argument("--ampliconmax", type=int, default=240, help='minimum amplicon length (def. 240)')

    # Use the T7 sequence as a default tag
    parser.add_argument('--tagprimer', type=str, help='tag for primers (def. T7 sequence, AATTCTAATACGACTCACTATAGGGAGAAGG)',
                        default='AATTCTAATACGACTCACTATAGGGAGAAGG')

    parser.add_argument("--ntest", type=int, default=None, help='for load tests. Default is None (--ntest 10 = test for 10 primer pairs and finish')
    parser.add_argument("--overwrite1", action='store_true', help='Force repeat single primer generation and background viruses scanning step', default=0)
    parser.add_argument("--overwrite2", action='store_true', help='Force repeat 1-2 primer pairs and secondary structure asssessment', default=0)

    parser.add_argument("-p", "--progressbar", action='store_true', default=False,
                        help='Show progress bar (not implemented in deployed version).')
    parser.add_argument("--checkothers", action='store_true', default=False)

    parser.add_argument('--fastaid', type=str, default='GCF_009858895.2_CDS', help='fastaid to use from input dir (def. GCF_009858895.2_CDS)')
    parser.add_argument('--linearfold', type=str, default='linearfold', help='path to linearfold if not declared in $PATH')

    parser.add_argument('--inputdir', type=str, default="input", help='input directory (def. ./input)')
    parser.add_argument('--outputdir', type=str, default="output", help='output directory (def. ./output)')

    opts = parser.parse_args()

    run(opts.pmin, opts.pmax, opts.gcmin, opts.gcmax, opts.tmmin, opts.ampliconmin, opts.ampliconmax,
        opts.fastaid, opts.inputdir, opts.outputdir, opts.linearfold, opts.checkothers, ntest=opts.ntest, tag=opts.tagprimer,
        overwrite1=opts.overwrite1, overwrite2=opts.overwrite2)