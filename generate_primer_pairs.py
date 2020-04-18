'''
Created on 3/26/2020

Fast analyzer for designing primer pairs with 2D-stable structure

@author: ignacio
'''
from difflib import SequenceMatcher
from itertools import combinations
import tempfile
from path_functions import *
from FastaAnalyzer import FastaAnalyzer
from SequenceMethods import SequenceMethods
from DataFrameAnalyzer import DataFrameAnalyzer
import pandas as pd
import numpy as np

def run(pmin, pmax, amplicon_min, amplicon_max, fasta_id, input_dir, output_dir,
        linearfold_bin, check_others, **kwargs):
    sequences = FastaAnalyzer.get_fastas(join(input_dir, "%s.fasta" % fasta_id), as_dict=True)

    T7 = 'AATTCTAATACGACTCACTATAGGGAGAAGG'
    overwrite = True

    # RULES 1-8 + VIRUSES COMPETITION
    bkp_path_df = join(input_dir, "%s.tsv.gz" % fasta_id)
    if not exists(bkp_path_df) or overwrite:
        table = []
        for h in sequences:
            s = sequences[h]
            for primer_len in range(pmin, pmax):
                print('scanning %s for primers primer of len %i' % (fasta_id, primer_len))
                for si in range(len(s) - primer_len + 1):
                    primer = s[si: si + primer_len]

                    for seq, direction in zip([primer, SequenceMethods.get_complementary_seq(primer)],
                                               ['+', '-']):
                        # GC content 40-60
                        gc = SequenceMethods.get_gc_content(seq)
                        rule1 = gc >= 40 and gc <= 60

                        # GC content first 6 nt
                        gc5prime = SequenceMethods.get_gc_content(seq[:6])
                        gc3prime = SequenceMethods.get_gc_content(seq[-6:])


                        # rule2
                        # Tm = 81.5 + 16.6(log10([Na+])) + .41*(%GC) - 600/length,
                        na_concentration = 0.1
                        Tm = 81.5 + 16.6*(np.log10(na_concentration)) + .41*(gc) - 600/len(seq)
                        rule2 = Tm > 41

                        # homopolymers of length 4
                        rule3 = SequenceMethods.has_homopolymer(seq, 4)

                        rule4 = seq[-1] == 'A'

                        if not rule1 or not rule2 or not rule3 or not rule4:
                            continue
                        table.append([fasta_id, h, si, primer_len, seq, direction, gc, Tm, rule1, rule2, rule3, rule4,
                                      gc5prime, gc3prime])


        df = pd.DataFrame(table, columns=['fasta.id', 'fa.name', 'fasta.position', 'primer.len', 'seq', 'direction',
                                          'GC', 'Tm', 'rule.1', 'rule.2', 'rule.3', 'rule.4', 'gc.5p6nt', 'gc.3p6n'])
        df['k'] = df['fasta.id'] + "_" + df['fasta.position'].astype(str) + "_" + df['primer.len'].astype(str)


        # scan whether primers intersect with other viruses
        other_viruses_dir = join(input_dir, "other_viruses")
        column_names_by_f = DataFrameAnalyzer.get_dict(DataFrameAnalyzer.read_tsv(join(input_dir, 'other_viruses',
                                                                                       'names.tsv')),
                                                       'FILENAME', 'VIRUS')

        if check_others:
            for f in listdir(other_viruses_dir):
                print('filtering negative genome sequences', f)
                best_hits = []
                fa = FastaAnalyzer.get_fastas(join(other_viruses_dir, f))
                for ri, r in df.iterrows():
                    a = r['seq']
                    if ri % 100 == 0:
                        print(ri, f, a)
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

        DataFrameAnalyzer.to_tsv(df, join(output_dir, "%s.tsv.gz" % fasta_id))
        df.to_excel(join(output_dir, "%s.xlsx" % fasta_id))
    df = DataFrameAnalyzer.read_tsv_gz(bkp_path_df)
    best_hit_by_k = DataFrameAnalyzer.get_dict(df, 'k', 'best.hit.others')


    # Analyze group primers by pairs and filter ones that are not good amplicon length min-max
    # and Run LinearFold to get score estimates
    bkp_path_df2 = join(input_dir, "%s_pairs.tsv.gz" % fasta_id)
    if not exists(bkp_path_df2) or overwrite:
        iloc_by_idx = {idx: df.iloc[idx] for idx in df.index}
        faname_by_idx = {idx: df.iloc[idx]['fa.name'] for idx in df.index}
        direction_by_idx = {idx: df.iloc[idx]['direction'] for idx in df.index}

        print('Generating pairs (within ORFs). This is the slowest step (please wait 2-3 min)')
        accepted_pairs = []
        for orf, grp in df.groupby('fa.name'):
            ntest = kwargs.get('ntest')
            print(orf, grp.shape[0])

            accepted_pairs += [[a, b] for a, b in combinations(grp.head(ntest if ntest is not None else grp.shape[0]).index, 2) if
                              (abs(iloc_by_idx[b]['fasta.position'] - iloc_by_idx[a]['fasta.position'] + 1) >= amplicon_min) and
                              (abs(iloc_by_idx[b]['fasta.position'] - iloc_by_idx[a]['fasta.position'] + 1) <= amplicon_max) and
                               (direction_by_idx[a] == '+') and (direction_by_idx[b] == '-')]


        df2 = pd.DataFrame([[a, iloc_by_idx[a]['k'], b, iloc_by_idx[b]['k']] for a, b in accepted_pairs], columns=['i', 'vi', 'j', 'vj'])

        for symbol in ['i', 'j']:
            df2['seq.%s' % symbol] = [iloc_by_idx[r[symbol]]['seq'] for ri, r in df2.iterrows()]
            df2['gc.%s' % symbol] =  df2['seq.%s' % symbol].apply(SequenceMethods.get_gc_content)
            df2['direction.%s' % symbol] = [direction_by_idx[r[symbol]] for ri, r in df2.iterrows()]
            df2['T7.seq.%s' % symbol] = T7 + df2['seq.%s' % symbol]


        assert sum(df2['i'].map(faname_by_idx) != df2['j'].map(faname_by_idx)) == 0

        df2['cds'] = df2['i'].map(faname_by_idx)
        df2['amplicon.len.idx'] = df2['vj'].str.split("_").str[-2].astype(int) - \
                                  df2['vi'].str.split("_").str[-2].astype(int) + 1
        df2 = df2[df2['amplicon.len.idx'] > 0]

        # find the strongest local match between primer pairs
        longest_local_match = []
        for ri, r in df2.iterrows():
            if ri % 100 == 0:
                print("# RUNNING FINDING PRIMERS HITS", ri, 'out of', df2.shape[0])
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
        df2['T7.amplicon.fwd'] = T7 + df2['amplicon.fwd']
        df2['T7.amplicon.rev'] = T7 + df2['amplicon.rev']

        df2['amplicon.len.str'] = df2['vj'].str.split("_").str[-2].astype(int) -\
                                  df2['vi'].str.split("_").str[-2].astype(int) + \
                                  df2['vi'].str.split("_").str[-1].astype(int) + 1

        tmppath = tempfile.mkstemp()[1]
        inpath = tmppath

        queries_linearfold = ['T7.seq.i', 'T7.seq.j',
                              'amplicon.fwd', 'amplicon.rev',
                              'T7.amplicon.fwd', 'T7.amplicon.rev']

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
        # do this twice (i) only for T7+primers and (ii) for T7+longer amplicons
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

        DataFrameAnalyzer.to_tsv_gz(df2, join(output_dir, "%s_pairs.tsv.gz" % fasta_id))
        df2.to_excel(join(output_dir, "%s_pairs.xlsx" % fasta_id), index=None)

if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument("--pmin", type=int, default=20)
    parser.add_argument("--pmax", type=int, default=24)
    parser.add_argument("--ampliconmin", type=int, default=120)
    parser.add_argument("--ampliconmax", type=int, default=240)

    parser.add_argument("--ntest", type=int, default=None)

    parser.add_argument("-p", "--progressbar", action='store_true', default=False)
    parser.add_argument("--others", action='store_true', default=False)

    parser.add_argument('--fastaid', type=str, default='GCF_009858895.2_CDS')
    parser.add_argument('--linearfold', type=str, default='linearfold')

    parser.add_argument('--inputdir', type=str, default="input")
    parser.add_argument('--outputdir', type=str, default="output")

    opts = parser.parse_args()

    run(opts.pmin, opts.pmax, opts.ampliconmin, opts.ampliconmax,
        opts.fastaid, opts.inputdir, opts.outputdir, opts.linearfold, opts.others, ntest=opts.ntest)