
from random import shuffle, randint, seed
from itertools import combinations, product
import pandas as pd
from functools import reduce

class SequenceMethods:

    @staticmethod
    def get_complementary_nt(nt):
        if nt == "A":
            return "T"
        if nt == "T":
            return "A"
        if nt == "C":
            return "G"
        if nt == "G":
            return "C"
        if nt == "a":
            return "t"
        if nt == "t":
            return "a"
        if nt == "c":
            return "g"
        if nt == "g":
            return "c"
        if nt == "R":
            return "Y"
        if nt == "Y":
            return "R"
        if nt == "K":
            return "M"
        if nt == "M":
            return "K"
        if nt == "S":
            return "S"
        if nt == "W":
            return "W"
        if nt == "B":
            return "V"
        if nt == "V":
            return "B"
        if nt == "D":
            return "H"
        if nt == "H":
            return "D"
        if nt == "N":
            return "N"
        if nt == "[":
            return ']'
        if nt == "]":
            return '['
        return nt

    @staticmethod
    def get_complementary_seq(s):
        return "".join([SequenceMethods.get_complementary_nt(nt) for nt in s][::-1])

    @staticmethod
    def parse_coordinate2range(df, series=None, names=['chr', 'start', 'end'], sep1=':', sep2='-'):
        # format example:
        # """chr:start-end"""
        if isinstance(df, pd.DataFrame):
            if series is None:
                series = df.index
            df[names[0]] = series.str.split(sep1).str[0]
            if sep1 != sep2:
                df[names[1]] = (series.str.split(sep1).str[1].str.split(sep2).str[0]).astype(int)
                df[names[2]] = (series.str.split(sep1).str[1].str.split(sep2).str[1]).astype(int)
            else:
                df[names[1]] = (series.str.split(sep1).str[1]).astype(int)
                df[names[2]] = (series.str.split(sep1).str[2]).astype(int)
            return df
        elif isinstance(df, str):
            return [df.split(sep1)[0]] + list(map(int, df.split(sep1)[1].split(sep2) if sep1 != sep2 else df.split(sep1)[1:]))


    @staticmethod
    def parse_range2coordinate(df, names=['chr', 'start', 'end'], k='range'):
        if isinstance(df, pd.DataFrame):
            df[k] = df[names[0]].astype(str) + ":" +\
                    df[names[1]].astype(int).astype(str) + "-" + df[names[2]].astype(int).astype(str)
            return df
        elif isinstance(df, list):
            return [t[0] + ":" + str(t[1]) + "-" + str(t[2]) for t in df]

    @staticmethod
    def get_kmers_list(s, k):
        return [s[si: si + k] for si in range(0, len(s) - k + 1)]
    @staticmethod
    def get_iupac_table():
        # all the IUPAC associations
        iupac_table = {"A": {"A"}, "C": {"C"}, "T": {"T"}, "G": {"G"},
                       "R": {"A", 'G'}, "Y": {"C", 'T'},
                       "W": {"A", 'T'}, "S": {"C", 'G'},
                       "K": {"G", 'T'}, "M": {"A", 'C'},
                       "B": {"C", 'G', 'T'},
                       "D": {"A", 'G', 'T'},
                       "H": {"A", 'C', 'T'},
                       "V": {"A", 'C', 'G'},
                       "N": {"C", 'G', 'T', 'A'}}
        return iupac_table

    @staticmethod
    def get_matches(query, target, get_mistmaches_pos=False, max_score=False, ignore_n=True):
        a, b = query, target

        if isinstance(b, tuple):
            b = list(b)
        if not isinstance(b, list) and not isinstance(b, tuple):
            assert len(a) <= len(b)

        if isinstance(b, list) or len(b) > len(a):
            from lib.SELEX.SELEXAnalyzer import SELEXAnalyzer
            selex = SELEXAnalyzer()
            queries = selex.mapped_motifs_in_sequences([b] if not isinstance(b, list) else b, len(a))
            # print 'scoring using kmer...'
            queries['n.matches'] = queries['seq'].apply(SequenceMethods.get_matches, args=[a])
            # print 'done...'
            if max_score:
                queries = queries.sort_values('n.matches', ascending=False)
                queries = queries.drop_duplicates('peak.id')
                queries = queries.sort_values('peak.id', ascending=True)
            return queries

        # all the IUPAC associations
        iupac_table = SequenceMethods.get_iupac_table()
        counter = 0
        i = 0
        mismatches_positions = set()
        for ai, bi in zip(a, b):
            if ignore_n and bi == 'N':
                continue
            match = False
            if ai == bi:
                counter += 1
                match = True
            else:
                if SequenceMethods.is_watson_crick_nt(ai) and SequenceMethods.is_watson_crick_nt(bi):
                    match = False
                elif SequenceMethods.is_watson_crick_nt(ai) and not SequenceMethods.is_watson_crick_nt(bi):
                    if ai in iupac_table[bi]:
                        counter += 1
                        match = True
                elif not SequenceMethods.is_watson_crick_nt(ai) and SequenceMethods.is_watson_crick_nt(bi):
                    if bi in iupac_table[ai]:
                        counter += 1
                        match = True
            i += 1
            if not match:
                mismatches_positions.add(i)
        if get_mistmaches_pos:
            return counter, mismatches_positions
        else:
            return counter

    @staticmethod
    def get_js_divergence(input_msa_fasta):
        import subprocess as sp
        import tempfile
        out_path = tempfile.mkstemp()[1]
        script_path = 'score_conservation.py'
        args = " ".join([script_path, input_msa_fasta, '>', out_path])
        print(args)
        system(args)
        cons= DataFrameAnalyzer.read_tsv(out_path, skiprows=1)
        cons['score.masked'] = np.where(cons['score'] == -1000, np.nan, cons['score'])
        remove(out_path)
        return cons

    @staticmethod
    def levenshtein_iterative(s, t, costs=(1, 1, 1)):
        """
            iterative_levenshtein(s, t) -> ldist
            ldist is the Levenshtein distance between the strings
            s and t.
            For all i and j, dist[i,j] will contain the Levenshtein
            distance between the first i characters of s and the
            first j characters of t

            costs: a tuple or a list with three integers (d, i, s)
                   where d defines the costs for a deletion
                         i defines the costs for an insertion and
                         s defines the costs for a substitution
        """
        rows = len(s) + 1
        cols = len(t) + 1
        deletes, inserts, substitutes = costs

        dist = [[0 for x in range(cols)] for x in range(rows)]
        # source prefixes can be transformed into empty strings
        # by deletions:
        for row in range(1, rows):
            dist[row][0] = row * deletes
        # target prefixes can be created from an empty source string
        # by inserting the characters
        for col in range(1, cols):
            dist[0][col] = col * inserts

        for col in range(1, cols):
            for row in range(1, rows):
                if s[row - 1] == t[col - 1]:
                    cost = 0
                else:
                    cost = substitutes
                dist[row][col] = min(dist[row - 1][col] + deletes,
                                     dist[row][col - 1] + inserts,
                                     dist[row - 1][col - 1] + cost)  # substitution
        # for r in range(rows):
        #     print(dist[r])

        return dist[row][col]

    @staticmethod
    def get_gc_content(seq):
        return sum([nt == 'G' or nt == 'C' for nt in seq]) / float(len(seq)) * 100

    @staticmethod
    def has_homopolymer(s, k):
        for nt in 'ACGT':
            next = nt * k
            if next in s:
                return True
        return False

    @staticmethod
    def get_sequence_combinations(length, dict_nt="ACGT"):
        for output in product(dict_nt, repeat=length):
            yield(''.join(output))

    @staticmethod
    def get_best_alignment(query, target, max_shift=None):
        k = 0
        table = []

        # left shift of the target sequence
        # print 'left'
        for i in range(len(query)):
            if max_shift is not None and i > max_shift:
                break
            a = [qi for qi in query]
            b = ['' for bi in range(i)] + [ti for ti in target]
            a, b = list(zip(*[[ai, bi] for ai, bi in zip(a, b) if ai != '' and bi != '']))
            a = "".join(a)
            b = "".join(b)
            t = [-i, a, b, SequenceMethods.get_matches(a, b)]
            table.append(t)
        # right shift of the target sequence
        # print 'right'
        for i in range(len(query) + len(target) - 1):
            if i == 0:
                continue
            a = ['' for xi in range(i)] + [qi for qi in query]
            b = [ti for ti in target] + ['' for xi in range(i)]
            ali = [[ai, bi] for ai, bi in zip(a, b) if ai != '' and bi != '']
            a = "".join([ai[0] for ai in ali])
            b = "".join([bi[1] for bi in ali])

            if len(a) == 0 or len(b) == 0:
                continue
            t = [i, a, b, SequenceMethods.get_matches(a, b)]
            table.append(t)

        res = pd.DataFrame(table, columns=['shift', 'query', 'target', 'matches'])
        return res

    @staticmethod
    def get_best_alignment_seq_vs_pwm(query, pwm, max_shift=None):
        k = 0
        table = []

        # left shift of the target sequence
        # print 'left'
        # print pwm
        for i in range(len(query)):
            if max_shift is not None and i > max_shift:
                break
            a = [qi for qi in query]
            b = [-1 for bi in range(i)] + [pwm[bi] for bi in range(len(pwm.columns))]
            a, b = list(zip(*[[ai, bi] for ai, bi in zip(a, b) if ai != '' and isinstance(bi, pd.Series)]))
            a = "".join(a)
            t = [-i, a, 'pwm', sum([bi[ai] for ai, bi in zip(a, b)])]
            table.append(t)

        # right shift of the target sequence
        # print 'right'
        for i in range(len(query) + len(pwm.columns) - 1):
            if i == 0:
                continue
            if max_shift is not None and i > max_shift:
                break
            a = ['' for xi in range(i)] + [qi for qi in query]
            b = [pwm[bi] for bi in range(len(pwm.columns))] + [-1 for xi in range(i)]
            ali = [[ai, bi] for ai, bi in zip(a, b) if ai != '' and isinstance(bi, pd.Series)]
            a = "".join([ai[0] for ai in ali])
            if len(a) == 0 or len(b) == 0:
                continue
            t = [i, a, 'pwm', sum([bi[ai] for ai, bi in zip(a, b)])]
            table.append(t)

        res = pd.DataFrame(table, columns=['shift', 'query', 'target', 'matches'])
        return res

    @staticmethod
    def get_best_alignment_ppm_vs_ppm(query, target, max_shift=None, get_max=False):
        k = 0
        table = []

        # left shift of the target sequence
        # print 'left'
        # print pwm
        for i in range(len(query.columns)):
            if max_shift is not None and abs(i) > max_shift:
                break
            # print i
            a = [query[ai] for ai in range(len(query.columns))]
            b = [-1 for bi in range(i)] + [target[bi] for bi in range(len(target.columns))]
            a, b = list(zip(*[[ai, bi] for ai, bi in zip(a, b) if isinstance(ai, pd.Series) and isinstance(bi, pd.Series)]))
            # distance = sum([sum(((ai - bi) ** 2)) for ai, bi in zip(a, b)])

            # this is wrong. It has to be MULTIPLIED by position
            distance = reduce(lambda x, y: x * y, [(ai * bi).sum() for ai, bi in zip(a, b)])
            t = [-i, 'ppm1', 'ppm2', distance / len(a)]
            table.append(t)
            # print table[-1]

        # right shift of the target sequence
        # print 'right'
        for i in range(len(query.columns) + len(target.columns) - 1):
            if i == 0:
                continue
            if max_shift is not None and abs(i) > max_shift:
                break
            # print i
            a = [-1 for xi in range(i)] + [query[ai] for ai in range(len(query.columns))]
            b = [target[bi] for bi in range(len(target.columns))] + [-1 for xi in range(i)]
            # print list(a)
            # print list(b)
            if len([ai for ai, bi in zip(a, b) if isinstance(ai, pd.Series) and isinstance(bi, pd.Series)]) == 0:
                continue

            a, b = list(zip(*[[ai, bi] for ai, bi in zip(a, b) if isinstance(ai, pd.Series) and isinstance(bi, pd.Series)]))
            if len(a) == 0 or len(b) == 0:
                continue
            distance = reduce(lambda x, y: x * y, [(ai * bi).sum() for ai, bi in zip(a, b)])
            t = [i, 'ppm1', 'ppm2', distance / len(a)]
            table.append(t)
            # print table[-1]

        if get_max:
            best = max([t[-1] for t in table])
            table = [t for t in table if t[-1] == best]
        res = pd.DataFrame(table, columns=['shift', 'query', 'target', 'distance'])
        return res



    @staticmethod
    def is_watson_crick_nt(nt):
        return nt in {'A', 'C', 'G', 'T'}

    @staticmethod
    def randomize_sequence(sequence):
        nucleotides = [nt for nt in sequence]
        shuffle(nucleotides)
        return "".join(nucleotides)

    @staticmethod
    def get_one_hot_encoding(seq):
        return [int(ci == k) for ci in seq for k in ['A', 'C', 'G', 'T']]


    @staticmethod
    def get_one_hot_encoding_dinucl(seq):
        dinucleotides = [a + b for a in ['A', 'C', 'G', 'T'] for b in ['A', 'C', 'G', 'T']]
        return [int(seq[i: i + 2] == k) for i in range(0, len(seq) - 1) for k in dinucleotides]

    @staticmethod
    def get_one_hot_encoding_trinucl(seq):
        nucl = ['A', 'C', 'G', 'T']
        trinucleotides = [a + b + c for a in nucl for b in nucl for c in nucl]
        return [int(seq[i: i + 3] == k) for i in range(0, len(seq) - 2) for k in trinucleotides]

    @staticmethod
    def get_len_regex(seq):
        seq = seq.replace("[", "").replace("]", '')
        while ',' in seq:
            pos = seq.find(',')
            seq = seq[:pos - 1] + "N" + seq[pos + 2:]
        return len(seq)

    @staticmethod
    def get_dna_code():
        return 'ACGT'

    @staticmethod
    def get_protein_code():
        return 'ACDEFGHIKLMNPQRSTVWY'

    @staticmethod
    def get_random_sequence(length, code=None):
        if code is None:
            code = SequenceMethods.get_dna_code()

        return "".join([code[randint(0, len(code) - 1)] for i in range(length)])

    @staticmethod
    def get_mutated_sequence(s, n=1):
        for ni in range(n):
            nt = "ATGC"[randint(0, 3)]
            pos = randint(0, len(s) - len(nt))
            s = s[:pos] + nt + s[pos + len(nt):]
        return s

    @staticmethod
    def get_mutated_sequence_prot(s, n=1, protected_aa=None):
        code = SequenceMethods.get_protein_code()

        old_s = s
        new_s = None
        for ni in range(n):
            changed = False
            new_s = s
            while old_s == new_s:
                aa = code[randint(0, len(code) - 1)]
                pos = randint(0, len(s) - len(s))
                aa_pos = old_s[pos]
                if protected_aa is not None and aa_pos in protected_aa:
                    continue
                new_s = old_s[:pos] + aa + old_s[pos + len(aa):]
            old_s = new_s
        return new_s

    @staticmethod
    def get_mutated_sequence_prob(s, mutation_prob):
        n = len(s) * mutation_prob
        return SequenceMethods.get_mutated_sequence(s, int(n))

    @staticmethod
    def mutate_at_positions(s, positions):
        new_s = ''
        for ci, c in enumerate(s):
            options = {nt for nt in "ATGC" if nt != c}
            new_s += list(options)[randint(0, 2)] if ci in positions else c
        return new_s

    @staticmethod
    def count_mismatches(s, other):
        return sum([ai != bi for ai, bi in zip(s, other)])

    @staticmethod
    def is_purine(nt):
        return nt == 'A' or nt == 'G'

    @staticmethod
    def is_pyrimidine(nt):
        return nt == 'T' or nt == 'C'

    @staticmethod
    def is_transition(a, b):
        assert a != b
        return (SequenceMethods.is_purine(a) and SequenceMethods.is_purine(b)) or\
               (SequenceMethods.is_pyrimidine(a) and SequenceMethods.is_pyrimidine(b))

    @staticmethod
    def is_transversion(a, b):
        assert a != b
        return (SequenceMethods.is_purine(a) and SequenceMethods.is_pyrimidine(b)) or\
               (SequenceMethods.is_pyrimidine(a) and SequenceMethods.is_purine(b))


    @staticmethod
    def get_highest_score_sequence(weights):
        '''
        Given a 4 * N sequence one hot encoding vector, return the highest scoring sequence
        :param weights:
        :return:
        '''

        assert len(weights) % 4 == 0
        seq = ""
        for i in range(len(weights) / 4):
            next_weights = weights[i * 4: (i + 1) * 4]
            options = [[w, b] for w, b in zip(next_weights, 'ACGT')]
            # print options
            selected_base = sorted(options, key=lambda x: -x[0])[0][1]
            seq += selected_base
        return seq

    @staticmethod
    def get_single_nucleotide_changes(query_seq):
        return [query_seq[:i] + b + query_seq[i + 1:]
                for i in range(len(query_seq)) for b in 'ACGT']

    @staticmethod
    def get_n_nucleotide_changes(s, d=2):
        N = len(s)
        letters = 'ACGT'
        pool = list(s)

        for indices in combinations(list(range(N)), d):
            for replacements in product(letters, repeat=d):
                skip = False
                for i, a in zip(indices, replacements):
                    if pool[i] == a: skip = True
                if skip: continue

                keys = dict(list(zip(indices, replacements)))
                yield ''.join([pool[i] if i not in indices else keys[i]
                               for i in range(N)])

    @staticmethod
    def is_palindrome(seq):
        return SequenceMethods.get_complementary_seq(seq) == seq

    @staticmethod
    def get_non_redundant_sequence(s, random_seed=None):
        iupac_table = SequenceMethods.get_iupac_table()

        final = ""
        if random_seed is not None:
            seed(random_seed)

        for si in s:
            if len(iupac_table[si]) == 1:
                final += list(iupac_table[si])[0]
            else:
                final += list(iupac_table[si])[randint(0, len(iupac_table[si]) - 1)]

        return final