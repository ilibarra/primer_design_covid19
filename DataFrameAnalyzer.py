

import pandas as pd
from os import listdir
from path_functions import *
from os.path import join, abspath
import numpy as np
import pickle

class DataFrameAnalyzer:
    @staticmethod
    def to_tsv_gz(df, path, sep='\t', index=None, **kwargs):
        df.to_csv(path, sep=sep, index=index, compression='gzip', **kwargs)
        print('dataframe saved at...')
        print(path)
        print((abspath(path)))

    @staticmethod
    def to_pickle(df, path):
        pickle.dump(df, open(path, 'wb'))
        print('pickle saved at...')
        print(path)
        print((abspath(path)))

    @staticmethod
    def split(dfm, chunk_size):
        def index_marks(nrows, chunk_size):
            return list(range(1 * chunk_size, (nrows // chunk_size + 1) * chunk_size, chunk_size))
        indices = index_marks(dfm.shape[0], chunk_size)
        return np.split(dfm, indices)

    @staticmethod
    def read_pickle(path, log=False, **kwargs):
        if log:
            print(('loading pickle path', path))
        return pickle.load(open(path, 'rb'), **kwargs)

    @staticmethod
    def to_tsv(df, path, sep='\t', index=None, log=True, **kwargs):
        df.to_csv(path, sep=sep, index=index, **kwargs)
        if log:
            print('dataframe saved at (path/abspath)...')
            print(path)
            print((abspath(path)))

    @staticmethod
    def read_tsv_gz(path, engine=None, header='infer', sep='\t', index_col=None, columns=None, log=False,
                    **kwargs):

        if log:
            print(('loading file', basename(path)))
        df = pd.read_csv(path, sep=sep, header=header,
                         compression='gzip', engine=engine, index_col=index_col,
                         **kwargs)
        if columns is not None:
            df.columns = columns
        if log:
            print(('done...', basename(path)))
        return df

    @staticmethod
    def read_tsv(path, engine=None, index_col=None, header='infer', columns=None, sep='\t', **kwargs):
        try:
            df = pd.read_csv(path, sep=sep, index_col=index_col, header=header,
                             engine=engine, **kwargs)
            if columns is not None:
                df.columns = columns
            return df
        except pd.errors.EmptyDataError:
            print('dataframe file was empty...')
            return pd.DataFrame(columns=columns)


    @staticmethod
    def read_mwords_tsv(path, engine=None, index_col=None, header=None, sep='\t', **kwargs):
        df = pd.read_csv(path, sep=sep, index_col=index_col, header=header,
                           engine=engine, **kwargs)
        df.columns = ['seq', 'rel.affinity', 'counts']
        return df

    @staticmethod
    def pkldump(obj, path):
        pickle.dump(obj, open(path, 'wb'))

    @staticmethod
    def pklload(path):
        return pickle.load(open(path, 'rb'))

    @staticmethod
    def norm_01(series):
        return (series - np.min(series)) / (np.max(series) - np.min(series))
    @staticmethod
    def read_multiple_tsv_gz(directory_path, stopat=None, query=None, column_filter=None, column_filter_thr=None, **kwargs):
        res = []
        counter = 0

        for tsv_gz_filename in listdir(directory_path):
            if not tsv_gz_filename.endswith('.tsv.gz') and not tsv_gz_filename.endswith('.gz'):
                continue
            if query is not None and not query in tsv_gz_filename:
                continue
            counter += 1
            if counter % 150 == 0:
                print(('files read so far: %i' % counter))
            p = join(directory_path, tsv_gz_filename)
            if filesize(p) == 0:
                print(('empty dataframe. Skip.', p))
                continue
            try:
                df = DataFrameAnalyzer.read_tsv_gz(p, **kwargs)
                df['filename'] = tsv_gz_filename.replace(".tsv.gz", '')
                if kwargs.get('columns', None):
                    df = df[kwargs.get('columns')]

                # filter
                if column_filter is not None:
                    assert column_filter_thr is not None
                    df = df[df[column_filter] >= column_filter_thr].reset_index(drop=True)
                res.append(df)
            except IOError:
                print('IOError when reading path... ')
                print(p)
            finally:
                if stopat is not None and len(res) >= stopat:
                    break

        # print len(res)
        if len(res) != 0:
            return pd.concat(res).reset_index(drop=True)
        else:
            return None

    @staticmethod
    def read_multiple_tsv(directory_path, stopat=None, keyword_filename=None, **kwargs):
        res = []

        for i, tsv_filename in enumerate(listdir(directory_path)):
            if keyword_filename is not None and not keyword_filename in tsv_filename:
                continue
            # print tsv_filename
            df = DataFrameAnalyzer.read_tsv(join(directory_path, tsv_filename), **kwargs)
            df['filename'] = tsv_filename.replace(".tsv", '')
            res.append(df)
            if stopat != None and len(res) >= stopat:
                break
            # print df.head()
        return pd.concat(res).reset_index(drop=True)

    @staticmethod
    def get_dict(df, a, b):
        return pd.Series(df[b].values if b is not None else df.index,
                         index=df[a].values if a is not None else df.index).to_dict()

    @staticmethod
    def write_list(g, output_path):
        writer = open(output_path, 'w')
        for g in list(g):
            writer.write(g + "\n")
        writer.close()
        print('output written to...')
        print((abspath(output_path)))
    @staticmethod
    def dataframe_to_matrix(df, not_found_value=None, sep='_and_'):
        '''
        Given a dataframe with three entries [key1, key2, VALUE] return a key1 x key2 matrix
        :param df:
        :return:
        '''
        k1, k2, v = df.columns
        print((k1, k2, v))
        df['tmp.k'] = df[k1] + sep + df[k2]
        print((df[v].values))
        value_by_key = DataFrameAnalyzer.get_dict(df, 'tmp.k', v)
        uniq_k1 = list(set(df[k1]))
        uniq_k2 = list(set(df[k2]))
        matrix = [[None for j in range(len(uniq_k2))] for i in range(len(uniq_k1))]
        n_found = 0
        n_not_found = 0
        for i, next_k1 in enumerate(uniq_k1):
            for j, next_k2 in enumerate(uniq_k2):
                k1 = next_k1 + sep + next_k2
                k2 = next_k2 + sep + next_k1

                # print i, j, 'wgEncodeAwgTfbsSydhImr90CebpbIggrabUniPk.fa_JUND_NRTGACTCATN_me' in value_by_key
                if k1 in value_by_key:
                    matrix[i][j] = value_by_key[k1]
                    n_found += 1
                elif k2 in value_by_key:
                    matrix[i][j] = value_by_key[k2]
                    n_found += 1
                else:
                    matrix[i][j] = not_found_value
                    n_not_found += 1
        return pd.DataFrame(matrix, index=uniq_k1, columns=uniq_k2)
