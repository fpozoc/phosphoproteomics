#!/usr/bin/env python
# -*- coding: utf-8 -*-

import re

import pandas as pd 
from Bio import SeqIO
from matplotlib import colors
import matplotlib.pyplot as plt
import numpy as np

import plotly.graph_objects as go
from plotly.subplots import make_subplots


def background_gradient(s, m, M, cmap='RdBu', low=0, high=0):
    '''
    Background gradient color function.
    '''
    rng = M - m
    norm = colors.Normalize(m - (rng * low),
                            M + (rng * high))
    normed = norm(s.values)
    c = [colors.rgb2hex(x) for x in plt.cm.get_cmap(cmap)(normed)]
    return ['background-color: %s' % color for color in c]


def calculate_structure_stats(df):
    df = df[~df['ss3'].isnull()]
    df =  df.sort_values(by=['proteinName', 'proteinSeq-loc'], ascending=[False, True]).drop_duplicates(
        ['proteinName', 'proteinSeq-loc']).reset_index(drop=True)
    df = df[['proteinName', 'proteinSeq-loc', 'res-ptm', 'proteinSeq', 'proteinLen', 'rsa', 'ss3', 'ss8', 'dis']]
    df = df.groupby(['proteinName', 'proteinSeq', 'proteinLen', 'rsa', 'ss3', 'ss8', 'dis']).agg(list).reset_index()
    df = df.loc[df['proteinLen'] == (df['rsa'].str.len() & df['ss3'].str.len() & df['dis'].str.len())]
    
    df['ph_rsa'] = [list(map(replaceitem,[int(k[i-1]) for i in ploc])) for ploc, k  in zip(df['proteinSeq-loc'], df['rsa'])]
    df['ph_rsa'] = [sum(map(int, k))/len(k) for k in df['ph_rsa']]
    df['nph_rsa'] = [remove(k, ploc) for ploc, k  in zip(df['proteinSeq-loc'], df['rsa'])]
    df['nph_rsa'] = [sum(list(map(replaceitem, list(map(int, k)))))/len(k) for k in df['nph_rsa'].str.split('').str[1:-1]]

    df['ph_dis'] = [list(map(replaceitem, [int(k[i-1]) for i in ploc])) for ploc, k  in zip(df['proteinSeq-loc'], df['dis'])]
    df['ph_dis'] = [sum(map(int, k))/len(k) for k in df['ph_dis']]
    df['nph_dis'] = [remove(k, ploc) for ploc, k  in zip(df['proteinSeq-loc'], df['dis'])]
    df['nph_dis'] = [sum(list(map(replaceitem, list(map(int, k)))))/len(k) for k in df['nph_dis'].str.split('').str[1:-1]]
    
    df['ph_ss3'] = [''.join([k[i-1] for i in ploc]) for ploc, k in zip(df['proteinSeq-loc'], df['ss3'])]
    df['ph_ss3_C'] = df['ph_ss3'].str.count('C')
    df['ph_ss3_H'] = df['ph_ss3'].str.count('H')
    df['ph_ss3_E'] = df['ph_ss3'].str.count('E')

    df['nph_ss3'] = [remove(k, ploc) for ploc, k  in zip(df['proteinSeq-loc'], df['ss3'])]
    df['nph_ss3_C'] = df['nph_ss3'].str.count('C')
    df['nph_ss3_H'] = df['nph_ss3'].str.count('H')
    df['nph_ss3_E'] = df['nph_ss3'].str.count('E')
    df = df[['ph_rsa', 'nph_rsa', 'ph_dis', 'nph_dis', 
            'ph_ss3_C', 'ph_ss3_H', 'ph_ss3_E', 'nph_ss3_C',
            'nph_ss3_H', 'nph_ss3_E']]
#     return round(df.sum(),2)
    return df


def disorder_structure(df_disprot, df_fastas):
    '''
    Function to map the disorder structure (available in disprot human structures) with uniprot aa sequence. 
    '''
    df = pd.merge(df_disprot, df_fastas[['proteinName', 'proteinLen', 'proteinSeq']], on='proteinName', how='left')
    df = df[~df['proteinSeq'].isnull()]
    df['disorder'] = [re.sub(rseq, '#'*(re.search(rseq, seq).end()-re.search(rseq, seq).start()), seq) if re.search(rseq, seq) else re.sub(seq[start-1:end-1], '#'*(end-start), seq) for seq, rseq, start, end in zip(df['proteinSeq'], df['region_sequence'], df['start'], df['end'])]
    df['disorder'] = [re.sub('[^#]', '-', disseq) for disseq in df['disorder']]
    return df[['proteinName', 'disorder']]


def disorder_counter(df, name):
    '''
    Counting ordered and disordered regions
    '''
    col = 'kmer-dis'
    df_counts = pd.DataFrame()
    df_counts['ordered'] = df[col].str.count('-')
    df_counts['disordered'] = df[col].str.count('#')

    counts= pd.DataFrame({
        'ordered': [df_counts['ordered'].sum()],
        'disordered': [df_counts['disordered'].sum()],
    }).T.rename(columns={0: f'{name}_n'})
    counts[f'{name}_%'] = counts[f'{name}_n']/counts[f'{name}_n'].sum()
    return counts
    

def generate_motifs(df, column, aa):
    '''
    Generating 13mers from pandas DataFrame with motifs.
    '''
    df_mot1 = pd.DataFrame()
    df_mot1['kmer-str'] = [[f"{pname}_{x.start()+1-6}_{x.start()+1+7}_{seq[x.start()-6:x.start()+7]}" for x in re.finditer(aa, seq)] for pname, seq in zip(df.proteinName, df.proteinSeq)]
    df_mot2 = df_mot1['kmer-str'].explode().reset_index()
    df_mot2 = df_mot2.loc[~df_mot2['kmer-str'].isnull()]
    df_mot2['proteinName'] = df_mot2['kmer-str'].str.split('_').str[0]
    df_mot2['kmer_start'] = df_mot2['kmer-str'].str.split('_').str[1].astype(int)
    df_mot2['kmer_end'] = df_mot2['kmer-str'].str.split('_').str[2].astype(int)
    df_mot2['kmer'] = df_mot2['kmer-str'].str.split('_').str[3]
    df_mot2 = df_mot2[df_mot2['kmer'].str.len() == 13].reset_index()
    df_mot2['res'] = aa
    return df_mot2[['proteinName', 'res', 'kmer_start', 'kmer_end', 'kmer']]


def get_percentages(df, axis_count='column'):
    '''
    Calculating percentages per row/column for a pandas DataFrame values.
    '''
    if axis_count == 'row':
        df = df[df.columns].div(df[df.columns].sum(axis=1), axis=0)
    elif axis_count == 'column':
        df = df[df.columns].div(df[df.columns].sum(axis=0), axis=1)
    return df


def load_fasta(filepath):
    """
    Fasta Loader

    It takes a fasta file from GENCODE annotation and returns a DataFrame.
    It uses SeqIO module from BioPython (https://biopython.org/wiki/SeqIO) to parse an standard fasta
    file and retrive id and seq and parse info from the id.

    Parameters
    ----------
    file: text file
        fasta file downloaded from source

    Returns
    -------
    df: pandas DataFrame
        DataFrame which contains fasta sequences, lengh and id from annotation source.
    """
    seqlist = list()
    with open(filepath, 'r') as handle:
        for record in SeqIO.parse(handle, "fasta"):
            seqlist.append({'proteinName': record.name.split('|')[1], # problems with synonyms
                            'proteinLen': int(len(str(record.seq))), 
                            'proteinSeq': str(record.seq), 
                            'id': str(record.id)})             
    df = pd.DataFrame(seqlist)
    print(f'Fasta shape = {df.shape}')            
    return df


def pandas_vis(df):
    '''
    Visualize a dataframe with gradient in their minimum and maximum values.
    '''
    df = df.style.apply(background_gradient,
               cmap='RdBu',
               m=-1.5,
               M=1.5,
               low=0,
               high=0.2)
    return df


def merge_dfs(df_motifs, df_exptypes, df_structural):
    df = pd.merge(df_motifs, df_structural, on='proteinName', how='left')
    df['kmer-ss'] = [ss[start-1:end-1] if ss is not np.nan else np.nan for start, end, ss in zip(df.kmer_start, df.kmer_end, df.secondary_structure)]
    df['kmer-dis'] = [dis[start-1:end-1] if dis is not np.nan else np.nan for start, end, dis in zip(df.kmer_start, df.kmer_end, df.disorder)]
    df = pd.merge(df, df_exptypes[['proteinName', 'kmer', 'type']], how='left', on=['proteinName', 'kmer'])
    df = df.drop_duplicates().reset_index(drop=True)
    return df 


def parse_ss_features(netsurp_filepath:str)->list:
    df = pd.read_csv(netsurp_filepath)
    df[df.select_dtypes(include='float').columns] = round(df.select_dtypes(include='float')*10,0).astype(int).astype(str).replace('10','9')
    df['p_ss3'] = tuple(zip(df['p[q3_H]'], df['p[q3_E]'], df['p[q3_C]']))
    df['p_ss3'] = df['p_ss3'].apply(max)
    df['p_ss8'] = tuple(zip(df['p[q8_G]'], df['p[q8_H]'], 
                            df['p[q8_I]'], df['p[q8_B]'], 
                            df['p[q8_E]'], df['p[q8_S]'],
                            df['p[q8_T]'], df['p[q8_B]']))
    df['p_ss8'] = df['p_ss8'].apply(max)
    ss_dict = dict()
    ss_dict['id'] = df['id'].values[0]
    ss_dict['seq'] = df['seq'].sum()
    ss_dict['len'] = len(df['seq'].sum())
    ss_dict['rsa'] = df['rsa'].sum()
    ss_dict['ss3'] = df['q3'].sum()
    ss_dict['p_ss3'] = df['p_ss3'].sum()
    ss_dict['ss8'] = df['q8'].sum()
    ss_dict['p_ss8'] = df['p_ss8'].sum()
    ss_dict['dis'] = df['disorder'].sum()
    ss_df = pd.DataFrame(ss_dict, index=[0])
    return ss_df


def plot_stacked_data(df:list, title:str, cmp:str='Paired'):
    stacked_data = df.apply(lambda x: x*100/sum(x), axis=1)
    stacked_data.plot.bar(stacked=True, rot=0,  fontsize=12, figsize=(8,6), colormap=cmp)
    plt.title(title)
    plt.ylabel("(%)")


def process_peptides(df, seqs, filetype):
    '''
    Processing peptides file and getting column of interest.
    Exploding ptms and Location probabilities is an important step here.
    '''
    df['peptide'] = df['peptide'].astype(str)
    df['proteinName'] = df['proteinName'].astype(str)
    df['geneName'] = df['geneName'].astype(str)
    df['ptm'] = df['ptm'].astype(str)
    

    if filetype == 'wt':
        df = df[['peptide', 'ac', 'geneName', 'ptm', 'Localization Probability', 'FTestQValue']]
        df = df.rename(columns={'Localization Probability': 'p-loc', 'FTestQValue': 'q-value', 'ac': 'proteinName'})
        # converting to proper types
        df['q-value'] = df['q-value'].astype(float)
        df['p-loc'] = df['p-loc'].astype(str)
        # exploding ptm and location probability
        s1 = df['ptm'].str.split('|', expand=True).stack().str.strip().reset_index(level=1, drop=True)
        s2 = df['p-loc'].str.split(',', expand=True).stack().str.strip().reset_index(level=1, drop=True)
        df1 = pd.merge(s1.reset_index(), s2.reset_index(), how='left', on='index').set_index('index').rename(
            columns={'0_x': 'ptm', '0_y': 'p-loc'})
        df = df.drop(['ptm','p-loc'], axis=1).join(df1).reset_index(drop=True)
        # parsing ptms and location probabilities
        df['p-loc'] = [float(col.split(' %')[0])/100 if '%' in col else float(col) for col in df['p-loc']]

    elif filetype == 'ko':
        df = df[['peptide', 'ac', 'geneName', 'ptm', 'FTestQValue']]
        df = df.rename(columns={'FTestQValue': 'q-value', 'ac': 'proteinName'})
        df['q-value'] = df['q-value'].astype(float)
        df = df.drop(['ptm'], axis=1).join(df['ptm'].str.split('|').explode()).reset_index(drop=True)

    df = pd.merge(df, seqs, how='left', on='proteinName')
    df['n-ptm'] = df['ptm'].str.split('[').str[1].str.split(']').str[0].astype(int)
    df['type-ptm'] = df['ptm'].str.split('] ').str[1].str.split('(').str[0]
    df['aa-ptm'] = df['ptm'].str.split('(').str[1].str.split(')').str[0]
    df['res-ptm'] = [pep[n-1] for pep, n in zip(df['peptide'], df['n-ptm'])]
    df['proteinSeq-loc'] = [str(seq).find(pep)+n for n, pep, seq in zip(df['n-ptm'], df['peptide'], df['proteinSeq'])]
    df = df.loc[~df['proteinSeq'].isnull()] # 'P08107' 'Q8N655' were removed (obsolete).
    df['kmer'] = [pep[n-1-6:n-1+7] if len(pep[0:n]) > 6 and len(pep[n-1:]) > 6 else seq[nseq-1-6:nseq-1+7] if nseq > 6 else '-' for n, pep, seq, nseq in zip(df['n-ptm'], df['peptide'], df['proteinSeq'], df['proteinSeq-loc'])]
    df['kmerLen'] = df['kmer'].str.len()     
    
    # filters 
    df = df[df['type-ptm'].str.contains('Phospho')] # filter 1 - only phosphorilated motifs and non-oxidation
    df = df.loc[df['q-value'] < 0.05] ## filter 2 - q value 
    df['type'] = filetype
    df = df.drop_duplicates().reset_index(drop=True) # there are duplicates. It removes duplications for the whole row
    return df


def remove(s, indx):
    return ''.join(x for n, x in enumerate(s) if n not in indx)


def replaceitem(x):
    if x > 5:
        return 1
    else:
        return 0


def replace_list_brackets(df):
    for feature in df.columns.values:
        df[f'{feature}'] = df[f'{feature}'].apply(lambda x: str(x).replace("['","").replace("']","").replace("""'""", ""))
    return df

def secondary_structure(df_uniprot_pdb, df_fastas):
    '''
    Function to map the secondary structure (available in pdb human structures) with uniprot aa sequence. 
    '''
    df = pd.merge(df_uniprot_pdb, df_fastas[['proteinName', 'proteinLen']], on='proteinName', how='left')
    df = df[~df['proteinLen'].isnull()]
    df.loc[df['ss'] == 'Helix', 'ss'] = 'H'
    df.loc[df['ss'] == 'Beta strand', 'ss'] = 'B'
    df.loc[df['ss'] == 'Turn', 'ss'] = 'T'
    df = df.sort_values(by=['proteinName', 'start'], ascending=True).reset_index(drop=True)
    df['spaces'] = df.groupby('proteinName').apply(lambda x: x['start'].sub(x['end'].shift())).reset_index().sort_values('level_1').reset_index(drop=True)[0]
    df.loc[df['spaces'].isnull(), 'spaces'] = df['start']
    df['spaces'] = df['spaces'].astype(int)
    df['ss-str'] = ['-'*sp + ss*(end-start) for ss, start, end, sp in zip(df['ss'], df['start'], df['end'], df['spaces'])]

    dfg = df.groupby('proteinName').apply(lambda x: x['ss-str'].sum()).reset_index().rename(columns={0:'secondary_structure'})
    dfg['ss-end'] = (df.groupby('proteinName').last()['proteinLen'] - df.groupby('proteinName').last()['end']).reset_index(drop=True).astype(int)
    dfg['secondary_structure'] = [ss+ssend*'-' for ss, ssend in zip(dfg['secondary_structure'], dfg['ss-end'])]
    return dfg.drop('ss-end', axis=1)


def secondary_structure_counter(df, name):
    '''
    Counting helix, beta strands and turns
    '''
    col = 'kmer-ss'
    df_counts = pd.DataFrame()
    df_counts['helix'] = df[col].str.count('H')
    df_counts['strand'] = df[col].str.count('B')
    df_counts['turn'] = df[col].str.count('T')

    counts= pd.DataFrame({
        'helix': [df_counts['helix'].sum()],
        'beta-strand': [df_counts['strand'].sum()],
        'turn': [df_counts['turn'].sum()]
    }).T.rename(columns={0: f'{name}_n'})
    counts[f'{name}_%'] = counts[f'{name}_n']/counts[f'{name}_n'].sum()
    return counts


def string_counts(df, col):
    '''
    Counting character occurence per DataFrame with same length strings.
    '''
    df_sum = pd.DataFrame()
    for i in range(0, 13):
        df_sum[i] = df[col].str[i]
    return df_sum.apply(pd.Series.value_counts)

## visualizations

def df_to_plotly(df):
    return {'z': df.values.tolist(),
            'x': [i+1 for i in df.columns.tolist()],
            'y': df.index.tolist(),
            'zmin': -1.2, 
            'zmax': 1.2,
            'xgap': 0.1,
            'ygap': 0.1,
            'colorscale': 'RdBu',
            'opacity': 1,
            'colorbar': dict(
                title='log10',
                dtick=0.2,
                tickvals=np.linspace(-1,1,11),
                ticks='outside',
                thickness=10)
           }


def single_heatmap(df):
    return go.Heatmap(df_to_plotly(df))
    

def plot_heatmaps(df1, df2, df3, title):
    fig = make_subplots(rows=1, cols=3, horizontal_spacing=0.05,
                        subplot_titles=("Serine", "Threonine", "Tyrosine"))
    fig.append_trace(single_heatmap(df1), 1, 1)
    fig.append_trace(single_heatmap(df2), 1, 2)
    fig.append_trace(single_heatmap(df3), 1, 3)
    for row, col in enumerate(range(1,4)):
        fig.update_xaxes(title_text="", tickmode='linear', row=1, col=col)
    fig.update_layout(
        title=title,
        height=500,
        width=850,
        font=dict(
            family="Arial",
            size=11,
            color="black"
        )
    )
    fig.show()
