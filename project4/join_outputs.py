#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import pandas as pd
import re

recname1 = re.compile('RecName: Full=(.*?);')
recname2 = re.compile('RecName: Full=(.*)')

def get_annot(s):
    m = recname1.match(s)
    if m:
        return m.group(1)
    m = recname2.match(s)
    if m:
        return m.group(1)


blast = pd.read_csv('rawdata/blast.tsv', header=None, sep='\t', 
                    names=['prot', 'sacc', 'evalue', 'bitscore', 'pident', 'ssciname', 'stitle', 'qcovs'])
psort = pd.read_csv('rawdata/psort.out', header=None, sep='\t', names=['prot', 'psort_loc'])
targetp = pd.read_csv('rawdata/targetp.txt', header=None, sep='\t', names=['prot', 'targetp_loc'])
pfam = pd.read_csv('rawdata/HHMER.out', header=None, sep='  ', names=['prot', 'pfam_dom'])

blast['annot'] = blast['stitle'].apply(get_annot)
blast0 = blast.sort_values(['prot', 'evalue']).drop_duplicates(['prot'])

res = blast0[['prot', 'evalue', 'annot']].merge(pfam, how='outer', on='prot')
res = res.merge(psort, how='outer', on='prot')
res = res.merge(targetp, how='outer', on='prot')
res['id'] = res['prot'].apply(lambda x: int(x[1:-3]))
res.sort_values(by='id', inplace=True)
#res.fillna(value='', inplace=True)
res.drop(columns='id', inplace=True)

res.to_csv('final.csv', index=False, sep='\t',
           header=['Protein', 'Blast e-value', 'Blast annotation', 'Pfam domains', 'PSORT localization', 'TargetP localization'])