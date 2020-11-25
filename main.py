"""
" Dataset preprocessing
" Handling protein sequence data.
"   Split data into train/test set
"""

import os
import re
import pandas as pd


def uniprot_to_dataframe(file_path, clust_file=''):
    inc_list = []
    if not clust_file == '':
        fin = open(clust_file)
        for line in fin:
            inc_list.append(line.split(' ')[0])
    n = -1
    df = pd.DataFrame(columns=['ID', 'SEQUENCE', 'LENGTH', 'LABEL'])
    all_data = re.split(
        r'^\/\/', ''.join(open(file_path).readlines()), flags=re.M)
    for data in all_data[:-1]:
        matches = re.findall(r'^AC   (\w+);', data, flags=re.M)
        fid = matches[0]

        if (not len(inc_list) == 0 and fid not in inc_list):
            continue

        n += 1
        matches = re.split(r'(^SQ   .*)', data, flags=re.M)
        seq = ''.join(matches[2].split())

        matches = re.findall(r'^KW   (.*)', data, flags=re.M)
        label = 'TP'
        for match in matches:
            if match.find('Electron transport') != -1:
                label = 'ET'
                break
        df.loc[n] = [fid, seq, len(seq), label]
    return df

def dataframe_to_fasta(df, output_folder='fasta'):
    for i in range(len(df)):
        fout = open(f"{output_folder}/{df.iloc[i,:].ID}.fasta",'w')
        fout.write(f"{df.iloc[i,:].SEQUENCE}")
    fout.close()

def dataframe_to_fasta_all(df, output_file='fasta.txt'):
    fout = open(output_file,'w')
    for i in range(len(df)):
        fout.write(f">sp|{df.iloc[i,:].ID}|length={df.iloc[i,:].LENGTH}|label={df.iloc[i,:].LABEL}\n{df.iloc[i,:].SEQUENCE}\n")
    fout.close()


df = uniprot_to_dataframe('uniprot/uniprot.txt', 'uniprot/uniprot.30.out')
dataframe_to_fasta(df)

# Split validation/test
from sklearn.model_selection import train_test_split
X_train, X_test, y_train, y_test = train_test_split(df[['SEQUENCE']],df[['LABEL']], test_size=0.167, random_state=42)

"""""
" Export train/test dataset for Deepfam
"
"""
# Check max seqlength
maxseq = 0
for i in range(len(df)):
    if len(df.iloc[i,:].SEQUENCE)>maxseq:
        maxseq = len(df.iloc[i,:].SEQUENCE)

ftrn = open('raw_train.txt','w')
for i in X_train.index:
    ftrn.write(f"{1 if df.iloc[i,:].LABEL=='ET' else 0}\t{df.iloc[i,:].SEQUENCE.ljust(maxseq,'_')}\n")
ftrn.close()

ftrn = open('raw_test.txt','w')
for i in y_train.index:
    ftrn.write(f"{1 if df.iloc[i,:].LABEL=='ET' else 0}\t{df.iloc[i,:].SEQUENCE.ljust(maxseq,'_')}\n")
ftrn.close()


"""""
" Export PSSM
"
"""

ftrn = open('pssm_train.txt','w')
m = 0
for i in X_train.index:
    m += 1
    pssm = open(f"pssm_all/{df.iloc[i,:].ID}.pssm").readlines()[3:-6]
    print(f"\r{m}/{len(X_train)}", end='')
    arr = []
    j = 0
    for line in pssm:
        j += 1
        arr.extend([float(k) for k in line.split()[2:22]])
    for c in range(j,maxseq):
        arr.extend([0]*20)
    ftrn.write(f"{1 if df.iloc[i,:].LABEL=='ET' else 0}\t{','.join(str(x) for x in arr)}\n")

ftrn = open('pssm_test.txt','w')
m = 0
for i in X_test.index:
    m += 1
    pssm = open(f"pssm_all/{df.iloc[i,:].ID}.pssm").readlines()[3:-6]
    print(f"\r{m}/{len(X_test)}", end='')
    arr = []
    j = 0
    for line in pssm:
        j += 1
        arr.extend([float(k) for k in line.split()[2:22]])
    for c in range(j,maxseq):
        arr.extend([0]*20)
    ftrn.write(f"{1 if df.iloc[i,:].LABEL=='ET' else 0}\t{','.join(str(x) for x in arr)}\n")

"""""
" Export PSFM
"
"""

ftrn = open('psfm_train.txt','w')
m = 0
for i in X_train.index:
    m += 1
    pssm = open(f"pssm_all/{df.iloc[i,:].ID}.pssm").readlines()[3:-6]
    print(f"\r{m}/{len(X_train)}", end='')
    arr = []
    j = 0
    for line in pssm:
        j += 1
        arr.extend([float(k) for k in line.split()[22:42]])
    for c in range(j,maxseq):
        arr.extend([0]*20)
    ftrn.write(f"{1 if df.iloc[i,:].LABEL=='ET' else 0}\t{','.join(str(x) for x in arr)}\n")

ftrn = open('psfm_test.txt','w')
m = 0
for i in X_test.index:
    m += 1
    pssm = open(f"pssm_all/{df.iloc[i,:].ID}.pssm").readlines()[3:-6]
    print(f"\r{m}/{len(X_test)}", end='')
    arr = []
    j = 0
    for line in pssm:
        j += 1
        arr.extend([float(k) for k in line.split()[22:42]])
    for c in range(j,maxseq):
        arr.extend([0]*20)
    ftrn.write(f"{1 if df.iloc[i,:].LABEL=='ET' else 0}\t{','.join(str(x) for x in arr)}\n")