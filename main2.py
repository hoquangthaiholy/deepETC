"""
" Dataset preprocessing
" Handling protein sequence data.
"   Split data into train/test set
"""

from sklearn.model_selection import train_test_split
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
    df = pd.DataFrame(columns=['ID', 'SEQUENCE', 'LENGTH', 'COMPLEX' , 'LABEL'])
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

        matches = re.findall(r'^DR   GO; GO:(\S*);', data, flags=re.M)
        comp = ''

        for match in matches:
            if match in ['0016651','0016655','0003954','0050136','0008137']:
                comp = 'I'
                break
            if match in ['0000104','0008177']:
                comp = 'II'
                break
            if match in ['0008121','0045153']:
                comp = 'III'
                break
            if match in ['0004129']:
                comp = 'IV'
                break
            if match in ['0046961']:
                comp = 'V'
                break

        matches = re.findall(r'^KW   (.*)', data, flags=re.M)
        label = 'TP'

        for match in matches:
            if match.find('Electron transport') != -1:
                label = 'ET'
                break
        
        # If a protein belongs to 5 complexes, it will be ET
        if comp != '':
            label = 'ET'

        df.loc[n] = [fid, seq, len(seq), comp, label]
    return df


def dataframe_to_fasta(df, output_folder='fasta'):
    for i in range(len(df)):
        fout = open(f"{output_folder}/{df.iloc[i,:].ID}.fasta", 'w')
        fout.write(f"{df.iloc[i,:].SEQUENCE}")
    fout.close()


def dataframe_to_fasta_all(df, output_file='fasta.txt'):
    fout = open(output_file, 'w')
    for i in range(len(df)):
        fout.write(
            f">sp|{df.iloc[i,:].ID}|length={df.iloc[i,:].LENGTH}|label={df.iloc[i,:].LABEL}\n{df.iloc[i,:].SEQUENCE}\n")
    fout.close()

# df = uniprot_to_dataframe('uniprot/u niprot.txt')
df = uniprot_to_dataframe('uniprot/uniprot.txt', 'uniprot/uniprot.30.out')
# dataframe_to_fasta(df)

# Check max seqlength
maxseq = 0
for i in range(len(df)):
    if len(df.iloc[i, :].SEQUENCE) > maxseq:
        maxseq = len(df.iloc[i, :].SEQUENCE)

# Split validation/test
X_train, X_test, y_train, y_test = train_test_split(
    df[['SEQUENCE']], df[['LABEL']], test_size=0.167, random_state=42)

from collections import Counter
print(Counter(df.LABEL))
print(Counter(y_train.LABEL))
print(Counter(y_test.LABEL))

"""""
" Export train/test dataset for Deepfam
"
"""
CHARSET = {'A': 0, 'C': 1, 'D': 2, 'E': 3, 'F': 4, 'G': 5, 'H': 6,
           'I': 7, 'K': 8, 'L': 9, 'M': 10, 'N': 11, 'P': 12, 'Q': 13,
           'R': 14, 'S': 15, 'T': 16, 'V': 17, 'W': 18, 'Y': 19, 'X': 20,
           'O': 20, 'U': 20,
           'B': (2, 11),
           'Z': (3, 13),
           'J': (7, 9)}
CHARLEN = 21


def encoding_seq_np(seq):
    arr = [0] * CHARLEN
    for i, c in enumerate(seq):
        if c == "_":
            # let them zero
            continue
        elif isinstance(CHARSET[c], int):
            idx = CHARLEN * i + CHARSET[c]
            arr[idx] = 1
        else:
            idx1 = CHARLEN * i + CHARSET[c][0]
            idx2 = CHARLEN * i + CHARSET[c][1]
            arr[idx1] = 0.5
            arr[idx2] = 0.5
        return arr


ftrn = open('bin_train.txt', 'w')
for i in X_train.index:
    arr = []
    arr.append(1 if df.iloc[i, :].LABEL == 'ET' else 0)
    seq = list(df.iloc[i, :].SEQUENCE)
    for j in range(maxseq):
        if j < len(seq):
            arr.extend(encoding_seq_np(seq[j]))
        else:
            arr.extend([0]*CHARLEN)
    ftrn.write(f"{','.join([ str(k) for k in arr])}\n")
ftrn.close()

ftrn = open('bin_test.txt', 'w')
for i in X_test.index:
    arr = []
    arr.append(1 if df.iloc[i, :].LABEL == 'ET' else 0)
    seq = list(df.iloc[i, :].SEQUENCE)
    for j in range(maxseq):
        if j < len(seq):
            arr.extend(encoding_seq_np(seq[j]))
        else:
            arr.extend([0]*CHARLEN)
    ftrn.write(f"{','.join([ str(k) for k in arr])}\n")
ftrn.close()

"""""
" Export PSSM
"
"""

ftrn = open('pssm_train.txt', 'w')
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
    for c in range(j, maxseq):
        arr.extend([0]*20)
    ftrn.write(
        f"{1 if df.iloc[i,:].LABEL=='ET' else 0},{','.join(str(x) for x in arr)}\n")

ftrn = open('pssm_test.txt', 'w')
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
    for c in range(j, maxseq):
        arr.extend([0]*20)
    ftrn.write(
        f"{1 if df.iloc[i,:].LABEL=='ET' else 0},{','.join(str(x) for x in arr)}\n")

"""""
" Export PSFM
"
"""

ftrn = open('psfm_train.txt', 'w')
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
    for c in range(j, maxseq):
        arr.extend([0]*20)
    ftrn.write(
        f"{1 if df.iloc[i,:].LABEL=='ET' else 0},{','.join(str(x) for x in arr)}\n")

ftrn = open('psfm_test.txt', 'w')
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
    for c in range(j, maxseq):
        arr.extend([0]*20)
    ftrn.write(
        f"{1 if df.iloc[i,:].LABEL=='ET' else 0},{','.join(str(x) for x in arr)}\n")


"""""
" Export PSSM + PSFM
"
"""

ftrn = open('pssm_psfm_train2.txt', 'w')
m = 0
for i in X_train.index:
    m += 1
    pssm = open(f"pssm_all/{df.iloc[i,:].ID}.pssm").readlines()[3:-6]
    print(f"\r{m}/{len(X_train)}", end='')
    arr = []
    j = 0
    for line in pssm:
        j += 1
        arr.extend([float(k) for k in line.split()[2:42]])
    for c in range(j, maxseq):
        arr.extend([0]*40)
    ftrn.write(
        f"{1 if df.iloc[i,:].LABEL=='ET' else 0},{','.join(str(x) for x in arr)}\n")

ftrn = open('pssm_psfm_test2.txt', 'w')
m = 0
for i in X_test.index:
    m += 1
    pssm = open(f"pssm_all/{df.iloc[i,:].ID}.pssm").readlines()[3:-6]
    print(f"\r{m}/{len(X_test)}", end='')
    arr = []
    j = 0
    for line in pssm:
        j += 1
        arr.extend([float(k) for k in line.split()[2:42]])
    for c in range(j, maxseq):
        arr.extend([0]*40)
    ftrn.write(
        f"{1 if df.iloc[i,:].LABEL=='ET' else 0},{','.join(str(x) for x in arr)}\n")

"""""
" Export Sum PSSM
"
"""
import math

ftrn = open('pssm_sum_train.csv', 'w')
m = 0
for i in X_train.index:
    m += 1
    arr = {'A': [0]*20, 'R': [0]*20, 'N': [0]*20, 'D': [0]*20, 'C': [0]*20, 'Q': [0]*20, 'E': [0]*20, 'G': [0]*20, 'H': [0]*20, 'I': [0]
       * 20, 'L': [0]*20, 'K': [0]*20, 'M': [0]*20, 'F': [0]*20, 'P': [0]*20, 'S': [0]*20, 'T': [0]*20, 'W': [0]*20, 'Y': [0]*20, 'V': [0]*20}
    arr2 = {'A': 0, 'R': 0, 'N': 0, 'D': 0, 'C': 0, 'Q': 0, 'E': 0, 'G': 0, 'H': 0, 'I': 0
        , 'L': 0, 'K': 0, 'M': 0, 'F': 0, 'P': 0, 'S': 0, 'T': 0, 'W': 0, 'Y': 0, 'V': 0}
    pssm = open(f"pssm_all/{df.iloc[i,:].ID}.pssm").readlines()[3:-6]
    print(f"\r{m}/{len(X_train)}", end='')
    for line in pssm:
        c = line.split()[1]
        if c in ['X','U', 'O']:
            continue
        for k in range(20):
            arr2[c] += 1
            arr[c][k] += float(line.split()[k+2])
    
    for j in arr:
        for k in range(20):
            arr[j][k] /= len(pssm)
            arr[j][k] = 1/(1 + math.exp(-arr[j][k]))
    
    _arr = []
    for j in arr:
        _arr.extend(arr[j])
    ftrn.write(
        f"{1 if df.iloc[i,:].LABEL=='ET' else 0},{','.join(str(x) for x in _arr)}\n")

ftrn = open('pssm_sum_test.csv', 'w')
m = 0
for i in X_test.index:
    m += 1
    arr = {'A': [0]*20, 'R': [0]*20, 'N': [0]*20, 'D': [0]*20, 'C': [0]*20, 'Q': [0]*20, 'E': [0]*20, 'G': [0]*20, 'H': [0]*20, 'I': [0]
       * 20, 'L': [0]*20, 'K': [0]*20, 'M': [0]*20, 'F': [0]*20, 'P': [0]*20, 'S': [0]*20, 'T': [0]*20, 'W': [0]*20, 'Y': [0]*20, 'V': [0]*20}
    arr2 = {'A': 0, 'R': 0, 'N': 0, 'D': 0, 'C': 0, 'Q': 0, 'E': 0, 'G': 0, 'H': 0, 'I': 0
        , 'L': 0, 'K': 0, 'M': 0, 'F': 0, 'P': 0, 'S': 0, 'T': 0, 'W': 0, 'Y': 0, 'V': 0}
    pssm = open(f"pssm_all/{df.iloc[i,:].ID}.pssm").readlines()[3:-6]
    print(f"\r{m}/{len(X_test)}", end='')
    for line in pssm:
        c = line.split()[1]
        if c in ['X','U', 'O']:
            continue
        for k in range(20):
            arr2[c] += 1
            arr[c][k] += float(line.split()[k+2])
    
    for j in arr:
        for k in range(20):
            arr[j][k] /= len(pssm)
            arr[j][k] = 1/(1 + math.exp(-arr[j][k]))
    
    _arr = []
    for j in arr:
        _arr.extend(arr[j])
    ftrn.write(
        f"{1 if df.iloc[i,:].LABEL=='ET' else 0},{','.join(str(x) for x in _arr)}\n")

#####
#
#  Complexes
#
#######
df_et = df.loc[(df.LABEL=='ET') & (df.COMPLEX!='')]
# df_comp = df_et.loc[df_et.COMPLEX!='']
# df_comp = df.loc[(df.COMPLEX!='')]
df_et = df_et.reset_index(drop=True)

# Factorize labels
df_et.COMPLEX = df_et.COMPLEX.map({ 'I': 0, 'II': 1, 'III': 2, 'IV': 3, 'V': 4 })

# Split validation/test
X_train, X_test, y_train, y_test = train_test_split(
    df_et[['SEQUENCE']], df_et[['COMPLEX']], test_size=0.25, random_state=42, stratify=df_et[['COMPLEX']])


from collections import Counter
print(Counter(df_et.COMPLEX))
print(Counter(y_train.COMPLEX))
print(Counter(y_test.COMPLEX))

"""""
" Export PSSM
"
"""

ftrn = open('pssm_comp_train.csv', 'w')
m = 0
for i in X_train.index:
    m += 1
    pssm = open(f"pssm_all/{df_et.iloc[i,:].ID}.pssm").readlines()[3:-6]
    print(f"\r{m}/{len(X_train)}", end='')
    arr = []
    j = 0
    for line in pssm:
        j += 1
        arr.extend([float(k) for k in line.split()[2:22]])
    for c in range(j, maxseq):
        arr.extend([0]*20)
    ftrn.write(
        f"{df_et.iloc[i,:].COMPLEX},{','.join(str(x) for x in arr)}\n")

ftrn = open('pssm_comp_test.csv', 'w')
m = 0
for i in X_test.index:
    m += 1
    pssm = open(f"pssm_all/{df_et.iloc[i,:].ID}.pssm").readlines()[3:-6]
    print(f"\r{m}/{len(X_test)}", end='')
    arr = []
    j = 0
    for line in pssm:
        j += 1
        arr.extend([float(k) for k in line.split()[2:22]])
    for c in range(j, maxseq):
        arr.extend([0]*20)
    ftrn.write(
        f"{df_et.iloc[i,:].COMPLEX},{','.join(str(x) for x in arr)}\n")
