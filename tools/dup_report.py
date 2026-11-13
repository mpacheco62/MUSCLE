#!/usr/bin/env python3
# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright (C) 2025 Matias Pacheco-Alarcon <matias.pacheco.a@gmail.com>
"""Repeated blocks of source (analysis of duplication, docs/fypp/ANALISIS_REPETICION.md).

    tools/dup_report.py <dir[,dir...]> [N]

Lists the N largest blocks of >= 6 normalized lines that appear more than once (comments, case,
spacing, numbers and the kind names sp/dp are normalized, so unrolled code that differs only in
indices counts as repeated: an upper bound of the textual copies).
"""
import collections
import glob
import re
import sys
W = 6
roots=sys.argv[1].split(','); topn=int(sys.argv[2])
files=[f for r in roots for f in glob.glob(r+'/**/*.fypp',recursive=True)+glob.glob(r+'/**/*.f90',recursive=True)]
def norm(l):
    if l.lstrip().startswith('#!'): return None
    if not l.lstrip().startswith('#'): l=re.sub(r'!.*$','',l)
    l=re.sub(r'\s+',' ',l).strip()
    if not l or l.startswith('#line'): return None
    l=re.sub(r'\b(sp|dp)\b','K',l); l=re.sub(r'\d+(\.\d+)?([eE][-+]?\d+)?(_\w+)?','N',l)
    return l.lower()
seqs={f:[(i,n) for i,l in enumerate(open(f,errors='ignore'),1) if (n:=norm(l))] for f in files}
idx=collections.defaultdict(list)
for f,ls in seqs.items():
    for j in range(len(ls)-W+1):
        key=tuple(x[1] for x in ls[j:j+W])
        if len(set(key))>=3: idx[key].append((f,j))
# merge windows into maximal runs keyed by the tuple of (file) occurrences shifted together
used=set(); clusters=[]
for f,ls in seqs.items():
    j=0
    while j<=len(ls)-W:
        key=tuple(x[1] for x in ls[j:j+W]); occ=idx.get(key,[])
        if len(occ)>1 and (f,j) not in used:
            others=[o for o in occ if o!=(f,j)]
            L=W
            while True:
                ok=j+L<len(ls) and all(o[1]+L<len(seqs[o[0]]) and seqs[o[0]][o[1]+L][1]==ls[j+L][1] for o in others)
                if not ok: break
                L+=1
            for o in occ:
                for k in range(L): used.add((o[0],o[1]+k))
            clusters.append((L*len(occ),L,sorted({(o[0],seqs[o[0]][o[1]][0]) for o in occ})))
            j+=L
        else: j+=1
clusters.sort(reverse=True)
intra=sum(c[0]-c[1] for c in clusters if len({a for a,_ in c[2]})==1)
inter=sum(c[0]-c[1] for c in clusters if len({a for a,_ in c[2]})>1)
print(f"redundant lines (extra copies): within a file {intra}, across files {inter}")
for tot,L,occ in clusters[:topn]:
    fs=collections.Counter(a.split('/')[-1] for a,_ in occ)
    print(f"{L:4d} lines x{len(occ)} ({tot-L} redundant)  "+", ".join(f"{k}{'' if v==1 else '*%d'%v}" for k,v in list(fs.items())[:5]))
