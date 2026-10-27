import sys
import gzip #altrimenti non riesce a leggere i vcf.gz!!
from collections import Counter

input_vcf = sys.argv[1]
output_txt = sys.argv[2]

sv_types = Counter() #conteggio varianti divise per svtype
supp_vecs = Counter() #conteggio concordanza dei tool
#SUPP_VEC=101 significa che quella variante è stata trovata da sniffles e SVIM, ma non da cutesv
# SUPP_VEC=111 significa trovata da tutti

# Apre il file con gzip se finisce per .gz altrimenti usa l'open normale
opener = gzip.open if input_vcf.endswith('.gz') else open

with opener(input_vcf, 'rt') as f:
    for line in f:
        if line.startswith('#'): continue #headder
        
        info = line.split('\t')[7]
        svtype = "UNKNOWN"
        supp_vec = "UNKNOWN"
        svlen = 0
        
        # Cerca SVTYPE e SUPP_VEC nella colonna INFO
        for item in info.split(';'):
            if item.startswith('SVLEN='):
                svlen = abs(int(item.split('=')[1])) #prende roba a destra dell'uguale
            if item.startswith('SVTYPE='):
                svtype = item.split('=')[1] 
            elif item.startswith('SUPP_VEC='):
                supp_vec = item.split('=')[1] 
        if svlen>=50:
            sv_types[svtype] += 1
            if supp_vec != "UNKNOWN":
                supp_vecs[supp_vec] += 1
        else: continue #sv troppo piccola

with open(output_txt, 'w') as out:
    out.write("=== STRUCTURAL VARIANTS STATS ===\n\n")

    out.write("Total number of SVs:\n")
    tot=0
    for sv, count in sv_types.items():
        tot += count
    out.write(f" {tot}\n")
    
    out.write("\n SV types:\n")
    for svt, count in sv_types.items():
        out.write(f"  - {svt}: {count}\n")
        
    out.write("\n Tool support (Order: Sniffles, cuteSV, SVIM):\n")
    for vec, count in sorted(supp_vecs.items(), key=lambda x: x[1], reverse=True):
        out.write(f"  - Vector {vec}: {count} variants\n")