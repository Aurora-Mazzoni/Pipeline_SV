import sys

# ===================== CONFIGURAZIONE =====================
# FILTER: "keep" mantiene il FILTER originale di ogni tool (dichiarando nell'header i filtri mancanti),
#         "pass" mette sempre PASS come nella versione precedente
FILTER_MODE = "keep"

# chiavi INFO specifiche dei singoli tool da eliminare (possono creare conflitti di header)
INFO_TO_REMOVE = {'CIPOS', 'CILEN', 'RE', 'AF', 'STD_SPAN', 'STD_POS', 'STD_POS1', 'STD_POS2', 'CUTPASTE'}

# DEL_NEGATIVE  True = negativo per le delezioni 
#               False = tutte le lunghezze positive 
DEL_NEGATIVE = True

# lunghezza minima delle SV da tenere (valore assoluto di SVLEN), 
# 0 = nessun filtro
MIN_SVLEN = 50


# True = per le inserzioni END diventa POS + SVLEN, così la regione ha una lunghezza
#        e si possono calcolare GC e coverage (stessa finestra che usa speedxx.py)
# False = END delle inserzioni lasciato com'è (standard VCF, END = POS)
INS_END_FROM_SVLEN = True

# True = tiene i BND (non validati nel benchmark), False = li elimina
KEEP_BND = False
# ==========================================================


def collect_filters(input_vcf):
    # legge il file una prima volta per sapere quali FILTER sono dichiarati nell'header e quali sono usati
    declared, used = set(), set()
    with open(input_vcf, 'r') as f:
        for line in f:
            if line.startswith('##FILTER=<ID='):
                declared.add(line.split('ID=', 1)[1].split(',', 1)[0].rstrip('>\n'))
            elif not line.startswith('#'):
                cols = line.split('\t')
                if len(cols) > 6:
                    for flt in cols[6].split(';'):
                        if flt not in ('PASS', '.'):
                            used.add(flt)
    return declared, used


def clean_vcf_info(input_vcf, output_vcf):

    missing_filters = set()
    if FILTER_MODE == "keep":
        declared, used = collect_filters(input_vcf)
        missing_filters = used - declared
    support_declared = False

    with open(input_vcf, 'r') as infile, open(output_vcf, 'w') as outfile:
        for line in infile:
            if line.startswith('#'): #header
                if line.startswith('##INFO=<ID=SUPPORT,'):
                    support_declared = True
                if line.startswith('#CHROM'): #prima della riga delle colonne aggiunge le dichiarazioni mancanti
                    for flt in sorted(missing_filters):
                        outfile.write(f'##FILTER=<ID={flt},Description="Filtro originale del caller">\n')
                    if not support_declared:
                        outfile.write('##INFO=<ID=SUPPORT,Number=1,Type=Integer,Description="Number of reads supporting the SV (RE for cuteSV)">\n')
                outfile.write(line)
                continue
            
            columns = line.strip().split('\t') #vcf line split
            if len(columns) < 10: #servono anche FORMAT e SAMPLE (colonne 9 e 10)
                outfile.write(line)
                continue
                
            pos = int(columns[1])
            info_field = columns[7]
            format_field = columns[8]
            sample_field = columns[9]
            chrom = columns[0]

            
            info_dict = {}
            for item in info_field.split(';'): #ciclo sulla colonna INFO
                if '=' in item:
                    key, val = item.split('=', 1) #INFO diventa dizionario
                    info_dict[key] = val
                else:
                    info_dict[item] = True # Tag booleani, ex: PRECISE
            
            
            # pulizia di SVTYPE (spostata prima perché serve per SVLEN ed END)
            if 'SVTYPE' in info_dict:
                if ',' in info_dict['SVTYPE']:
                    info_dict['SVTYPE'] = info_dict['SVTYPE'].split(',')[0] #prende sempre il primo valore nel caso sia una lista
            svtype = info_dict.get('SVTYPE', '')

            #Pulizia di SVLEN 
            if 'SVLEN' in info_dict:
                val = info_dict['SVLEN']
                if ',' in val:   #potrebbe esserci una lista "-500,505" data dalle varie lunghezze chiamate dai tool
                    first_val = val.split(',')[0] #prende il primo valore come riferimento
                    clean_len = abs(int(float(first_val)))
                else:
                    clean_len = abs(int(float(val)))
                if DEL_NEGATIVE and svtype == 'DEL': #delezioni negative come da convenzione VCF 4.2
                    clean_len = -clean_len
                info_dict['SVLEN'] = str(clean_len)
            
            #ricalcolo END 
            if 'END' in info_dict:
                if ',' in info_dict['END']: #anche END potrebbbe essere una lista
                    info_dict['END'] = info_dict['END'].split(',')[0] #prende sempre il primo valore
            if svtype != 'BND' and 'SVLEN' in info_dict: #nei BND END è la posizione sull'altro cromosoma, non va toccato
                if (INS_END_FROM_SVLEN and svtype == 'INS') or 'END' not in info_dict or int(info_dict['END']) < pos: #inserzione, manca END oppure end non può essere prima di start
                    info_dict['END'] = str(pos + abs(int(info_dict['SVLEN'])))
            
            # supporto in read uniformato: Sniffles e SVIM usano SUPPORT, cuteSV usa RE
            if 'RE' in info_dict and 'SUPPORT' not in info_dict:
                info_dict['SUPPORT'] = info_dict['RE']
            
            # Ricrea la stringa INFO 
            new_info_parts = []
            for k, v in info_dict.items():
                if k in INFO_TO_REMOVE: #info specifiche dei tool da eliminare (lista in alto)
                    continue
                else:
                    if v is True: #se è un booleano mette solo la key
                        new_info_parts.append(k)
                    else:
                        new_info_parts.append(f"{k}={v}") #inserisce key = value
            
            columns[7] = ';'.join(new_info_parts) #riscrive la colonna info separata da ;
            
            new_format = format_field.split(':')[0] #prende solo il genotipo (primo blocco GT) per evitare confusione
            columns[8] = new_format #riscrive la colonna format

            new_sample = sample_field.split(':')[0] #prende solo il genotipo (primo blocco GT) per evitare confusione
            columns[9] = new_sample #riscrive la colonna sample

            if FILTER_MODE == "pass":
                columns[6] = 'PASS' #sempre pass, niente codici strani di SVIM

            if svtype == 'BND' and not KEEP_BND:
                continue #BND esclusi

            #controllo lunghezza (solo se la variante ha SVLEN e non è un BND)
            if svtype != 'BND' and 'SVLEN' in info_dict:
                if abs(int(info_dict['SVLEN'])) < MIN_SVLEN:
                    continue #sv troppo piccola, non viene scritta

            #controllo CHROM
            if chrom != 'chrM':
                outfile.write('\t'.join(columns) + '\n') # riga pulita nel nuovo file


# se riceve input e output da riga di comando usa quelli, altrimenti i nomi di sempre
if len(sys.argv) == 3:
    clean_vcf_info(sys.argv[1], sys.argv[2])
else:
    clean_vcf_info("jasmine_union_3tools.vcf", "jasmine_union_3tools_clean.vcf")
