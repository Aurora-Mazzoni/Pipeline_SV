import pysam, statistics, csv, math
pysam.set_verbosity(0)

from multiprocessing import Pool, cpu_count


import argparse
import warnings
warnings.filterwarnings("ignore", category=FutureWarning)

parser = argparse.ArgumentParser(description="Input for extrat features")

parser.add_argument("--cram_file", type=str, required=True, help="absolute path of cram file")
parser.add_argument("--vcf_file", type=str, required=True, help="absolute path of vcf ")
parser.add_argument("--out_file", type=str, required=True, help="output file")
parser.add_argument("--threads", type=int, default=None, help="numero di processi in parallelo (default: fino a 12)")
parser.add_argument("--checkpoint_dir", type=str, default=None, help="cartella dove salvare i blocchi già fatti, per riprendere dopo un'interruzione")

args = parser.parse_args()
cram_file = args.cram_file
vcf_file  = args.vcf_file
out_file  = args.out_file


# cram_file = "/shared/archive/ngsbo/migrated-from-ngsra/cram/genome_iit/HG002.100_reads.markdup.recal.cram"
# vcf_file  = "/shared/work/PI-tommaso.pippucci/RF-WGS/benchmark_sv/files_scr_paper/HG002/vcf/merged_jasmine2.vcf"
##vcf_file = "/shared/work/PI-tommaso.pippucci/RF-WGS/benchmark_sv/files_scr_paper/RealData/HG002_SVs_Tier1_v0.6.lft.hg38.INFO.highconf.fn.vcf.gz"
# out_file  = "sv_feature_matrix_speed.tsv"
#out_file = "snfn_features_matrix_speed.tsv"
# === PARAMETRI ===
NPROC     = args.threads if args.threads else min(12, cpu_count())
FLANK     = 50
MAX_TLEN  = 10000
MAX_READS_PER_REGION = 50000

def get_info_scalar(rec, key, default=None):
    v = rec.info.get(key, default)
    if isinstance(v, (list, tuple)) and v:
        v = v[0]
    return v

def normalize_svtype(rec):
    svt = get_info_scalar(rec, "SVTYPE", None)
    if isinstance(svt, bytes):
        svt = svt.decode()
    if svt is None and rec.alts and rec.alts[0].startswith("<"):
        svt = rec.alts[0].strip("<>")
    return svt or "NA"

def normalize_svlen(rec, svt):
    svlen = get_info_scalar(rec, "SVLEN", None)
    if svlen is None:
        svlen = (rec.stop - rec.pos) if svt != "INS" else 0
    try:
        return int(svlen)
    except Exception:
        return 0

class OnlineStats:
    def __init__(self):
        self.n = 0
        self.mean = 0.0
        self.M2 = 0.0
    def add(self, x):
        self.n += 1
        d = x - self.mean
        self.mean += d / self.n
        self.M2 += d * (x - self.mean)
    def result(self):
        if self.n == 0:
            return (None, None, 0)
        if self.n == 1:
            return (self.mean, 0.0, 1)
        return (self.mean, math.sqrt(self.M2 / (self.n - 1)), self.n)

# il BAM viene aperto UNA volta per ogni processo in parallelo (invece che a ogni regione)
_cram = None
def _init_worker():
    global _cram
    _cram = pysam.AlignmentFile(cram_file, "r")  # "r" riconosce da solo BAM o CRAM

# FEATURES ON REGIONS ===
def collect_features_region(args):
    chrom, s0, e0 = args
    cram = _cram
    n_reads = n_clipped = n_split = n_discordant = 0
    mapq_stats = OnlineStats()
    tlen_stats = OnlineStats()
    count = 0

    for r in cram.fetch(chrom, s0, e0):
        if r.is_secondary:
            continue
        n_reads += 1
        mapq_stats.add(r.mapping_quality or 0)

        # Clipping: conta la read se ha almeno un soft (S) o hard (H) clip
        # get_cigar_stats() conta le operazioni del CIGAR in C, molto più veloce che scorrere tutta la lista in Python 
        blocks = r.get_cigar_stats()[1]
        if blocks[4] > 0 or blocks[5] > 0:  # S/H
            n_clipped += 1

        if r.is_supplementary:
            n_split += 1

        if r.is_paired:
            if not r.is_proper_pair:
                n_discordant += 1
            t = abs(r.template_length)
            if 0 < t < MAX_TLEN:
                tlen_stats.add(t)

        count += 1
        if count >= MAX_READS_PER_REGION:
            break

    mean_insert, sd_insert, _ = tlen_stats.result()
    mean_mapq, _, _ = mapq_stats.result()
    return {
        "coverage": n_reads,
        "mean_insert": mean_insert,
        "sd_insert": sd_insert,
        "mean_mapq": mean_mapq,
        "n_clipped": n_clipped,
        "n_split": n_split,
        "n_discordant": n_discordant,
    }

def build_windows(rec):
    chrom = rec.chrom
    BPl = int(rec.pos)
    BPr = int(rec.stop)
    svt = normalize_svtype(rec)
    svlen = normalize_svlen(rec, svt)

    left_start0  = max(0, (BPl - FLANK) - 1)
    left_end0    = max(left_start0, BPl - 1)

    if svt == "INS":
        L = max(1, abs(svlen))
        inside_start0 = BPl - 1
        inside_end0   = inside_start0 + L
        right_start0  = inside_end0
        right_end0    = right_start0 + FLANK
    else:
        inside_start0 = BPl - 1
        inside_end0   = BPr
        right_start0  = BPr
        right_end0    = BPr + FLANK

    # --- nuovi campi dal VCF ---
    filt = ";".join(list(rec.filter.keys())) if rec.filter else "PASS"
    supp_vec = str(get_info_scalar(rec, "SUPP_VEC", "NA"))
    precision = "IMPRECISE" if "IMPRECISE" in rec.info else ("PRECISE" if "PRECISE" in rec.info else "NA")

    # --- parsing SUPP_VEC per i tre caller ---
    # ordine SUPP_VEC: Sniffles, cuteSV, SVIM (manta = Sniffles, delly = cuteSV, smoove = SVIM)
    manta, delly, smoove = (0, 0, 0)
    if supp_vec.isdigit() and len(supp_vec) >= 3:
        manta  = int(supp_vec[0])
        delly  = int(supp_vec[1])
        smoove = int(supp_vec[2])

    return (chrom, BPl, BPr, svt, svlen, filt, supp_vec, manta, delly, smoove, precision,
            (chrom, left_start0, left_end0),
            (chrom, inside_start0, inside_end0),
            (chrom, right_start0, right_end0))

# === PREPARA JOBS ===
vcf = pysam.VariantFile(vcf_file, "r")
jobs, meta = [], []
for rec in vcf:
    w = build_windows(rec)
    if not w:
        continue
    chrom, BPl, BPr, svt, svlen, filt, supp_vec, manta, delly, smoove, precision, L, I, R = w
    start_idx = len(jobs)
    jobs.extend([L, I, R])
    meta.append((chrom, BPl, BPr, svt, svlen, filt, supp_vec, manta, delly, smoove, precision, start_idx))
vcf.close()

# === ESEGUI PARALLEL ===
# Le regioni vengono analizzate a blocchi di 5000
# Se è indicata --checkpoint_dir, ogni blocco finito viene salvato su disco così se l'esecuzione si interrompe i blocchi già fatti vengono riecuperati
import os, pickle, hashlib
BLOCK = 5000
print(f"[speedxx] {len(meta)} varianti, {len(jobs)} regioni da analizzare", flush=True)

ckpt = None
if args.checkpoint_dir:
    # firma serve così se il VCF cambia cambia la sottocartella e i vecchi blocchi non vengono riciclati
    firma = hashlib.md5(repr(jobs).encode()).hexdigest()[:12]
    ckpt = os.path.join(args.checkpoint_dir, f"regioni_{firma}")
    os.makedirs(ckpt, exist_ok=True)

feats = []
with Pool(processes=NPROC, initializer=_init_worker) as pool:
    for b in range(0, len(jobs), BLOCK):
        fatte = min(b + BLOCK, len(jobs))
        fn = os.path.join(ckpt, f"blocco_{b // BLOCK:05d}.pkl") if ckpt else None
        if fn and os.path.exists(fn):
            with open(fn, "rb") as fh:
                part = pickle.load(fh)
            print(f"[speedxx] {fatte}/{len(jobs)} regioni (recuperate dal checkpoint)", flush=True)
        else:
            part = pool.map(collect_features_region, jobs[b:b + BLOCK], chunksize=50)
            if fn:
                tmp = fn + ".tmp"  # scrive prima un file temporaneo, così un'interruzione non lascia blocchi a metà
                with open(tmp, "wb") as fh:
                    pickle.dump(part, fh)
                os.replace(tmp, fn)
            print(f"[speedxx] {fatte}/{len(jobs)} regioni fatte", flush=True)
        feats.extend(part)

# === SCRIVI TSV ===
with open(out_file, "w", newline="") as tsv:
    wtr = csv.writer(tsv, delimiter="\t")
    wtr.writerow([
        "CHROM","BPl","BPr","SVTYPE","SVLEN","FILTER","SUPP_VEC","SNIFFLES","CUTESV","SVIM","PRECISION",
        "coverage_inside","mean_insert_inside","sd_insert_inside","mean_mapq_inside",
        "n_clipped_inside","n_split_inside","n_discordant_inside",
        "coverage_left","coverage_right","mean_flank_insert","Δinsert",
        "mean_flank_mapq","Δmapq",
        "clipped_ratio","split_ratio","discordant_ratio"
    ])

    for chrom, BPl, BPr, svt, svlen, filt, supp_vec, manta, delly, smoove, precision, idx in meta:
        left, inside, right = feats[idx], feats[idx+1], feats[idx+2]

        # flanks mean
        flank_ins = [x for x in [left["mean_insert"], right["mean_insert"]] if x is not None]
        flank_map = [x for x in [left["mean_mapq"], right["mean_mapq"]] if x is not None]
        mean_flank_insert = statistics.mean(flank_ins) if flank_ins else None
        mean_flank_mapq   = statistics.mean(flank_map) if flank_map else None

        delta_insert = (inside["mean_insert"] - mean_flank_insert) if (inside["mean_insert"] is not None and mean_flank_insert is not None) else None
        delta_mapq   = (inside["mean_mapq"] - mean_flank_mapq) if (inside["mean_mapq"] is not None and mean_flank_mapq is not None) else None

        cov_in = inside["coverage"] or 0
        clipped_ratio    = (inside["n_clipped"]/cov_in)    if cov_in>0 else None
        split_ratio      = (inside["n_split"]/cov_in)      if cov_in>0 else None
        discordant_ratio = (inside["n_discordant"]/cov_in) if cov_in>0 else None

        wtr.writerow([
            chrom, BPl, BPr, svt, svlen, filt, supp_vec, manta, delly, smoove, precision,
            inside["coverage"], inside["mean_insert"], inside["sd_insert"], inside["mean_mapq"],
            inside["n_clipped"], inside["n_split"], inside["n_discordant"],
            left["coverage"], right["coverage"],
            mean_flank_insert, delta_insert,
            mean_flank_mapq, delta_mapq,
            clipped_ratio, split_ratio, discordant_ratio
        ])

print(f"\n✅ TSV done: {out_file}")

