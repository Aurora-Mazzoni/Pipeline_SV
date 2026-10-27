import os
import sys
import argparse
import gzip
import pandas as pd
import pysam
pysam.set_verbosity(0)

import pyranges as pr

from sklearn.preprocessing import OneHotEncoder
from sklearn.compose import ColumnTransformer

import warnings
warnings.filterwarnings("ignore", category=FutureWarning)


def log_error(msg):
    print(f"Error: {msg}", file=sys.stderr)
    sys.exit(1)


def load_repeat_masker(file_path):
    df = pd.read_csv(
        file_path,
        sep="\t",
        comment="#",
        header=None,
        usecols=[0, 1, 2],  # BED del T2T: cromosoma, inizio, fine nelle prime 3 colonne
        names=["Chromosome", "Start", "End"]
    )
    df["Start"] = df["Start"].astype(int)
    df["End"] = df["End"].astype(int)
    return pr.PyRanges(df).merge()


def load_segmental_duplications(file_path):
    df = pd.read_csv(
        file_path,
        sep="\t",
        comment="#",
        header=None,
        usecols=[0, 1, 2],  # legge solo le prime 3 colonne anche se il file ne ha di più
        names=["Chromosome", "Start", "End"]
    )
    df["Start"] = df["Start"].astype(int)
    df["End"] = df["End"].astype(int)
    return pr.PyRanges(df).merge()


def load_gc_content(file_path):
    df = pd.read_csv(file_path, sep="\t", comment="#", header=None)
    df = df.iloc[:, [0, 1, 2, 6]].copy()
    df.columns = ["Chromosome", "Start", "End", "GC_CONTENT"]
    df["Start"] = df["Start"].astype(int)
    df["End"] = df["End"].astype(int)
    df["GC_CONTENT"] = pd.to_numeric(df["GC_CONTENT"], errors="coerce")
    return df


def load_coverage(file_path):
    with gzip.open(file_path, "rt") as f:
        df = pd.read_csv(
            f,
            sep="\t",
            header=None,
            usecols=[0, 1, 2, 4],
            names=["Chromosome", "Start", "End", "COVERAGE_MOSDEPTH"]
        )
    df["Start"] = df["Start"].astype(int)
    df["End"] = df["End"].astype(int)
    df["COVERAGE_MOSDEPTH"] = pd.to_numeric(df["COVERAGE_MOSDEPTH"], errors="coerce")
    return df


def load_feature_matrix(file_path):
    df = pd.read_csv(file_path, sep="\t")
    df = df.rename(columns={"CHROM": "Chromosome", "BPl": "Start", "BPr": "End"})

    required = ["Chromosome", "Start", "End"]
    missing = [c for c in required if c not in df.columns]
    if missing:
        log_error(
            "Missing required columns in feature matrix:\n" +
            "\n".join(f"- {c}" for c in missing)
        )

    df["Start"] = df["Start"].astype(int)
    df["End"] = df["End"].astype(int)
    return df


def process_vcf(vcf_path):
    rows = []
    with pysam.VariantFile(vcf_path) as vcf:
        sample = list(vcf.header.samples)[0] if len(vcf.header.samples) else "NA"

        for i, rec in enumerate(vcf):
            chrom = rec.chrom
            start = int(rec.pos)
            end = rec.stop or rec.info.get("END") or start
            end = int(end)

            svtype = rec.info.get("SVTYPE", "NA")
            svlen = rec.info.get("SVLEN", "NA")

            # chiave unica tecnica da portare fino al secondo script
            variant_key = f"{chrom}_{start}_{end}_{svtype}_{i}"

            supp = rec.info.get("SUPP_VEC", "000")
            supp = supp if isinstance(supp, str) and len(supp) == 3 else "000"
            sniffles, cutesv, svim = supp[0], supp[1], supp[2]  # ordine SUPP_VEC: Sniffles, cuteSV, SVIM

            # QUAL della variante (colonna 6 del VCF) al posto del GQ, che non c'è più dopo la pulizia
            qual = float(rec.qual) if rec.qual is not None else "."

            # tool da cui viene la riga scelta da Jasmine (serve per normalizzare il QUAL per tool):
            # legge dall'ID (ex "Sniffles2.DEL.12")
            rec_id = (rec.id or "").lower()
            if "sniffles" in rec_id:
                qual_source = "SNIFFLES"
            elif "cutesv" in rec_id:
                qual_source = "CUTESV"
            elif "svim" in rec_id:
                qual_source = "SVIM"
            elif supp.count("1") == 1:
                qual_source = ["SNIFFLES", "CUTESV", "SVIM"][supp.index("1")]
            else:
                qual_source = "UNKNOWN"

            # FILTER del caller: 1 = PASS, 0 = segnata come problematica, vuoto = nessun filtro 
            filtri = list(rec.filter.keys())
            filter_pass = (1 if "PASS" in filtri else 0) if filtri else "."

            # numero di read che supportano la variante (RE di cuteSV è già stato copiato in SUPPORT dalla pulizia)
            support = rec.info.get("SUPPORT")
            if isinstance(support, (list, tuple)):
                support = support[0] if support else None
            support = int(support) if support is not None else "."

            # genotipo come numero di copie alternative: 0/0 -> 0, 0/1 -> 1, 1/1 -> 2, mancante -> vuoto
            gt_tuple = rec.samples[sample].get("GT") or ()
            alleli = [a for a in gt_tuple if a is not None]
            gt = sum(1 for a in alleli if a > 0) if alleli else "."

            rows.append({
                "VARIANT_KEY": variant_key,
                "SAMPLE_NAME": sample,
                "Chromosome": chrom,
                "Start": start,
                "End": end,
                "SVTYPE_CALLER": svtype,
                "SVLEN_CALLER": svlen,
                "SNIFFLES_feature": int(sniffles) if str(sniffles).isdigit() else 0,
                "CUTESV_feature": int(cutesv) if str(cutesv).isdigit() else 0,
                "SVIM_feature": int(svim) if str(svim).isdigit() else 0,
                "QUAL": qual,
                "QUAL_SOURCE": qual_source,
                "FILTER_PASS": filter_pass,
                "SUPPORT": support,
                "GT": gt
            })

    df = pd.DataFrame(rows)

    # QUAL normalizzato: percentile del QUAL tra le chiamate dello stesso tool (da 0 a 100), per renderli confrontabili
    qual_num = pd.to_numeric(df["QUAL"], errors="coerce")
    df["QUAL_NORM"] = qual_num.groupby(df["QUAL_SOURCE"]).rank(pct=True) * 100
    print("[DEBUG] provenienza del QUAL:", df["QUAL_SOURCE"].value_counts().to_dict())

    return df


def add_overlap_flag(sv_pr, anno_pr, out_col):
    joined = sv_pr.join(anno_pr).as_df()
    sv_df = sv_pr.as_df()

    # se nessuna SV si sovrappone all'annotazione, la tabella del join è vuota e senza colonne: in quel caso il flag è False per tutte 
    if joined.empty:
        sv_df[out_col] = False
        return pr.PyRanges(sv_df)

    hits = joined[["Chromosome", "Start", "End"]].drop_duplicates()

    sv_df[out_col] = sv_df.merge(
        hits.assign(_hit=True),
        on=["Chromosome", "Start", "End"],
        how="left"
    )["_hit"].fillna(False)

    return pr.PyRanges(sv_df)


def _mean_on_overlaps_chunk(sv_chunk_df, anno_df, value_col, max_mean_span):
    sv_chunk_df = sv_chunk_df.copy()
    sv_chunk_df["_span"] = (sv_chunk_df["End"] - sv_chunk_df["Start"]).abs()

    short = sv_chunk_df[sv_chunk_df["_span"] <= max_mean_span]
    long_ = sv_chunk_df[sv_chunk_df["_span"] > max_mean_span]

    out = sv_chunk_df[["Chromosome", "Start", "End"]].copy()
    out[value_col] = pd.NA

    if not short.empty:
        chroms = short["Chromosome"].unique().tolist()
        min_s = int(short["Start"].min())
        max_e = int(short["End"].max())

        a = anno_df[
            (anno_df["Chromosome"].isin(chroms)) &
            (anno_df["End"] >= min_s) &
            (anno_df["Start"] <= max_e)
        ].copy()

        if not a.empty:
            sv_pr = pr.PyRanges(short[["Chromosome", "Start", "End"]])
            a_pr = pr.PyRanges(a[["Chromosome", "Start", "End", value_col]])

            j = sv_pr.join(a_pr).as_df()
            if not j.empty:
                means = (
                    j.groupby(["Chromosome", "Start", "End"], as_index=False)[value_col]
                    .mean()
                )
                out = out.merge(
                    means,
                    on=["Chromosome", "Start", "End"],
                    how="left",
                    suffixes=("", "_m")
                )
                out[value_col] = out[value_col].fillna(out[f"{value_col}_m"])
                out.drop(columns=[f"{value_col}_m"], inplace=True)

    if not long_.empty:
        for _, r in long_.iterrows():
            chrom = r["Chromosome"]
            pos = int(r["Start"])

            hit = anno_df[
                (anno_df["Chromosome"] == chrom) &
                (anno_df["Start"] <= pos) &
                (anno_df["End"] > pos)
            ]

            if not hit.empty:
                val = hit.iloc[0][value_col]
                out.loc[
                    (out["Chromosome"] == chrom) &
                    (out["Start"] == r["Start"]) &
                    (out["End"] == r["End"]),
                    value_col
                ] = val

    return out


def add_numeric_annotation_ram_safe(
    sv_df,
    anno_df,
    value_col,
    out_col,
    chunk_size=20000,
    max_mean_span=5_000_000
):
    results = []

    for chrom, sub in sv_df.groupby("Chromosome", sort=False):
        sub = sub.sort_values("Start").reset_index(drop=True)

        anno_chr = anno_df[anno_df["Chromosome"] == chrom].copy()
        if anno_chr.empty:
            tmp = sub[["Chromosome", "Start", "End"]].copy()
            tmp[out_col] = pd.NA
            results.append(tmp)
            continue

        for i in range(0, len(sub), chunk_size):
            chunk = sub.iloc[i:i + chunk_size]
            tmp = _mean_on_overlaps_chunk(chunk, anno_chr, value_col, max_mean_span)
            tmp = tmp.rename(columns={value_col: out_col})
            results.append(tmp)

    ann = pd.concat(results, ignore_index=True)
    merged = sv_df.merge(ann, on=["Chromosome", "Start", "End"], how="left")
    return merged


def main():
    parser = argparse.ArgumentParser()

    parser.add_argument("--repeat_masker", required=True)
    parser.add_argument("--segmental_duplications", required=True)
    parser.add_argument("--vcf_file", required=True)
    parser.add_argument("--gc_content", required=True)
    parser.add_argument("--coverage_mosdepth", required=True)
    parser.add_argument("--sv_feature_matrix", required=True)
    parser.add_argument("--output_path", required=True)

    parser.add_argument("--chunk_size", type=int, default=20000)
    parser.add_argument(
        "--max_mean_span",
        type=int,
        default=5_000_000,
        help="SV più lunghi di questo usano lookup su Start"
    )

    args = parser.parse_args()

    sv_df = process_vcf(args.vcf_file)
    if sv_df.empty:
        log_error("VCF produced 0 variants.")

    sv_df_raw = sv_df.copy()

    if not sv_df_raw["VARIANT_KEY"].is_unique:
        log_error("VARIANT_KEY is not unique in VCF.")

    print(f"[DEBUG] VCF rows: {len(sv_df_raw)}")
    print(f"[DEBUG] Unique VARIANT_KEY: {sv_df_raw['VARIANT_KEY'].nunique()}")

    sv_pr = pr.PyRanges(sv_df[["Chromosome", "Start", "End"]])
    repeats_pr = load_repeat_masker(args.repeat_masker)
    segdup_pr = load_segmental_duplications(args.segmental_duplications)

    sv_pr2 = add_overlap_flag(sv_pr, repeats_pr, "OVERLAPS_REPEATS")
    sv_pr3 = add_overlap_flag(sv_pr2, segdup_pr, "OVERLAPS_SEG_DUP")

    sv_df = sv_pr3.as_df().merge(
        sv_df,
        on=["Chromosome", "Start", "End"],
        how="left"
    )

    gc_df = load_gc_content(args.gc_content)
    cov_df = load_coverage(args.coverage_mosdepth)

    sv_df = add_numeric_annotation_ram_safe(
        sv_df,
        gc_df,
        value_col="GC_CONTENT",
        out_col="CG_CONTENT",
        chunk_size=args.chunk_size,
        max_mean_span=args.max_mean_span
    )

    sv_df = add_numeric_annotation_ram_safe(
        sv_df,
        cov_df,
        value_col="COVERAGE_MOSDEPTH",
        out_col="COVERAGE_MOSDEPTH",
        chunk_size=args.chunk_size,
        max_mean_span=args.max_mean_span
    )

    encoder = ColumnTransformer(
        [("ohe", OneHotEncoder(sparse_output=False), ["SVTYPE_CALLER"])],
        remainder="passthrough",
        verbose_feature_names_out=False
    ).set_output(transform="pandas")

    encoded_df = encoder.fit_transform(sv_df)

    feat_df = load_feature_matrix(args.sv_feature_matrix)

    if len(feat_df) != len(sv_df_raw):
        log_error(
            "Feature matrix and VCF do not have the same number of rows."
        )

    # assegna la chiave tecnica alla feature matrix mantenendo l'ordine
    feat_df["VARIANT_KEY"] = sv_df_raw["VARIANT_KEY"].values

    compare_df = sv_df_raw[["VARIANT_KEY", "Chromosome", "Start", "End"]].merge(
        feat_df[["VARIANT_KEY", "Chromosome", "Start", "End"]],
        on="VARIANT_KEY",
        how="inner",
        suffixes=("_vcf", "_feat")
    )

    same_order = (
        (compare_df["Chromosome_vcf"] == compare_df["Chromosome_feat"]) &
        (compare_df["Start_vcf"] == compare_df["Start_feat"]) &
        (compare_df["End_vcf"] == compare_df["End_feat"])
    )

    if (~same_order).sum() > 0:
        log_error("VCF and feature matrix row-order mismatch detected.")

    final_df = encoded_df.merge(
        feat_df,
        on="VARIANT_KEY",
        how="left",
        suffixes=("_encoded", "_feat")
    )

    rename_map = {}
    if "Chromosome_encoded" in final_df.columns:
        rename_map["Chromosome_encoded"] = "CHROM_CALLER"
    if "Start_encoded" in final_df.columns:
        rename_map["Start_encoded"] = "POS_CALLER"
    if "End_encoded" in final_df.columns:
        rename_map["End_encoded"] = "END_CALLER"

    final_df = final_df.rename(columns=rename_map)
    final_df = final_df.drop_duplicates(subset=["VARIANT_KEY"], keep="first")

    base = os.path.splitext(args.output_path)[0]
    final_df.to_csv(f"{base}.tsv", sep="\t", index=False)
    final_df.to_csv(f"{base}.csv", sep=",", index=False)

    print(f"[INFO] Written:")
    print(f"  {base}.tsv")
    print(f"  {base}.csv")


if __name__ == "__main__":
    main()
