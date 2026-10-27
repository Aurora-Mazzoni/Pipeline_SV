import argparse
import pandas as pd

SVTYPES = ["DEL", "DUP", "INS", "INV"]

# colonne finali per VarTrustML (long read): tolte le 6 colonne su insert size, coppie discordanti clipped_ratio/n_clipped_inside 
#  aggiunte le sovrapposizioni con repeatmasker e segdup, QUAL normalizzato, FILTER, GT e frazione di read a supporto (divisa per coverage)
FINAL_COLS = [
    "SNIFFLES","CUTESV","SVIM","QUAL_NORM","FILTER_PASS","SUPPORT_RATIO","GT","CG_CONTENT","COVERAGE_MOSDEPTH",
    "OVERLAPS_REPEATS","OVERLAPS_SEG_DUP",
    "coverage_inside","mean_mapq_inside",
    "n_split_inside",
    "coverage_left","coverage_right",
    "mean_flank_mapq","DELTA_mapq",
    "split_ratio",
    "SVLEN_CALLER"
]


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Prepare final datasets for VarTrustML"
    )
    parser.add_argument("input_file", help="Input TSV file (tab-separated).")
    parser.add_argument("-o", "--out_prefix", default="dataset",
                        help="Output prefix (default: dataset).")
    args = parser.parse_args()

    df = pd.read_csv(args.input_file, sep="\t", low_memory=False)

    df = df.replace({".": pd.NA, "": pd.NA, "NA": pd.NA, "NaN": pd.NA, "nan": pd.NA})

    if "SVLEN_CALLER" in df.columns:
        s = (
            df["SVLEN_CALLER"]
            .astype(str)
            .str.replace(r"[()\s,]", "", regex=True)
            .replace({"": pd.NA, "None": pd.NA, "nan": pd.NA})
        )
        df["SVLEN_CALLER"] = pd.to_numeric(s, errors="coerce")

    keep = [
        "VARIANT_KEY",
        "SNIFFLES_feature", "CUTESV_feature", "SVIM_feature",
        "QUAL_NORM", "FILTER_PASS", "SUPPORT", "GT", "CG_CONTENT", "COVERAGE_MOSDEPTH",
        "OVERLAPS_REPEATS", "OVERLAPS_SEG_DUP",
        "coverage_inside", "mean_mapq_inside",
        "n_split_inside",
        "coverage_left", "coverage_right",
        "mean_flank_mapq", "Δmapq",
        "split_ratio",
        "SVLEN_CALLER",
        "SVTYPE",
    ]

    missing = [c for c in keep if c not in df.columns]
    if missing:
        raise SystemExit(
            "Colonne mancanti nel file:\n" + "\n".join(f"- {c}" for c in missing)
        )

    df = df[keep].rename(columns={
        "SNIFFLES_feature": "SNIFFLES",
        "CUTESV_feature": "CUTESV",
        "SVIM_feature": "SVIM",
        "Δmapq": "DELTA_mapq",
    })

    if df["VARIANT_KEY"].isna().any():
        raise SystemExit("VARIANT_KEY with missing values")

    if not df["VARIANT_KEY"].is_unique:
        raise SystemExit("VARIANT_KEY duplicated")

    # sovrapposizioni: da True/False a 1/0, come le altre feature numeriche
    for c in ["OVERLAPS_REPEATS", "OVERLAPS_SEG_DUP"]:
        df[c] = df[c].astype(str).map({"True": 1, "False": 0})

    # SUPPORT diviso la coverage media ai lati della SV (vuota se i fianchi non hanno read)
    fianchi = (pd.to_numeric(df["coverage_left"], errors="coerce") +
               pd.to_numeric(df["coverage_right"], errors="coerce")) / 2
    df["SUPPORT_RATIO"] = pd.to_numeric(df["SUPPORT"], errors="coerce") / fianchi.where(fianchi > 0)

    for c in FINAL_COLS:
        df[c] = pd.to_numeric(df[c], errors="coerce")

    missing_counts = df[FINAL_COLS].isna().sum()
    print("\nMissing:")
    print(missing_counts[missing_counts > 0])

    print("\nNumber of missing:")
    print(df[FINAL_COLS].isna().any(axis=1).sum())

    df["SVTYPE"] = df["SVTYPE"].astype(str).str.upper()
    # sottotipi (DUP:INT, DUP:TANDEM) -> tipo principale (DUP), così finiscono nel file delle DUP
    df["SVTYPE"] = df["SVTYPE"].str.split(":").str[0]

    # VARIANT_KEY diventa row_id
    df = df.set_index("VARIANT_KEY")
    df.index.name = "row_id"

    for sv in SVTYPES:
        out = f"{args.out_prefix}_{sv}.csv"
        sub = df[df["SVTYPE"] == sv][FINAL_COLS]
        sub.to_csv(out, index=True, index_label="row_id")

    print("Files:")
    for sv in SVTYPES:
        print(f"  {args.out_prefix}_{sv}.csv")


if __name__ == "__main__":
    main()
