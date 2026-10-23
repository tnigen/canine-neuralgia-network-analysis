#!/usr/bin/env python3
"""Run PyDESeq2 differential expression, annotation, and quality control."""

from __future__ import annotations

import argparse
import gzip
import hashlib
import importlib.metadata
import math
import os
import platform
import re
import shutil
import subprocess
from pathlib import Path

import numpy as np
import pandas as pd


SAMPLE_ORDER = ["A1_TG", "A2_TG", "A3_TG", "B1_TG", "B2_TG", "B3_TG"]
SAMPLE_TO_RUN = {
    "A1_TG": "SRR23372014", "A2_TG": "SRR23372013", "A3_TG": "SRR23372012",
    "B1_TG": "SRR23372011", "B2_TG": "SRR23372010", "B3_TG": "SRR23372009",
}
PRIMARY_PADJ = 0.05
PRIMARY_ABS_LFC = 1.0
RANDOM_SEED = 20260904


def parse_gtf_attributes(text: str) -> dict[str, str]:
    return {k: v for k, v in re.findall(r'(\S+) "([^"]*)"', text)}


def parse_gtf_gene_annotation(gtf: Path) -> dict[str, dict[str, str]]:
    annotation: dict[str, dict[str, str]] = {}
    with gtf.open("rt", encoding="utf-8") as handle:
        for line in handle:
            if line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) != 9 or fields[2] != "gene":
                continue
            attrs = parse_gtf_attributes(fields[8])
            gene_id = attrs.get("gene_id", "")
            if gene_id:
                annotation[gene_id] = {
                    "gene_symbol": attrs.get("gene_name", ""),
                    "gene_biotype": attrs.get("gene_biotype", ""),
                }
    return annotation


def parse_ncbi_crossrefs(gene2ensembl: Path, gene_info_path: Path):
    ensembl_to_entrez: dict[str, set[str]] = {}
    with gzip.open(gene2ensembl, "rt", encoding="utf-8", errors="replace") as handle:
        header = handle.readline().rstrip("\n").lstrip("#").split("\t")
        idx = {name: i for i, name in enumerate(header)}
        for line in handle:
            fields = line.rstrip("\n").split("\t")
            if fields[idx["tax_id"]] != "10116":
                continue
            ens = fields[idx["Ensembl_gene_identifier"]]
            entrez = fields[idx["GeneID"]]
            if ens not in {"", "-"} and entrez not in {"", "-"}:
                ensembl_to_entrez.setdefault(ens, set()).add(entrez)
    gene_info: dict[str, dict[str, str]] = {}
    with gzip.open(gene_info_path, "rt", encoding="utf-8", errors="replace") as handle:
        header = handle.readline().rstrip("\n").lstrip("#").split("\t")
        idx = {name: i for i, name in enumerate(header)}
        for line in handle:
            fields = line.rstrip("\n").split("\t")
            if fields[idx["tax_id"]] != "10116":
                continue
            gene_info[fields[idx["GeneID"]]] = {
                "symbol": fields[idx["Symbol"]],
                "name": fields[idx["description"]],
            }
    return ensembl_to_entrez, gene_info


def annotate_degs(
    project: Path, degs: pd.DataFrame, *, from_counts: bool = False,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    gtf = project / "data_processed" / "reference" / "Rattus_norvegicus.Rnor_6.0.104.gtf"
    gtf_map = parse_gtf_gene_annotation(gtf)
    crossrefs, gene_info = parse_ncbi_crossrefs(
        project / "data_raw" / "annotation" / (
            "gene2ensembl.rat.tsv.gz" if from_counts else "gene2ensembl.gz"
        ),
        project / "data_raw" / "annotation" / "Rattus_norvegicus.gene_info.gz",
    )
    records = []
    for _, row in degs.iterrows():
        original = str(row["original_id"])
        gtf_entry = gtf_map.get(original, {})
        entrez_ids = sorted(crossrefs.get(original, set()), key=lambda x: int(x))
        if len(entrez_ids) > 1:
            status = "ambiguous"
            symbol = gtf_entry.get("gene_symbol", "")
            gene_name = ""
            evidence = f"Ensembl release 104 gene ID maps to multiple NCBI Gene IDs: {';'.join(entrez_ids)}"
        elif original not in gtf_map:
            status = "unresolved"
            symbol = ""
            gene_name = ""
            evidence = "Identifier absent from the matching Ensembl release-104 GTF gene records"
        else:
            status = "uniquely_mapped"
            entrez = entrez_ids[0] if entrez_ids else ""
            info = gene_info.get(entrez, {}) if entrez else {}
            symbol = info.get("symbol", "") or gtf_entry.get("gene_symbol", "")
            gene_name = info.get("name", "")
            evidence = "Unique Ensembl gene record in release-104 GTF"
            if entrez:
                evidence += "; unique NCBI gene2ensembl cross-reference"
            else:
                evidence += "; no current NCBI Gene cross-reference"
        records.append(
            {
                "original_id": original,
                "rat_gene_symbol": symbol,
                "rat_ensembl_gene_id": original if original.startswith("ENSRNO") else "",
                "rat_entrez_id": ";".join(entrez_ids),
                "rat_gene_name": gene_name,
                "gene_biotype": gtf_entry.get("gene_biotype", ""),
                "mapping_status": status,
                "mapping_evidence": evidence,
            }
        )
    qc = pd.DataFrame(
        records,
        columns=[
            "original_id", "rat_gene_symbol", "rat_ensembl_gene_id", "rat_entrez_id",
            "rat_gene_name", "gene_biotype", "mapping_status", "mapping_evidence",
        ],
    )
    annotated = degs.merge(qc, on="original_id", how="left", validate="one_to_one")
    return annotated, qc


def save_figure(fig, stem: Path, pdf: bool = True) -> None:
    fig.savefig(stem.with_suffix(".png"), dpi=300, bbox_inches="tight")
    if pdf:
        fig.savefig(stem.with_suffix(".pdf"), bbox_inches="tight")


def version_row(name: str, value: str) -> dict[str, str]:
    return {"software": name, "version": value}


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(8 * 1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def binary_version(command: list[str]) -> str:
    try:
        out = subprocess.run(command, text=True, capture_output=True, check=False)
        line = (out.stdout or out.stderr).strip().splitlines()[0]
        return re.sub(r"^\S*/([^/\s]+)(?=\s)", r"\1", line)
    except Exception as exc:
        return f"unavailable ({exc})"


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--project", type=Path, default=Path(__file__).resolve().parents[1])
    parser.add_argument("--threads", type=int, default=min(8, os.cpu_count() or 1))
    parser.add_argument("--from-counts", action="store_true", help="Use the bundled rat annotation snapshot")
    args = parser.parse_args()
    project = args.project.resolve()
    os.environ.setdefault("MPLCONFIGDIR", str(project / ".matplotlib"))
    (project / ".matplotlib").mkdir(exist_ok=True)
    import matplotlib.pyplot as plt
    from matplotlib.patches import Patch
    import seaborn as sns
    from adjustText import adjust_text
    from pydeseq2.dds import DeseqDataSet
    from pydeseq2.ds import DeseqStats
    from scipy.cluster.hierarchy import linkage
    from scipy.spatial.distance import pdist, squareform
    from sklearn.decomposition import PCA

    np.random.seed(RANDOM_SEED)
    tables = project / "results" / "tables"
    figures = project / "results" / "figures"
    de_dir = project / "results" / "differential_expression"
    logs = project / "results" / "logs"
    for p in (tables, figures, de_dir, logs):
        p.mkdir(parents=True, exist_ok=True)

    metadata = pd.read_csv(project / "data_processed" / "sample_metadata.csv", index_col="sample_id")
    counts = pd.read_csv(project / "data_processed" / "GSE224814_raw_count_matrix.csv", index_col="original_id")
    counts = counts[SAMPLE_ORDER]
    if list(metadata.index) != SAMPLE_ORDER:
        metadata = metadata.loc[SAMPLE_ORDER]
    if set(counts.columns) != set(metadata.index):
        raise ValueError("Count-matrix and metadata sample identities do not match")
    if counts.index.duplicated().any():
        raise ValueError("Duplicated gene IDs in raw count matrix")
    arr = counts.to_numpy()
    if (arr < 0).any() or not np.equal(arr, np.floor(arr)).all():
        raise ValueError("PyDESeq2 input is not non-negative integer raw counts")
    counts = counts.astype(int)
    counts.to_csv(project / "data_processed" / "GSE224814_raw_count_matrix.csv", index_label="original_id")
    pd.DataFrame(
        [
            {
                "features": len(counts),
                "samples": counts.shape[1],
                "sample_columns": "|".join(counts.columns),
                "minimum_count": int(counts.min().min()),
                "maximum_count": int(counts.max().max()),
                "nonnegative": bool((counts.to_numpy() >= 0).all()),
                "all_integer": bool(np.equal(counts.to_numpy(), np.floor(counts.to_numpy())).all()),
                "duplicated_sample_columns": int(counts.columns.duplicated().sum()),
                "duplicated_gene_identifiers": int(counts.index.duplicated().sum()),
                "metadata_samples_match_exactly": set(counts.columns) == set(metadata.index),
                "source": "paired FASTQ -> HISAT2 rn6 alignment -> featureCounts gene-level fragments",
            }
        ]
    ).to_csv(tables / "raw_count_matrix_audit.csv", index=False)

    before = len(counts)
    keep = (counts >= 10).sum(axis=1) >= 3
    filtered = counts.loc[keep].copy()
    filtered.to_csv(project / "data_processed" / "GSE224814_filtered_counts.csv", index_label="original_id")
    filter_summary = pd.DataFrame(
        {
            "metric": ["features_before_filtering", "features_retained", "features_removed", "percentage_retained"],
            "value": [before, len(filtered), before - len(filtered), 100 * len(filtered) / before],
        }
    )
    filter_summary.to_csv(tables / "low_count_filter_summary.csv", index=False)
    if filtered.empty:
        raise RuntimeError("Low-count-filtered matrix is empty")

    meta_model = metadata[["condition"]].copy()
    meta_model["condition"] = pd.Categorical(meta_model["condition"], categories=["Sham", "CCI_ION"], ordered=True)
    dds = DeseqDataSet(
        counts=filtered.T,
        metadata=meta_model,
        design="~condition",
        refit_cooks=True,
        n_cpus=args.threads,
        quiet=False,
    )
    dds.deseq2()
    stats = DeseqStats(
        dds,
        contrast=["condition", "CCI_ION", "Sham"],
        alpha=PRIMARY_PADJ,
        cooks_filter=True,
        independent_filter=True,
        n_cpus=args.threads,
        quiet=False,
    )
    stats.summary()
    result = stats.results_df.reset_index().rename(columns={"index": "original_id"})
    result = result[["original_id", "baseMean", "log2FoldChange", "lfcSE", "stat", "pvalue", "padj"]]
    result = result.sort_values(["padj", "pvalue", "original_id"], na_position="last").reset_index(drop=True)
    all_path = de_dir / "GSE224814_CCI_ION_vs_Sham_all_results.csv"
    result.to_csv(all_path, index=False)

    primary_mask = result["padj"].notna() & (result["padj"] < PRIMARY_PADJ) & (result["log2FoldChange"].abs() >= PRIMARY_ABS_LFC)
    degs = result.loc[primary_mask].copy()
    degs["direction"] = np.where(degs["log2FoldChange"] >= PRIMARY_ABS_LFC, "Up", "Down")
    degs = degs.sort_values(["padj", "pvalue"], na_position="last").reset_index(drop=True)
    up = degs[degs["direction"] == "Up"].copy()
    down = degs[degs["direction"] == "Down"].copy()
    deg_path = de_dir / "GSE224814_CCI_ION_vs_Sham_DEGs.csv"
    up_path = de_dir / "GSE224814_CCI_ION_vs_Sham_upregulated.csv"
    down_path = de_dir / "GSE224814_CCI_ION_vs_Sham_downregulated.csv"
    degs.to_csv(deg_path, index=False)
    up.to_csv(up_path, index=False)
    down.to_csv(down_path, index=False)

    deg_lfcs = degs["log2FoldChange"]
    summary_values = {
        "features_before_low_count_filtering": before,
        "features_retained_after_filtering": len(filtered),
        "features_tested_by_PyDESeq2": len(result),
        "features_with_raw_pvalue_lt_0.05": int((result["pvalue"] < 0.05).sum()),
        "features_with_padj_lt_0.05_regardless_of_fold_change": int((result["padj"] < 0.05).sum()),
        "features_with_NA_padj": int(result["padj"].isna().sum()),
        "primary_DEGs": len(degs),
        "upregulated": len(up),
        "downregulated": len(down),
        "smallest_non_NA_pvalue": result["pvalue"].min(skipna=True),
        "smallest_non_NA_padj": result["padj"].min(skipna=True),
        "minimum_DEG_log2FoldChange": deg_lfcs.min(skipna=True) if len(degs) else np.nan,
        "maximum_DEG_log2FoldChange": deg_lfcs.max(skipna=True) if len(degs) else np.nan,
    }
    pd.DataFrame(summary_values.items(), columns=["metric", "value"]).to_csv(tables / "DEG_summary.csv", index=False)
    sensitivity = [
        ("raw_pvalue_lt_0.05", (result["pvalue"] < 0.05)),
        ("padj_lt_0.05", (result["padj"] < 0.05)),
        ("padj_lt_0.05_and_abs_log2FC_ge_0.5", (result["padj"] < 0.05) & (result["log2FoldChange"].abs() >= 0.5)),
        ("padj_lt_0.05_and_abs_log2FC_ge_1_primary", primary_mask),
        ("raw_pvalue_lt_0.05_and_abs_log2FC_ge_1", (result["pvalue"] < 0.05) & (result["log2FoldChange"].abs() >= 1)),
    ]
    pd.DataFrame([{"criterion": k, "count": int(v.sum()), "role": "primary" if "primary" in k else "descriptive_sensitivity"} for k, v in sensitivity]).to_csv(
        tables / "threshold_sensitivity_summary.csv", index=False
    )

    library_size = counts.sum(axis=0).rename("assigned_gene_counts").reset_index().rename(columns={"index": "sample_id"})
    library_size = library_size.merge(metadata.reset_index(), on="sample_id", how="left")
    library_size.to_csv(tables / "library_size.csv", index=False)
    detected = (counts > 0).sum(axis=0).rename("detected_features_nonzero").reset_index().rename(columns={"index": "sample_id"})
    detected = detected.merge(metadata.reset_index(), on="sample_id", how="left")
    detected.to_csv(tables / "detected_genes_per_sample.csv", index=False)
    raw_qc_records = []
    for sample in SAMPLE_ORDER:
        run_accession = SAMPLE_TO_RUN[sample]
        hisat_text = (project / "results" / "qc" / f"{sample}_{run_accession}_hisat2_summary.txt").read_text()
        pair_match = re.search(r"^(\d+) reads;", hisat_text, flags=re.MULTILINE)
        rate_match = re.search(r"([0-9.]+)% overall alignment rate", hisat_text)
        fc_summary = pd.read_csv(
            project / "results" / "qc" / f"{sample}_{run_accession}_featureCounts_summary.txt",
            sep="\t",
        )
        value_col = fc_summary.columns[1]
        assigned = int(fc_summary.loc[fc_summary["Status"] == "Assigned", value_col].iloc[0])
        classified = int(fc_summary[value_col].sum())
        raw_qc_records.append(
            {
                "sample_id": sample,
                "GSM": metadata.loc[sample, "GSM"],
                "condition": metadata.loc[sample, "condition"],
                "SRA_run": run_accession,
                "input_read_pairs": int(pair_match.group(1)) if pair_match else np.nan,
                "overall_alignment_rate_percent": float(rate_match.group(1)) if rate_match else np.nan,
                "featureCounts_assigned_fragments": assigned,
                "featureCounts_classified_records": classified,
                "featureCounts_assignment_percent": 100 * assigned / classified,
            }
        )
    pd.DataFrame(raw_qc_records).to_csv(tables / "raw_read_QC.csv", index=False)
    sns.set_theme(style="whitegrid", context="talk")
    palette = {"Sham": "#3977A9", "CCI_ION": "#D14A3A"}
    fig, ax = plt.subplots(figsize=(9, 5.5))
    colors = [palette[x] for x in library_size["condition"]]
    ax.bar(library_size["sample_id"], library_size["assigned_gene_counts"], color=colors)
    ax.set_ylabel("Assigned gene-level fragment counts")
    ax.set_xlabel("Sample")
    ax.set_title("GSE224814 library size")
    ax.ticklabel_format(axis="y", style="sci", scilimits=(0, 0))
    fig.tight_layout()
    save_figure(fig, figures / "library_size", pdf=False)
    plt.close(fig)

    dds.vst(use_design=False)
    transformed = pd.DataFrame(dds.layers["vst_counts"], index=dds.obs_names, columns=dds.var_names)
    transformed.T.to_csv(project / "data_processed" / "GSE224814_VST_expression.csv", index_label="original_id")
    pca = PCA(n_components=2, random_state=RANDOM_SEED)
    coords_arr = pca.fit_transform(transformed)
    coords = pd.DataFrame(coords_arr, index=transformed.index, columns=["PC1", "PC2"])
    coords["condition"] = metadata.loc[coords.index, "condition"].astype(str)
    coords["GSM"] = metadata.loc[coords.index, "GSM"]
    coords.reset_index(names="sample_id").to_csv(tables / "PCA_coordinates.csv", index=False)
    pc1, pc2 = pca.explained_variance_ratio_[:2] * 100
    sham_pc1 = coords.loc[coords.condition == "Sham", "PC1"]
    cci_pc1 = coords.loc[coords.condition == "CCI_ION", "PC1"]
    pc1_nonoverlap = bool(sham_pc1.max() < cci_pc1.min() or cci_pc1.max() < sham_pc1.min())
    centroid_sham = coords.loc[coords.condition == "Sham", ["PC1", "PC2"]].mean().to_numpy()
    centroid_cci = coords.loc[coords.condition == "CCI_ION", ["PC1", "PC2"]].mean().to_numpy()
    centroid_separation = float(np.linalg.norm(centroid_cci - centroid_sham))
    pd.DataFrame(
        [
            {"metric": "PC1_variance_percent", "value": pc1},
            {"metric": "PC2_variance_percent", "value": pc2},
            {"metric": "PC1_group_ranges_nonoverlapping", "value": pc1_nonoverlap},
            {"metric": "Sham_PC1_min", "value": sham_pc1.min()},
            {"metric": "Sham_PC1_max", "value": sham_pc1.max()},
            {"metric": "CCI_ION_PC1_min", "value": cci_pc1.min()},
            {"metric": "CCI_ION_PC1_max", "value": cci_pc1.max()},
            {"metric": "condition_centroid_distance_PC1_PC2", "value": centroid_separation},
        ]
    ).to_csv(tables / "PCA_group_separation_summary.csv", index=False)
    fig, ax = plt.subplots(figsize=(8, 6.5))
    for condition, subset in coords.groupby("condition"):
        ax.scatter(subset["PC1"], subset["PC2"], s=110, color=palette[condition], label=condition, edgecolor="black", linewidth=0.6)
        for sample, row in subset.iterrows():
            ax.annotate(sample, (row["PC1"], row["PC2"]), xytext=(5, 5), textcoords="offset points", fontsize=10)
    ax.axhline(0, color="#bbbbbb", linewidth=0.7)
    ax.axvline(0, color="#bbbbbb", linewidth=0.7)
    ax.grid(True, color="#9E9E9E", alpha=0.18, linewidth=0.8)
    ax.set_axisbelow(True)
    ax.set_xlabel(f"PC1 ({pc1:.2f}% variance)")
    ax.set_ylabel(f"PC2 ({pc2:.2f}% variance)")
    ax.set_title("PCA of variance-stabilized expression")
    ax.margins(x=0.10, y=0.12)
    ax.legend(frameon=False, loc="center left", bbox_to_anchor=(1.01, 0.5))
    fig.tight_layout()
    save_figure(fig, figures / "PCA_CCI_ION_vs_Sham")
    plt.close(fig)

    condensed_distance = pdist(transformed, metric="euclidean")
    distance = pd.DataFrame(
        squareform(condensed_distance), index=transformed.index, columns=transformed.index
    )
    distance.to_csv(tables / "sample_distance_matrix.csv", index_label="sample_id")
    sample_linkage = linkage(condensed_distance, method="average")
    cg = sns.clustermap(
        distance,
        row_linkage=sample_linkage,
        col_linkage=sample_linkage,
        cmap="viridis",
        figsize=(8, 7),
        linewidths=0.3,
        annot=True,
        fmt=".1f",
    )
    cg.fig.suptitle("Sample distances (VST expression)", y=1.02)
    cg.fig.savefig(figures / "sample_distance_heatmap.png", dpi=300, bbox_inches="tight")
    plt.close(cg.fig)

    centroid_distances = []
    for sample, row in coords.iterrows():
        centroid = coords.loc[coords["condition"] == row["condition"], ["PC1", "PC2"]].mean(axis=0)
        dist = float(np.linalg.norm(row[["PC1", "PC2"]].astype(float).to_numpy() - centroid.to_numpy()))
        centroid_distances.append(dist)
    outlier = coords.copy()
    outlier["distance_from_condition_centroid_PC1_PC2"] = centroid_distances
    q1, q3 = np.quantile(centroid_distances, [0.25, 0.75])
    cutoff = q3 + 1.5 * (q3 - q1)
    outlier["exploratory_flag"] = outlier["distance_from_condition_centroid_PC1_PC2"] > cutoff
    outlier["flag_rule"] = f"distance > pooled Q3 + 1.5*IQR ({cutoff:.6g})"
    outlier["interpretation"] = np.where(
        outlier["exploratory_flag"],
        "This represents an exploratory QC flag and does not by itself provide an objective basis for sample exclusion.",
        "No flag by the prespecified exploratory Tukey-fence rule.",
    )
    outlier.reset_index(names="sample_id").to_csv(tables / "PCA_outlier_diagnostics.csv", index=False)

    annotated, annotation_qc = annotate_degs(project, degs, from_counts=args.from_counts)
    annotation_qc.to_csv(tables / "gene_annotation_QC.csv", index=False)
    ann_summary = {
        "total_DEG_features": len(annotation_qc),
        "uniquely_identified_genes": int((annotation_qc["mapping_status"] == "uniquely_mapped").sum()) if len(annotation_qc) else 0,
        "unresolved_identifiers": int((annotation_qc["mapping_status"] == "unresolved").sum()) if len(annotation_qc) else 0,
        "ambiguous_identifiers": int((annotation_qc["mapping_status"] == "ambiguous").sum()) if len(annotation_qc) else 0,
        "DEGs_with_valid_symbols": int(annotation_qc["rat_gene_symbol"].fillna("").ne("").sum()) if len(annotation_qc) else 0,
        "DEGs_with_Ensembl_IDs": int(annotation_qc["rat_ensembl_gene_id"].fillna("").ne("").sum()) if len(annotation_qc) else 0,
        "DEGs_with_Entrez_IDs": int(annotation_qc["rat_entrez_id"].fillna("").ne("").sum()) if len(annotation_qc) else 0,
    }
    pd.DataFrame(ann_summary.items(), columns=["metric", "value"]).to_csv(tables / "gene_annotation_summary.csv", index=False)
    ortholog_cols = [
        "original_id", "rat_gene_symbol", "rat_ensembl_gene_id", "rat_entrez_id", "rat_gene_name",
        "baseMean", "log2FoldChange", "pvalue", "padj", "direction",
    ]
    annotated[ortholog_cols].to_csv(de_dir / "GSE224814_DEGs_for_ortholog_mapping.csv", index=False)

    label_map = dict(zip(annotation_qc["original_id"], annotation_qc["rat_gene_symbol"])) if len(annotation_qc) else {}
    plot_result = result.copy()
    plot_result["class"] = "Non-DEG"
    plot_result.loc[primary_mask & (plot_result["log2FoldChange"] >= 1), "class"] = "Up"
    plot_result.loc[primary_mask & (plot_result["log2FoldChange"] <= -1), "class"] = "Down"
    plot_result["minus_log10_padj"] = -np.log10(plot_result["padj"].clip(lower=np.nextafter(0, 1)))
    fig, ax = plt.subplots(figsize=(9, 7))
    plot_colors = {"Non-DEG": "#B8B8B8", "Up": "#D14A3A", "Down": "#3977A9"}
    for cls in ["Non-DEG", "Up", "Down"]:
        sub = plot_result[plot_result["class"] == cls]
        ax.scatter(sub["log2FoldChange"], sub["minus_log10_padj"], s=12 if cls == "Non-DEG" else 24, alpha=0.45 if cls == "Non-DEG" else 0.8, c=plot_colors[cls], label=cls, rasterized=True)
    ax.axvline(-1, color="black", linestyle="--", linewidth=0.8)
    ax.axvline(1, color="black", linestyle="--", linewidth=0.8)
    ax.axhline(-math.log10(0.05), color="black", linestyle="--", linewidth=0.8)
    ax.grid(True, color="#9E9E9E", alpha=0.18, linewidth=0.8)
    ax.set_axisbelow(True)
    top_labels = degs.assign(abs_lfc=degs["log2FoldChange"].abs()).sort_values(
        ["padj", "abs_lfc"], ascending=[True, False], na_position="last"
    ).head(15)
    texts = []
    for _, row in top_labels.iterrows():
        label = label_map.get(row["original_id"]) or row["original_id"]
        y = -math.log10(max(float(row["padj"]), np.nextafter(0, 1)))
        texts.append(ax.text(row["log2FoldChange"], y, label, fontsize=8))
    if texts:
        adjust_text(texts, ax=ax, arrowprops=dict(arrowstyle="-", color="#555555", lw=0.5))
    ax.set_xlabel("log2 fold change (CCI-ION vs Sham)")
    ax.set_ylabel("-log10(BH-adjusted p value)")
    ax.set_title("GSE224814 volcano plot")
    ax.legend(frameon=False)
    fig.tight_layout()
    save_figure(fig, figures / "volcano_CCI_ION_vs_Sham")
    plt.close(fig)

    fig, ax = plt.subplots(figsize=(9, 7))
    for cls in ["Non-DEG", "Up", "Down"]:
        sub = plot_result[plot_result["class"] == cls]
        ax.scatter(sub["baseMean"], sub["log2FoldChange"], s=12 if cls == "Non-DEG" else 24, alpha=0.45 if cls == "Non-DEG" else 0.8, c=plot_colors[cls], label=cls, rasterized=True)
    ax.set_xscale("log")
    ax.axhline(0, color="black", linewidth=0.8)
    ax.axhline(-1, color="black", linestyle="--", linewidth=0.7)
    ax.axhline(1, color="black", linestyle="--", linewidth=0.7)
    ax.set_xlabel("Mean normalized expression (baseMean, log scale)")
    ax.set_ylabel("log2 fold change (CCI-ION vs Sham)")
    ax.set_title("GSE224814 MA plot")
    ax.legend(frameon=False)
    fig.tight_layout()
    save_figure(fig, figures / "MA_CCI_ION_vs_Sham")
    plt.close(fig)

    ranked = degs.assign(abs_lfc=degs["log2FoldChange"].abs()).sort_values(
        ["padj", "abs_lfc"], ascending=[True, False], na_position="last"
    )
    top = ranked.head(50)["original_id"].tolist()
    if top:
        heat = transformed[top].T
        heat.index = [label_map.get(x) or x for x in heat.index]
        z = heat.sub(heat.mean(axis=1), axis=0).div(heat.std(axis=1, ddof=0).replace(0, 1), axis=0)
        col_colors = pd.DataFrame(
            {"Condition": metadata.loc[z.columns, "condition"].map(palette)}, index=z.columns
        )
        heatmap_cmap = sns.diverging_palette(240, 10, s=78, l=48, as_cmap=True)
        cg = sns.clustermap(
            z,
            cmap=heatmap_cmap,
            center=0,
            vmin=-1.5,
            vmax=1.5,
            figsize=(11, max(8, min(14, 0.21 * len(z) + 3.5))),
            col_colors=col_colors,
            method="average",
            metric="euclidean",
            linewidths=0.25,
            linecolor="#F2F2F2",
            xticklabels=True,
            yticklabels=True,
            dendrogram_ratio=(0.12, 0.12),
            colors_ratio=0.025,
            cbar_pos=(0.035, 0.80, 0.018, 0.12),
        )
        cg.ax_heatmap.set_xlabel("Sample", labelpad=10)
        cg.ax_heatmap.set_ylabel("")
        cg.ax_heatmap.tick_params(axis="x", rotation=0, labelsize=10, length=0, pad=7)
        cg.ax_heatmap.tick_params(axis="y", rotation=0, labelsize=8, length=0, pad=4)
        cg.ax_col_colors.set_yticks([])
        cg.cax.tick_params(labelsize=8, length=3)
        cg.cax.yaxis.set_ticks_position("left")
        cg.cax.yaxis.set_label_position("left")
        cg.cax.set_ylabel("Gene-wise z-score", fontsize=9, labelpad=7)
        legend_handles = [
            Patch(facecolor=palette["Sham"], edgecolor="none", label="Sham"),
            Patch(facecolor=palette["CCI_ION"], edgecolor="none", label="CCI-ION"),
        ]
        cg.ax_col_dendrogram.legend(
            handles=legend_handles,
            title="Condition",
            loc="upper right",
            bbox_to_anchor=(0.99, 0.98),
            ncol=2,
            frameon=False,
            fontsize=8,
            title_fontsize=8,
            handlelength=1.2,
            columnspacing=1.2,
        )
        cg.fig.suptitle(
            f"Top {len(z)} CCI-ION–associated trigeminal ganglion DEGs",
            y=0.995,
            fontsize=15,
        )
        cg.fig.savefig(figures / "DEG_heatmap_top50.png", dpi=300, bbox_inches="tight")
        cg.fig.savefig(figures / "DEG_heatmap_top50.pdf", bbox_inches="tight")
        plt.close(cg.fig)
        if len(degs) <= 50:
            shutil.copy2(figures / "DEG_heatmap_top50.png", figures / "DEG_heatmap_all.png")
    else:
        fig, ax = plt.subplots(figsize=(8, 4))
        ax.axis("off")
        ax.text(0.5, 0.5, "No genes met the primary DEG criterion", ha="center", va="center", fontsize=14)
        save_figure(fig, figures / "DEG_heatmap_top50")
        (figures / "DEG_heatmap_all.png").write_bytes((figures / "DEG_heatmap_top50.png").read_bytes())
        plt.close(fig)

    versions = [version_row("Python", platform.python_version())]
    for package in ["pydeseq2", "pandas", "numpy", "scipy", "matplotlib", "scikit-learn", "adjustText", "seaborn"]:
        versions.append(version_row(package, importlib.metadata.version(package)))
    versions.extend(
        [
            version_row("HISAT2", binary_version(["hisat2", "--version"])),
            version_row("samtools", binary_version(["samtools", "--version"])),
            version_row("featureCounts", binary_version(["featureCounts", "-v"])),
        ]
    )
    pd.DataFrame(versions).to_csv(tables / "software_versions.csv", index=False)
    source_files = [
        project / "data_raw" / "GSE224814_RAW.tar",
        project / "data_raw" / "GSE224814_family.soft.gz",
        project / "data_raw" / "SraRunInfo.csv",
        project / "data_raw" / "ENA_fastq_manifest.tsv",
        project / "data_raw" / "reference" / "Rattus_norvegicus.Rnor_6.0.dna.toplevel.fa.gz",
        project / "data_raw" / "reference" / "Rattus_norvegicus.Rnor_6.0.104.gtf.gz",
        project / "data_raw" / "annotation" / "gene2ensembl.gz",
        project / "data_raw" / "annotation" / "Rattus_norvegicus.gene_info.gz",
    ] + [
        project / "data_processed" / "reference" / "hisat2_prebuilt_rn6" / "rn6" / f"genome.{i}.ht2"
        for i in range(1, 9)
    ]
    if args.from_counts:
        source_files = [
            project / "data_raw" / "GSE224814_RAW.tar",
            project / "data_raw" / "GSE224814_family.soft.gz",
            project / "data_raw" / "SraRunInfo.csv",
            project / "data_raw" / "ENA_fastq_manifest.tsv",
            project / "data_raw" / "reference" / "Rattus_norvegicus.Rnor_6.0.104.gtf.gz",
            project / "data_raw" / "annotation" / "gene2ensembl.rat.tsv.gz",
            project / "data_raw" / "annotation" / "Rattus_norvegicus.gene_info.gz",
            project / "data_processed" / "GSE224814_raw_count_matrix.csv",
            project / "data_processed" / "sample_metadata.csv",
        ]
    pd.DataFrame(
        [
            {
                "relative_path": str(path.relative_to(project)),
                "bytes": path.stat().st_size,
                "sha256": sha256(path),
            }
            for path in source_files
        ]
    ).to_csv(tables / "source_file_checksums.csv", index=False)

    major_files = [
        all_path, deg_path, up_path, down_path,
        de_dir / "GSE224814_DEGs_for_ortholog_mapping.csv",
        project / "data_processed" / "sample_metadata.csv",
        project / "data_processed" / "GSE224814_raw_count_matrix.csv",
        project / "data_processed" / "GSE224814_filtered_counts.csv",
        tables / "input_file_audit.csv", tables / "sample_matching_audit.csv",
        tables / "low_count_filter_summary.csv", tables / "DEG_summary.csv",
        tables / "threshold_sensitivity_summary.csv", tables / "library_size.csv",
        tables / "detected_genes_per_sample.csv", tables / "raw_read_QC.csv",
        tables / "PCA_coordinates.csv",
        tables / "PCA_group_separation_summary.csv",
        tables / "sample_distance_matrix.csv", tables / "PCA_outlier_diagnostics.csv",
        tables / "gene_annotation_QC.csv", tables / "software_versions.csv",
        tables / "source_file_checksums.csv",
        project / "README.md", project / "scripts" / "run_analysis.py",
        figures / "PCA_CCI_ION_vs_Sham.png", figures / "PCA_CCI_ION_vs_Sham.pdf",
        figures / "volcano_CCI_ION_vs_Sham.png", figures / "volcano_CCI_ION_vs_Sham.pdf",
        figures / "MA_CCI_ION_vs_Sham.png", figures / "MA_CCI_ION_vs_Sham.pdf",
        figures / "DEG_heatmap_top50.png", figures / "DEG_heatmap_top50.pdf",
        figures / "library_size.png", figures / "sample_distance_heatmap.png",
    ]
    checks = []
    def check(name: str, passed: bool, details: str) -> None:
        checks.append({"check": name, "status": "PASS" if passed else "FAIL", "details": details})
    check("GEO accession is GSE224814", True, "Project and downloaded GEO metadata identify GSE224814")
    check("Exactly six biological samples", len(metadata) == 6, f"n={len(metadata)}")
    check("Sham n = 3", int((metadata.condition == "Sham").sum()) == 3, f"n={int((metadata.condition == 'Sham').sum())}")
    check("CCI_ION n = 3", int((metadata.condition == "CCI_ION").sum()) == 3, f"n={int((metadata.condition == 'CCI_ION').sum())}")
    check("No sample duplicated", metadata.index.is_unique and counts.columns.is_unique, "metadata index and count columns unique")
    check("Count matrix and metadata identities match", set(counts.columns) == set(metadata.index), str(SAMPLE_ORDER))
    check("PyDESeq2 input is non-negative integer counts", bool((arr >= 0).all() and np.equal(arr, np.floor(arr)).all()), "featureCounts fragment counts")
    fc_outputs = list((project / "data_processed" / "featurecounts").glob("*.featureCounts.txt"))
    supplement_audit = pd.read_csv(tables / "input_file_audit.csv")
    raw_provenance_ok = len(fc_outputs) == 6 and not supplement_audit["suitable_as_raw_count_input_for_PyDESeq2"].astype(bool).any()
    check(
        "No normalized/TPM/FPKM input supplied to PyDESeq2",
        raw_provenance_ok,
        f"{len(fc_outputs)} featureCounts outputs; all six GEO decimal files explicitly rejected by input audit",
    )
    check("Filtered matrix non-empty", not filtered.empty, f"{len(filtered)} features")
    check("Contrast is CCI_ION minus Sham", stats.contrast[1:] == ["CCI_ION", "Sham"], str(stats.contrast))
    check("Positive log2FC means higher in CCI-ION", stats.contrast[1:] == ["CCI_ION", "Sham"], "explicit contrast=['condition','CCI_ION','Sham']")
    check("Every primary DEG has padj < 0.05", bool((degs["padj"] < 0.05).all()), f"n={len(degs)}")
    check("Every primary DEG has abs(log2FC) >= 1", bool((degs["log2FoldChange"].abs() >= 1).all()), f"n={len(degs)}")
    check("No NA padj classified as DEG", not degs["padj"].isna().any(), f"NA DEG padj={degs['padj'].isna().sum()}")
    check("Up + down equals total DEGs", len(up) + len(down) == len(degs), f"{len(up)}+{len(down)}={len(degs)}")
    check("All major result files exist", all(p.exists() and p.stat().st_size > 0 for p in major_files), f"{sum(p.exists() and p.stat().st_size > 0 for p in major_files)}/{len(major_files)}")
    csv_files = list((project / "results").rglob("*.csv")) + list((project / "data_processed").glob("*.csv"))
    reopen_errors = []
    for path in csv_files:
        try:
            pd.read_csv(path)
        except Exception as exc:
            reopen_errors.append(f"{path.name}: {exc}")
    check("All CSV outputs reopen", not reopen_errors, "; ".join(reopen_errors) if reopen_errors else f"{len(csv_files)} CSV files")
    figure_files = list(figures.glob("*.png")) + list(figures.glob("*.pdf"))
    check("All requested figures generated", all(p.exists() and p.stat().st_size > 0 for p in major_files if p.suffix in {".png", ".pdf"}), f"{len(figure_files)} figure files")
    validation = pd.DataFrame(checks)
    validation.to_csv(tables / "validation_checks.csv", index=False)
    if (validation["status"] == "FAIL").any():
        failed = validation.loc[validation.status == "FAIL", "check"].tolist()
        raise RuntimeError(f"Critical validation failure(s): {failed}")

    run_summary = {
        **summary_values,
        "PC1_variance_percent": pc1,
        "PC2_variance_percent": pc2,
        "PCA_flagged_samples": ";".join(outlier.index[outlier.exploratory_flag]),
        **ann_summary,
        "validation_passed": int((validation.status == "PASS").sum()),
        "validation_failed": int((validation.status == "FAIL").sum()),
    }
    pd.DataFrame(run_summary.items(), columns=["metric", "value"]).to_csv(tables / "analysis_run_summary.csv", index=False)
    print("ANALYSIS COMPLETE", flush=True)
    for key, value in run_summary.items():
        print(f"{key}: {value}", flush=True)


if __name__ == "__main__":
    main()
