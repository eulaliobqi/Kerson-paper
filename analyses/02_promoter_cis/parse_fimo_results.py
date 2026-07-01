#!/usr/bin/env python3
"""
parse_fimo_results.py
Converte output do FIMO (MEME Suite) para matriz gene × motivo compatível
com o R script 02_plot_plantcare_heatmap.R.

Uso:
    python3 parse_fimo_results.py fimo_out_49genes/fimo.tsv
    python3 parse_fimo_results.py fimo_out_49genes/fimo.tsv --qval 0.05 --out fimo_parsed

Entrada (FIMO TSV colunas):
    motif_id | motif_alt_id | sequence_name | start | stop | strand |
    score | p-value | q-value | matched_sequence

Saída:
    fimo_parsed_49genes.tsv   (hit por linha, com gene_id resolvido)
    fimo_counts_49genes.csv   (matriz 49 genes × N motivos; input direto para o R)
    fimo_tf_families.tsv      (anotação de família por TF; input para col_annotation no R)
"""

import sys
import re
import csv
import argparse
from pathlib import Path
from collections import defaultdict

# IDs dos 49 genes na ordem cromossômica (mesma ordem do PlantCARE e domain arch)
GENE_IDS = [
    "Solyc01g005730","Solyc01g005760","Solyc01g005780","Solyc01g005990",
    "Solyc01g006550","Solyc01g008390","Solyc01g009690","Solyc01g009700",
    "Solyc01g073680","Solyc01g087510","Solyc01g098370","Solyc01g098680",
    "Solyc01g098690","Solyc01g099250","Solyc01g106500",
    "Solyc02g021770","Solyc02g072250","Solyc02g092040",
    "Solyc03g082780","Solyc03g083510","Solyc03g112680",
    "Solyc04g014400",
    "Solyc05g009990","Solyc05g054900","Solyc05g055190",
    "Solyc06g008270","Solyc06g008300","Solyc06g033920",
    "Solyc07g005150","Solyc07g008590","Solyc07g008600",
    "Solyc07g008620","Solyc07g008630","Solyc07g008640",
    "Solyc08g016270","Solyc08g077740",
    "Solyc09g005090",
    "Solyc10g007830","Solyc10g076500",
    "Solyc11g011180",
    "Solyc12g006020","Solyc12g009510","Solyc12g009520","Solyc12g013680",
    "Solyc12g042760","Solyc12g049190","Solyc12g099870",
    "Solyc12g099950","Solyc12g100030",
]

GENE_SET = set(GENE_IDS)

# Mapeamento de família TF JASPAR → categoria funcional biológica
# Fonte: famílias JASPAR2024 e literatura de TFs em plantas
TF_FAMILY_TO_CATEGORY = {
    # Luz
    "bHLH":    "Luz / Desenvolvimento",
    "MYB":     "ABA / Seca",
    "MYB-related": "ABA / Seca",
    # ABA/Seca
    "bZIP":    "ABA / Defesa",
    # JA
    "C2H2":    "JA / Desenvolvimento",
    # SA/Defesa
    "WRKY":    "SA / Defesa",
    # GA
    "GRAS":    "GA / Desenvolvimento",
    # Etileno/Defesa
    "AP2/ERF": "Etileno / Defesa",
    "ERF":     "Etileno / Defesa",
    # Desenvolvimento/Geral
    "NAC":     "Desenvolvimento",
    "BBR-BPC": "Desenvolvimento",
    "BES1/BZR1-homolog": "Brassinoesteroides",
    "B3":      "ABA / Desenvolvimento",
    "C3H":     "Desenvolvimento",
    "Dof":     "Luz / GA",
    "E2F/DP":  "Ciclo celular",
    "EIL":     "Etileno",
    "G2-like": "GA",
    "GeBP":    "Desenvolvimento",
    "GRF":     "Crescimento",
    "HSF":     "Estresse térmico",
    "LBD":     "Desenvolvimento lateral",
    "LFY":     "Floração",
    "LHY/CCA1/RevREV": "Circadiano",
    "MADS":    "Desenvolvimento floral",
    "NF-Y":    "Desenvolvimento",
    "RAV":     "ABA / Estresse",
    "TALE":    "Desenvolvimento",
    "TCP":     "Desenvolvimento",
    "Trihelix": "Luz",
    "VOZ":     "Floração",
    "WOX":     "Meristema",
    "ZF-HD":   "Desenvolvimento",
}


def resolve_gene_from_header(sequence_name: str) -> str | None:
    """Extrai gene_id da coluna sequence_name do FIMO TSV.
    Formato esperado do header do FASTA:
        >Solyc01g005730::SL4.0ch01:123-2123(+) upstream_2000bp
    FIMO usa tudo antes do primeiro espaço como sequence_name.
    """
    # Extrair a parte antes de '::' ou ':' ou espaço
    candidate = sequence_name.split("::")[0].split(":")[0].strip()
    # Remover sufixo de isoforma se presente (e.g. .1.1 → .1 → base)
    base = re.sub(r"\.\d+$", "", candidate)
    if base in GENE_SET:
        return base
    # Tentar com sufixo .1
    if candidate in GENE_SET:
        return candidate
    return None


def parse_fimo_tsv(fimo_file: Path, qval_cutoff: float) -> list[dict]:
    """Lê fimo.tsv e retorna lista de hits com gene_id resolvido."""
    rows = []
    unknown = set()
    with open(fimo_file, newline="", encoding="utf-8") as f:
        for line in f:
            if line.startswith("#") or not line.strip():
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 9:
                continue
            # Ignorar header se presente
            if parts[0] == "motif_id":
                continue
            motif_id, motif_alt, seq_name, start, stop, strand, score, pval, qval = parts[:9]
            try:
                if float(qval) > qval_cutoff:
                    continue
            except ValueError:
                continue
            gene_id = resolve_gene_from_header(seq_name)
            if gene_id is None:
                unknown.add(seq_name)
                continue
            rows.append({
                "gene_id":   gene_id,
                "motif_id":  motif_id,
                "tf_name":   motif_alt,
                "start":     int(start),
                "stop":      int(stop),
                "strand":    strand,
                "score":     float(score),
                "pvalue":    float(pval),
                "qvalue":    float(qval),
            })
    if unknown:
        n = len(unknown)
        print(f"  AVISO: {n} seq_name não mapeados para gene_id (exemplos: {sorted(unknown)[:3]})")
    return rows


def main():
    parser = argparse.ArgumentParser(description="Parse FIMO TSV → gene × motif matrix")
    parser.add_argument("fimo_tsv", help="Arquivo fimo.tsv gerado pelo MEME FIMO")
    parser.add_argument("--qval",   type=float, default=0.05,  help="FDR cutoff (default: 0.05)")
    parser.add_argument("--out",    default=None, help="Prefixo de saída (default: fimo_parsed)")
    args = parser.parse_args()

    fimo_file = Path(args.fimo_tsv)
    if not fimo_file.exists():
        sys.exit(f"ERRO: {fimo_file} não encontrado")

    out_prefix = Path(args.out) if args.out else fimo_file.parent / "fimo_parsed"
    out_tsv  = Path(str(out_prefix) + "_49genes.tsv")
    out_csv  = Path(str(out_prefix) + "_counts.csv")
    out_ann  = Path(str(out_prefix) + "_tf_families.tsv")

    print(f"Lendo: {fimo_file}")
    print(f"FDR cutoff: {args.qval}")

    rows = parse_fimo_tsv(fimo_file, args.qval)
    print(f"Hits filtrados (q ≤ {args.qval}): {len(rows)}")

    # ── TSV normalizado ──────────────────────────────────────────────────────
    with open(out_tsv, "w", newline="", encoding="utf-8") as f:
        w = csv.writer(f, delimiter="\t")
        w.writerow(["gene_id","motif_id","tf_name","start","stop","strand","score","pvalue","qvalue"])
        for r in rows:
            w.writerow([r["gene_id"], r["motif_id"], r["tf_name"],
                        r["start"], r["stop"], r["strand"],
                        f"{r['score']:.3f}", f"{r['pvalue']:.2e}", f"{r['qvalue']:.4f}"])
    print(f"TSV normalizado: {out_tsv}")

    # ── Matriz gene × TF (contagem de hits) ─────────────────────────────────
    counts = defaultdict(lambda: defaultdict(int))
    all_tfs = set()
    for r in rows:
        key = r["tf_name"] if r["tf_name"] else r["motif_id"]
        counts[r["gene_id"]][key] += 1
        all_tfs.add(key)

    all_tfs_sorted = sorted(all_tfs)
    with open(out_csv, "w", newline="", encoding="utf-8") as f:
        w = csv.writer(f)
        w.writerow(["gene"] + all_tfs_sorted)
        for gid in GENE_IDS:
            w.writerow([gid] + [counts[gid].get(tf, 0) for tf in all_tfs_sorted])
    print(f"Matriz de contagens (49 × {len(all_tfs_sorted)} TFs): {out_csv}")

    # ── Anotação de família por TF ───────────────────────────────────────────
    # Tentar inferir família a partir do motif_id JASPAR (formato MA0xxx.x)
    # e do nome do TF
    tf_family: dict[str, str] = {}
    for r in rows:
        key = r["tf_name"] if r["tf_name"] else r["motif_id"]
        if key not in tf_family:
            # Heurística simples por nome
            name_upper = key.upper()
            if "WRKY" in name_upper:
                fam = "WRKY"
            elif any(x in name_upper for x in ["ABF","ABI","HY5","GBF","OBF"]):
                fam = "bZIP"
            elif any(x in name_upper for x in ["MYB","MYC"]):
                fam = "MYB"
            elif any(x in name_upper for x in ["ERF","EBP","ORA","RAP"]):
                fam = "AP2/ERF"
            elif any(x in name_upper for x in ["PIF","BES","BIM","AIF"]):
                fam = "bHLH"
            elif any(x in name_upper for x in ["NAC","VND","NST"]):
                fam = "NAC"
            elif "DOF" in name_upper:
                fam = "Dof"
            elif any(x in name_upper for x in ["GRAS","SCR","SHR","HAM"]):
                fam = "GRAS"
            elif any(x in name_upper for x in ["TCP","CYC"]):
                fam = "TCP"
            elif any(x in name_upper for x in ["MADS","SEP","SOC"]):
                fam = "MADS"
            elif "WRK" in name_upper:
                fam = "WRKY"
            elif "HSF" in name_upper:
                fam = "HSF"
            else:
                fam = "Outros"
            tf_family[key] = fam

    with open(out_ann, "w", newline="", encoding="utf-8") as f:
        w = csv.writer(f, delimiter="\t")
        w.writerow(["tf_name","family","category"])
        for tf in all_tfs_sorted:
            fam = tf_family.get(tf, "Outros")
            cat = TF_FAMILY_TO_CATEGORY.get(fam, "Outros")
            w.writerow([tf, fam, cat])
    print(f"Anotação TF×família: {out_ann}")

    print()
    print("Próximo passo no servidor:")
    print(f"  Rscript 04_plot_fimo_heatmap.R {out_csv} {out_ann}")


if __name__ == "__main__":
    main()
