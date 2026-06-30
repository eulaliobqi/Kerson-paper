#!/usr/bin/env python3
"""
Busca dados de expressão dos 49 genes LRR-RLP no TFGD (Cornell) e alternativas.
Salva expression_matrix.csv para plotagem com 02_plot_expression_heatmap.R

Fontes tentadas em ordem:
  1. TFGD – ted.bti.cornell.edu (microarray/RNA-Seq, ~50 condições)
  2. SGN Expression Viewer (solgenomics.net)

Uso:
  python3 01_fetch_expression.py
  python3 01_fetch_expression.py --demo   # gera matriz simulada para testar o plot
"""

import sys, os, time, csv, json
import urllib.request
import urllib.parse
import argparse

# ── 49 genes ordenados por cromossomo ────────────────────────────────────────
GENES = [
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

# Rótulos: 7 genes focais com nomes funcionais; demais com ID curto
GENE_LABELS = {
    "Solyc01g005730": "SlRLP01","Solyc01g005760": "SlRLP02","Solyc01g005780": "SlRLP03",
    "Solyc01g005990": "SlRLP04","Solyc01g006550": "SlRLP05","Solyc01g008390": "SlRLP06",
    "Solyc01g009690": "SlRLP07","Solyc01g009700": "SlRLP08","Solyc01g073680": "SlRLP09",
    "Solyc01g087510": "SlRLP10","Solyc01g098370": "SlRLP11","Solyc01g098680": "SlRLP12",
    "Solyc01g098690": "SlRLP13","Solyc01g099250": "SlRLP14","Solyc01g106500": "SlRLP15",
    "Solyc02g021770": "SlRLP16",
    "Solyc02g072250": "SlRLP17*(SNC2/3)",  # focal
    "Solyc02g092040": "SlRLP18*",           # focal
    "Solyc03g082780": "SlRLP19","Solyc03g083510": "SlRLP20",
    "Solyc03g112680": "SlRLP21*",           # focal
    "Solyc04g014400": "SlRLP22",
    "Solyc05g009990": "SlRLP23*(RIC7)",     # focal
    "Solyc05g054900": "SlRLP24",
    "Solyc05g055190": "SlRLP25*(CLV2)",     # focal
    "Solyc06g008270": "SlRLP26","Solyc06g008300": "SlRLP27","Solyc06g033920": "SlRLP28",
    "Solyc07g005150": "SlRLP29","Solyc07g008590": "SlRLP30","Solyc07g008600": "SlRLP31",
    "Solyc07g008620": "SlRLP32","Solyc07g008630": "SlRLP33","Solyc07g008640": "SlRLP34",
    "Solyc08g016270": "SlRLP35","Solyc08g077740": "SlRLP36",
    "Solyc09g005090": "SlRLP37",
    "Solyc10g007830": "SlRLP38*(ToYSV/Seca)",  # focal
    "Solyc10g076500": "SlRLP39",
    "Solyc11g011180": "SlRLP40",
    "Solyc12g006020": "SlRLP41","Solyc12g009510": "SlRLP42","Solyc12g009520": "SlRLP43",
    "Solyc12g013680": "SlRLP44",
    "Solyc12g042760": "SlRLP45*(TMM)",      # focal
    "Solyc12g049190": "SlRLP46","Solyc12g099870": "SlRLP47",
    "Solyc12g099950": "SlRLP48","Solyc12g100030": "SlRLP49",
}

OUTFILE = os.path.join(os.path.dirname(__file__), "expression_matrix.csv")


# ── Fonte 1: TFGD ─────────────────────────────────────────────────────────────
def fetch_tfgd(gene_id: str) -> dict | None:
    base   = "http://ted.bti.cornell.edu/cgi-bin/TFGD/expression/input_gene_expression.cgi"
    params = urllib.parse.urlencode({"gene_id": gene_id, "submit": "Go"})
    url    = f"{base}?{params}"
    try:
        req = urllib.request.Request(url, headers={"User-Agent": "Mozilla/5.0"})
        with urllib.request.urlopen(req, timeout=15) as r:
            html = r.read().decode("utf-8", errors="ignore")
        import re
        rows = re.findall(r'<tr[^>]*>(.*?)</tr>', html, re.DOTALL)
        data = {}
        for row in rows:
            cells = re.findall(r'<t[dh][^>]*>(.*?)</t[dh]>', row, re.DOTALL)
            cells = [re.sub(r'<.*?>', '', c).strip() for c in cells]
            if len(cells) >= 2:
                try:
                    data[cells[0]] = float(cells[1])
                except ValueError:
                    pass
        return data if data else None
    except Exception as e:
        print(f"  TFGD falhou para {gene_id}: {e}")
        return None


# ── Fonte 2: SGN Expression API ───────────────────────────────────────────────
def fetch_sgn(gene_id: str) -> dict | None:
    url = f"https://solgenomics.net/api/v1/expression/{gene_id}"
    try:
        req = urllib.request.Request(url, headers={"Accept": "application/json",
                                                    "User-Agent": "Mozilla/5.0"})
        with urllib.request.urlopen(req, timeout=15) as r:
            data = json.loads(r.read())
        if isinstance(data, dict) and "expression" in data:
            return data["expression"]
        return None
    except Exception as e:
        print(f"  SGN API falhou para {gene_id}: {e}")
        return None


# ── Modo DEMO ─────────────────────────────────────────────────────────────────
def generate_demo_matrix() -> list[dict]:
    """Gera matriz de expressão simulada biologicamente plausível para 49 genes."""
    import random
    random.seed(2024)

    conditions = [
        "Raiz", "Caule", "Folha jovem", "Folha adulta", "Flor", "Fruto verde",
        "Fruto maduro", "Seca TRA65%", "Seca TRA50%", "Seca TRA45%", "Seca TRA40%",
        "ToYSV 7dpi", "ToYSV 15dpi", "ToYSV 21dpi",
        "Mock vírus", "Mock seca",
    ]

    # Grupos funcionais com padrões de expressão distintos
    group_patterns = {
        # Genes focais — alta expressão em seca e vírus
        "stress_focal": {
            "genes": {"Solyc02g072250","Solyc05g055190","Solyc10g007830"},
            "Seca": (2.0, 4.5), "ToYSV": (2.5, 5.0), "default": (-0.5, 1.5),
        },
        "virus_focal": {
            "genes": {"Solyc02g092040","Solyc03g112680"},
            "Seca": (0.5, 2.0), "ToYSV": (3.0, 6.0), "default": (-1.0, 1.0),
        },
        "dev_focal": {
            "genes": {"Solyc05g009990","Solyc12g042760"},
            "Flor": (1.5, 3.5), "Fruto": (2.0, 4.0), "default": (-0.5, 1.0),
        },
        # Chr01 cluster — expressão constitutiva + seca moderada
        "chr01_cluster": {
            "genes": {g for g in GENES if g.startswith("Solyc01")},
            "Seca": (0.5, 2.0), "ToYSV": (0.2, 1.5), "default": (0.0, 2.0),
        },
        # Chr07 cluster — defesa
        "chr07_cluster": {
            "genes": {"Solyc07g008590","Solyc07g008600","Solyc07g008620",
                      "Solyc07g008630","Solyc07g008640"},
            "ToYSV": (1.5, 4.0), "Seca": (0.2, 1.5), "default": (-0.5, 1.5),
        },
        # Chr12 cluster — fruto
        "chr12_cluster": {
            "genes": {"Solyc12g006020","Solyc12g009510","Solyc12g009520",
                      "Solyc12g099870","Solyc12g099950","Solyc12g100030"},
            "Fruto": (2.5, 5.0), "default": (-1.0, 0.5),
        },
    }

    def get_pattern(gene_id, cond):
        for pat in group_patterns.values():
            if gene_id in pat["genes"]:
                for key, rng in pat.items():
                    if key == "genes":
                        continue
                    if key in cond:
                        return rng
                return pat.get("default", (-1.0, 1.5))
        return (-0.5, 1.5)  # fallback genérico

    rows = []
    for gene_id in GENES:
        label = GENE_LABELS[gene_id]
        row = {"gene_id": gene_id, "label": label}
        for cond in conditions:
            lo, hi = get_pattern(gene_id, cond)
            row[cond] = round(random.uniform(lo, hi), 3)
        rows.append(row)
    return rows


# ── Main ──────────────────────────────────────────────────────────────────────
def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--demo", action="store_true",
                        help="Gera matriz simulada sem acessar internet")
    args = parser.parse_args()

    if args.demo:
        print("[DEMO] Gerando matriz simulada para 49 genes...")
        rows = generate_demo_matrix()
    else:
        print(f"Buscando dados de expressão para {len(GENES)} genes (TFGD → SGN → DEMO)...")
        rows = []
        for gene_id in GENES:
            print(f"  {gene_id}...")
            data = fetch_tfgd(gene_id) or fetch_sgn(gene_id)
            if data:
                row = {"gene_id": gene_id, "label": GENE_LABELS[gene_id]}
                row.update(data)
                rows.append(row)
                print(f"    OK: {len(data)} condições")
            else:
                print(f"    Sem dados online — adicionando linha vazia")
                rows.append({"gene_id": gene_id, "label": GENE_LABELS[gene_id]})
            time.sleep(0.5)

        if not any(len(r) > 2 for r in rows):
            print("\nNenhum dado obtido online. Gerando DEMO para testar visualização...")
            rows = generate_demo_matrix()

    all_keys = sorted(set(k for r in rows for k in r if k not in ("gene_id","label")))
    fieldnames = ["gene_id", "label"] + all_keys

    with open(OUTFILE, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames, extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)

    print(f"\nMatriz salva: {OUTFILE}")
    print(f"Genes: {len(rows)} | Condições: {len(all_keys)}")
    print(f"\nPróximo passo:")
    print(f"  Rscript 02_plot_expression_heatmap.R")


if __name__ == "__main__":
    main()
