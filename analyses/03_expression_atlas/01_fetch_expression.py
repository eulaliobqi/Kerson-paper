#!/usr/bin/env python3
"""
Busca dados de expressão dos 7 genes LRR-RLP no TFGD (Cornell) e alternativas.
Salva expression_matrix.csv para plotagem com 02_plot_expression_heatmap.R

Fontes tentadas em ordem:
  1. TFGD – ted.bti.cornell.edu (microarray/RNA-Seq, ~50 condições)
  2. TomExpresso – tomexpress.toulouse.inra.fr
  3. SGN Expression Viewer (solgenomics.net)

Uso:
  python3 01_fetch_expression.py
  python3 01_fetch_expression.py --demo   # gera matriz simulada para testar o plot
"""

import sys, os, time, csv, json
import urllib.request
import urllib.parse
import argparse

GENES = [
    "Solyc05g055190", "Solyc03g112680", "Solyc05g009990",
    "Solyc12g042760", "Solyc02g072250", "Solyc02g092040", "Solyc10g007830"
]

GENE_LABELS = {
    "Solyc05g055190": "SlRLP1 (CLV2)",
    "Solyc03g112680":  "SlRLP2",
    "Solyc05g009990":  "SlRLP3 (RIC7)",
    "Solyc12g042760":  "SlRLP4 (TMM)",
    "Solyc02g072250":  "SlRLP5 (SNC2/3)",
    "Solyc02g092040":  "SlRLP6",
    "Solyc10g007830":  "SlRLP7",
}

OUTFILE = os.path.join(os.path.dirname(__file__), "expression_matrix.csv")


# ── Fonte 1: TFGD ─────────────────────────────────────────────────────────────
def fetch_tfgd(gene_id: str) -> dict | None:
    """Tenta buscar expressão via TFGD (Cornell). Retorna dict {condition: log2FC} ou None."""
    base = "http://ted.bti.cornell.edu/cgi-bin/TFGD/expression/input_gene_expression.cgi"
    params = urllib.parse.urlencode({"gene_id": gene_id, "submit": "Go"})
    url = f"{base}?{params}"
    try:
        req = urllib.request.Request(url, headers={"User-Agent": "Mozilla/5.0"})
        with urllib.request.urlopen(req, timeout=15) as r:
            html = r.read().decode("utf-8", errors="ignore")
        # Parse da tabela HTML (simplificado – adaptar ao formato real do TFGD)
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
    """Tenta buscar expressão via SGN (solgenomics.net)."""
    # SGN tem um viewer mas API pública limitada; URL experimental:
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
    """Gera matriz de expressão simulada para testar o pipeline de visualização."""
    import random
    random.seed(2024)

    # Condições relevantes para o artigo (seca + vírus + tecidos)
    conditions = [
        "Raiz", "Caule", "Folha jovem", "Folha adulta", "Flor", "Fruto verde",
        "Fruto maduro", "Seca TRA65%", "Seca TRA50%", "Seca TRA45%", "Seca TRA40%",
        "ToYSV 7dpi", "ToYSV 15dpi", "ToYSV 21dpi",
        "Mock vírus", "Mock seca",
    ]

    rows = []
    for gene_id, label in GENE_LABELS.items():
        row = {"gene_id": gene_id, "label": label}
        for cond in conditions:
            # Simula padrões biologicamente plausíveis
            if "Seca" in cond and gene_id in ("Solyc02g072250", "Solyc05g055190"):
                val = round(random.uniform(1.5, 4.0), 3)
            elif "ToYSV" in cond and gene_id in ("Solyc02g072250", "Solyc02g092040"):
                val = round(random.uniform(2.0, 5.0), 3)
            elif cond in ("Raiz", "Caule"):
                val = round(random.uniform(-1.0, 1.0), 3)
            else:
                val = round(random.gauss(0, 1.2), 3)
            row[cond] = val
        rows.append(row)
    return rows


# ── Main ──────────────────────────────────────────────────────────────────────
def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--demo", action="store_true",
                        help="Gera matriz simulada sem acessar internet")
    args = parser.parse_args()

    if args.demo:
        print("[DEMO] Gerando matriz simulada...")
        rows = generate_demo_matrix()
    else:
        print("Buscando dados de expressão (TFGD → SGN → DEMO fallback)...")
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
                print(f"    Sem dados online — adicionando linha vazia (completar manualmente)")
                rows.append({"gene_id": gene_id, "label": GENE_LABELS[gene_id]})
            time.sleep(1)

        if not any(len(r) > 2 for r in rows):
            print("\nNenhum dado obtido online. Gerando DEMO para testar visualização...")
            rows = generate_demo_matrix()

    # Salvar CSV
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
    print(f"\nSe os dados online estiverem vazios, preencha {OUTFILE}")
    print("  manualmente com valores de expressão do TFGD ou de RNA-Seq público (GEO).")


if __name__ == "__main__":
    main()
