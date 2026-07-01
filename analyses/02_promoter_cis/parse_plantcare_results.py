#!/usr/bin/env python3
"""
parse_plantcare_results.py
Converte o output do PlantCARE para o formato TSV esperado pelo R script.

Uso:
    python3 parse_plantcare_results.py plantcare_results_49genes.txt

Saída:
    plantcare_parsed_49genes.tsv  (tab-sep, colunas: Sequence Signal Location Strand Seq Function)
    plantcare_counts_49genes.csv  (matriz gene × motivo, input direto para o R)

O PlantCARE pode entregar o resultado em dois formatos:
  A) Tab-separado com header: Sequence\tSignal\tLocation\tStrand\tSequence\tFunction
  B) HTML copiado/colado (colunas separadas por espaços irregulares)
  C) CSV com vírgula

O script detecta automaticamente o formato.
"""

import sys
import re
import csv
from pathlib import Path
from collections import defaultdict

# IDs dos 49 genes na ordem cromossômica
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

# Elementos cis de interesse (correspondem ao dicionário no R script)
CIS_ELEMENTS = [
    "ACE","Box 4","G-box","GT1-motif","I-box","Sp1","LAMP-element",
    "ABRE","CACGTG-motif","DRE core","MBS","MYB recognition",
    "CGTCA-motif","TGACG-motif",
    "SARE","TCA-element","W-box","TC-rich repeats","TGA-element",
    "GARE-motif","P-box","TATC-box",
    "AuxRR-core","TGA-element",
    "GCC-box","ERE",
    "CAT-box","CCAAT-box","as-2-element",
    "circadian","Evening Element",
]
CIS_SET = set(e.upper() for e in CIS_ELEMENTS)


def detect_format(lines):
    """Detecta o formato do arquivo PlantCARE."""
    for line in lines[:20]:
        line = line.strip()
        if not line or line.startswith('#'):
            continue
        if '\t' in line:
            return 'tsv'
        if ',' in line and line.count(',') >= 4:
            return 'csv'
    return 'space'  # colunas separadas por múltiplos espaços (cópia HTML)


def parse_tsv(lines):
    rows = []
    reader = csv.reader(lines, delimiter='\t')
    header = None
    for row in reader:
        if not row or not row[0].strip():
            continue
        if header is None:
            header = [c.strip().lower() for c in row]
            continue
        if len(row) >= 5:
            rows.append({
                'sequence': row[0].strip(),
                'signal':   row[1].strip(),
                'location': row[2].strip(),
                'strand':   row[3].strip(),
                'function': row[5].strip() if len(row) > 5 else '',
            })
    return rows


def parse_csv(lines):
    rows = []
    reader = csv.reader(lines)
    header = None
    for row in reader:
        if not row or not row[0].strip():
            continue
        if header is None:
            header = [c.strip().lower() for c in row]
            continue
        if len(row) >= 5:
            rows.append({
                'sequence': row[0].strip(),
                'signal':   row[1].strip(),
                'location': row[2].strip(),
                'strand':   row[3].strip(),
                'function': row[5].strip() if len(row) > 5 else '',
            })
    return rows


def parse_space(lines):
    """Parseia saída PlantCARE copiada do HTML (espaços variáveis entre colunas)."""
    rows = []
    # Padrão: gene_id  motif_name  position  strand  seq_motif  função
    # A posição é um número inteiro; strand é + ou -
    pat = re.compile(
        r'^(\S+)\s+'          # sequence id
        r'(.+?)\s{2,}'        # signal (pode ter espaços internos — separador = ≥2 espaços)
        r'(\-?\d+)\s+'        # location (inteiro, pode ser negativo)
        r'([+\-])\s+'         # strand
        r'(\S+)\s*'           # seq
        r'(.*)?$'             # function (opcional)
    )
    for line in lines:
        line = line.rstrip()
        if not line or line.startswith('#') or line.lower().startswith('sequence'):
            continue
        m = pat.match(line)
        if m:
            rows.append({
                'sequence': m.group(1),
                'signal':   m.group(2).strip(),
                'location': m.group(3),
                'strand':   m.group(4),
                'function': m.group(6).strip() if m.group(6) else '',
            })
    return rows


def main():
    infile = Path(sys.argv[1]) if len(sys.argv) > 1 else Path("plantcare_results_49genes.txt")
    if not infile.exists():
        print(f"ERRO: arquivo não encontrado: {infile}")
        sys.exit(1)

    lines = infile.read_text(encoding='utf-8', errors='replace').splitlines()
    fmt   = detect_format(lines)
    print(f"Formato detectado: {fmt}")

    if fmt == 'tsv':
        rows = parse_tsv(lines)
    elif fmt == 'csv':
        rows = parse_csv(lines)
    else:
        rows = parse_space(lines)

    if not rows:
        print("ERRO: nenhuma linha parseable encontrada. Verifique o formato do arquivo.")
        sys.exit(1)

    print(f"Total de hits: {len(rows)}")

    # ── Escrever TSV normalizado ─────────────────────────────────────────────
    out_tsv = infile.parent / "plantcare_parsed_49genes.tsv"
    with open(out_tsv, 'w', newline='', encoding='utf-8') as f:
        w = csv.writer(f, delimiter='\t')
        w.writerow(["Sequence","Signal","Location","Strand","Function"])
        for r in rows:
            w.writerow([r['sequence'], r['signal'], r['location'], r['strand'], r['function']])
    print(f"TSV normalizado: {out_tsv}")

    # ── Construir matriz gene × motivo ───────────────────────────────────────
    counts = defaultdict(lambda: defaultdict(int))
    unknown_genes = set()
    for r in rows:
        gid = r['sequence']
        # Aceitar tanto Solyc ID direto quanto Solyc ID com sufixo de isoforma
        gid_base = re.sub(r'\.\d+$', '', gid)
        if gid_base not in GENE_IDS:
            unknown_genes.add(gid)
            continue
        counts[gid_base][r['signal']] += 1

    if unknown_genes:
        print(f"AVISO: {len(unknown_genes)} IDs não reconhecidos no resultado PlantCARE:")
        for g in sorted(unknown_genes)[:10]:
            print(f"  {g}")

    # Coletar todos os motivos encontrados
    all_motifs = sorted({sig for gid in counts for sig in counts[gid]})
    print(f"Motivos únicos encontrados: {len(all_motifs)}")

    out_csv = infile.parent / "plantcare_counts_49genes.csv"
    with open(out_csv, 'w', newline='', encoding='utf-8') as f:
        w = csv.writer(f)
        w.writerow(['gene'] + all_motifs)
        for gid in GENE_IDS:
            row = [gid] + [counts[gid].get(m, 0) for m in all_motifs]
            w.writerow(row)
    print(f"Matriz de contagens (49 × {len(all_motifs)}): {out_csv}")
    print()
    print("Próximo passo no servidor:")
    print("  git pull")
    print("  cd ~/kerson-paper/analyses/02_promoter_cis")
    print("  Rscript 02_plot_plantcare_heatmap.R plantcare_results_49genes.txt")


if __name__ == "__main__":
    main()
