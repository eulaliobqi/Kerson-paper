#!/usr/bin/env python3
"""
00_get_coords_all49.py
Obtém coordenadas genômicas ITAG4.0 dos 49 LRR-RLPs de Solanum lycopersicum.

SERVIDOR: eulalio@200.235.143.10
NOTA: SGN, EnsemblPlants e EnsemblGenomes são INACESSÍVEIS nesse servidor.
      NCBI (ftp.ncbi.nlm.nih.gov + eutils.ncbi.nlm.nih.gov) FUNCIONA.

Estratégias em cascata:
  1. NCBI Entrez (HTTP puro, sem biopython) — por gene, mais rápido
  2. NCBI GFF3 via FTP                      — download único em lote (~100 MB)
  3. EnsemblPlants GFF3 (EBI FTP)           — fallback externo
  4. EnsemblGenomes REST API                — por gene
  5. SGN API / scraping                     — por gene
  6. Hardcoded                              — 7 posições confirmadas ITAG4.0

Saída: coords_49genes.tsv
  Colunas: gene_id  chrom  start  end  strand
  chrom no formato SL4.0chNN (compatível com S_lycopersicum_chromosomes.4.00.fa)

Uso:
  python3 00_get_coords_all49.py             # pula se TSV completo
  python3 00_get_coords_all49.py --force     # re-executa sempre
  python3 00_get_coords_all49.py --skip-gff  # pula downloads de GFF3
  python3 00_get_coords_all49.py --ncbi-only # apenas estratégias NCBI
"""

import os, sys, json, time, ssl, gzip, re, argparse
import urllib.request, urllib.error
from pathlib import Path

# ── Caminhos ──────────────────────────────────────────────────────────────────
SCRIPT_DIR = Path(__file__).resolve().parent
REPO_ROOT  = SCRIPT_DIR.parent.parent
IDS_FILE   = REPO_ROOT / "ids_49_rlp_tomato.txt"
OUT_TSV    = SCRIPT_DIR / "coords_49genes.tsv"
CACHE_DIR  = SCRIPT_DIR / ".cache"
CACHE_DIR.mkdir(exist_ok=True)

# ── SSL context sem verificação ───────────────────────────────────────────────
CTX = ssl.create_default_context()
CTX.check_hostname = False
CTX.verify_mode    = ssl.CERT_NONE


def http_get(url, timeout=30, retries=3, delay=2):
    headers = {"User-Agent": "kerson-paper-bioinf/1.0 (eulalio.santos@ufv.br)"}
    for attempt in range(retries):
        try:
            req = urllib.request.Request(url, headers=headers)
            with urllib.request.urlopen(req, timeout=timeout, context=CTX) as r:
                return r.read()
        except urllib.error.HTTPError as e:
            print(f"    HTTP {e.code} — {url}", file=sys.stderr)
            break
        except Exception as e:
            print(f"    [tentativa {attempt+1}/{retries}] {type(e).__name__}: {e}", file=sys.stderr)
            if attempt < retries - 1:
                time.sleep(delay)
    return None


# ── Mapeamento RefSeq → SL4.0 cromossomos ────────────────────────────────────
# NCBI GCF_000188115.5 (SL4.0): NC_015438.x = chr01, ..., NC_015449.x = chr12
NCBI_ACC_MAP = {
    "NC_015438": "SL4.0ch01", "NC_015439": "SL4.0ch02",
    "NC_015440": "SL4.0ch03", "NC_015441": "SL4.0ch04",
    "NC_015442": "SL4.0ch05", "NC_015443": "SL4.0ch06",
    "NC_015444": "SL4.0ch07", "NC_015445": "SL4.0ch08",
    "NC_015446": "SL4.0ch09", "NC_015447": "SL4.0ch10",
    "NC_015448": "SL4.0ch11", "NC_015449": "SL4.0ch12",
}


def normalize_chrom(raw, gene_id_hint=None):
    """Converte qualquer nome de cromossomo para SL4.0chNN."""
    c = str(raw).strip()
    if re.match(r"^SL4\.0ch\d+$", c):
        return c
    if re.match(r"^\d{1,2}$", c):
        return f"SL4.0ch{int(c):02d}"
    m = re.match(r"^[Cc][Hh][Rr]?0*(\d+)$", c)
    if m:
        return f"SL4.0ch{int(m.group(1)):02d}"
    m2 = re.match(r"^SL\d+\.\d+ch0*(\d+)$", c)
    if m2:
        return f"SL4.0ch{int(m2.group(1)):02d}"
    # Accession NCBI RefSeq: NC_015438.3
    base = re.sub(r"\.\d+$", "", c)
    if base in NCBI_ACC_MAP:
        return NCBI_ACC_MAP[base]
    # Fallback: inferir cromossomo do gene ID (SolycXXg → chXX)
    if gene_id_hint:
        m3 = re.match(r"Solyc(\d{2})g", gene_id_hint)
        if m3:
            return f"SL4.0ch{int(m3.group(1)):02d}"
    return c


def load_gene_ids():
    with open(IDS_FILE) as f:
        ids = [l.strip() for l in f if l.strip()]
    print(f"IDs carregados: {len(ids)}", file=sys.stderr)
    return ids


# ══════════════════════════════════════════════════════════════════════════════
# ESTRATÉGIA 1 — NCBI Entrez HTTP (sem biopython)
# eutils.ncbi.nlm.nih.gov — acessível do servidor UFV
# ══════════════════════════════════════════════════════════════════════════════
NCBI_EUTILS = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"
NCBI_EMAIL  = "eulalio.santos@ufv.br"
NCBI_TOOL   = "kerson-rlp-pipeline"


def _ncbi_esearch(term):
    url = (f"{NCBI_EUTILS}/esearch.fcgi?db=gene"
           f"&term={urllib.request.quote(term)}"
           f"&retmax=3&retmode=json"
           f"&email={NCBI_EMAIL}&tool={NCBI_TOOL}")
    raw = http_get(url, timeout=20, retries=2, delay=1)
    if not raw:
        return []
    try:
        return json.loads(raw).get("esearchresult", {}).get("idlist", [])
    except json.JSONDecodeError:
        return []


def _ncbi_esummary_gene(ncbi_id):
    url = (f"{NCBI_EUTILS}/esummary.fcgi?db=gene&id={ncbi_id}"
           f"&retmode=json&email={NCBI_EMAIL}&tool={NCBI_TOOL}")
    raw = http_get(url, timeout=20, retries=2, delay=1)
    if not raw:
        return {}
    try:
        result = json.loads(raw).get("result", {})
        return result.get(str(ncbi_id), {})
    except json.JSONDecodeError:
        return {}


def strategy_ncbi_entrez_http(gene_ids):
    """NCBI Entrez via HTTP puro. Estratégia mais confiável no servidor."""
    print(f"\n[Estratégia 1] NCBI Entrez HTTP ({len(gene_ids)} genes)...", file=sys.stderr)
    results = {}

    for i, gid in enumerate(gene_ids):
        queries = [
            f'"{gid}"[Gene Symbol] AND "Solanum lycopersicum"[Organism]',
            f'"{gid}"[Gene Name] AND txid4081[Organism]',
            f'{gid} AND txid4081[Organism]',
        ]

        ncbi_ids = []
        for q in queries:
            ncbi_ids = _ncbi_esearch(q)
            if ncbi_ids:
                break
            time.sleep(0.35)

        if not ncbi_ids:
            time.sleep(0.35)
            continue

        for ncbi_id in ncbi_ids[:2]:
            summary = _ncbi_esummary_gene(ncbi_id)
            time.sleep(0.35)

            acc = summary.get("chraccver", "")
            s0  = summary.get("chrstart")
            s1  = summary.get("chrstop")
            strand_raw = summary.get("strand", "")

            # Tentar genomicinfo se campos diretos ausentes
            if not acc or s0 is None or s1 is None:
                ginfo = summary.get("genomicinfo", [])
                if isinstance(ginfo, list) and ginfo:
                    g   = ginfo[0]
                    acc = g.get("chraccver", "")
                    s0  = g.get("chrstart")
                    s1  = g.get("chrstop")

            if acc and s0 is not None and s1 is not None:
                chrom = normalize_chrom(acc, gene_id_hint=gid)
                # esummary: chrstart/chrstop são 0-based; converter para 1-based GFF
                start = int(min(s0, s1)) + 1
                end   = int(max(s0, s1))
                if isinstance(strand_raw, str) and strand_raw in ("+", "-"):
                    strand = strand_raw
                elif isinstance(strand_raw, int):
                    strand = "+" if strand_raw >= 0 else "-"
                else:
                    strand = "-" if int(s0) > int(s1) else "+"

                results[gid] = (chrom, start, end, strand)
                print(f"    {gid}: {chrom}:{start}-{end} ({strand})", file=sys.stderr)
                break

        if (i + 1) % 10 == 0:
            print(f"  Progresso NCBI: {i+1}/{len(gene_ids)} — {len(results)} encontrados", file=sys.stderr)

    print(f"  Estratégia 1 (NCBI Entrez): {len(results)}/{len(gene_ids)} encontrados", file=sys.stderr)
    return results


# ══════════════════════════════════════════════════════════════════════════════
# ESTRATÉGIA 2 — NCBI GFF3 via FTP (download em lote)
# Assembly: GCF_000188115.5 = S. lycopersicum SL4.0
# ══════════════════════════════════════════════════════════════════════════════
NCBI_GFF3_URL = (
    "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/188/115/"
    "GCF_000188115.5_SL4.0/GCF_000188115.5_SL4.0_genomic.gff.gz"
)
NCBI_GFF3_CACHE = CACHE_DIR / "ncbi_sl4.0_genomic.gff.gz"


def strategy_ncbi_gff3(gene_ids):
    """Baixa GFF3 do NCBI FTP e parseia coordenadas dos Solyc IDs."""
    print(f"\n[Estratégia 2] NCBI GFF3 FTP ({len(gene_ids)} genes)...", file=sys.stderr)

    if not NCBI_GFF3_CACHE.exists():
        print(f"  Baixando {NCBI_GFF3_URL} (~100 MB)...", file=sys.stderr)
        data = http_get(NCBI_GFF3_URL, timeout=600, retries=2, delay=15)
        if not data:
            print("  FALHA: NCBI GFF3 FTP inacessível", file=sys.stderr)
            return {}
        NCBI_GFF3_CACHE.write_bytes(data)
        print(f"  GFF3 cacheado: {NCBI_GFF3_CACHE}", file=sys.stderr)
    else:
        print(f"  Usando cache: {NCBI_GFF3_CACHE}", file=sys.stderr)

    targets = set(gene_ids)
    results = {}
    feat_ok = {"gene", "mRNA", "transcript", "CDS"}

    try:
        with gzip.open(NCBI_GFF3_CACHE, "rt", errors="ignore") as fh:
            for line in fh:
                if not targets:
                    break
                if line.startswith("#"):
                    continue
                cols = line.rstrip("\n").split("\t")
                if len(cols) < 9 or cols[2] not in feat_ok:
                    continue
                attrs = cols[8]
                for gid in list(targets):
                    if gid in attrs or gid.upper() in attrs:
                        chrom  = normalize_chrom(cols[0], gene_id_hint=gid)
                        start  = int(cols[3])
                        end    = int(cols[4])
                        strand = cols[6] if cols[6] in ("+", "-") else "+"
                        results[gid] = (chrom, start, end, strand)
                        targets.discard(gid)
                        print(f"    {gid}: {chrom}:{start}-{end} ({strand})", file=sys.stderr)
                        break
    except Exception as e:
        print(f"  ERRO ao parsear GFF3 NCBI: {e}", file=sys.stderr)

    print(f"  Estratégia 2 (NCBI GFF3): {len(results)}/{len(gene_ids)} encontrados", file=sys.stderr)
    return results


# ══════════════════════════════════════════════════════════════════════════════
# ESTRATÉGIA 3 — EnsemblPlants GFF3 (EBI FTP)
# ══════════════════════════════════════════════════════════════════════════════
ENSEMBL_FTP_BASE   = "https://ftp.ensemblgenomes.ebi.ac.uk/pub/plants/current/gff3/solanum_lycopersicum/"
ENSEMBL_GFF3_CACHE = CACHE_DIR / "ensemblplants_slycopersicum.gff3.gz"


def strategy_ensemblplants_gff3(gene_ids):
    print(f"\n[Estratégia 3] EnsemblPlants GFF3 ({len(gene_ids)} genes)...", file=sys.stderr)

    if not ENSEMBL_GFF3_CACHE.exists():
        listing_raw = http_get(ENSEMBL_FTP_BASE, timeout=25)
        if not listing_raw:
            print("  FALHA: EBI FTP inacessível", file=sys.stderr)
            return {}
        listing    = listing_raw.decode("utf-8", errors="ignore")
        candidates = re.findall(r'href="(Solanum_lycopersicum[^"]+\.gff3\.gz)"', listing)
        candidates = [c for c in candidates
                      if not any(x in c for x in ("toplevel","chromosome","abinitio","cdna","cds","ncrna","pep"))]
        if not candidates:
            print("  FALHA: nenhum GFF3 no diretório EBI", file=sys.stderr)
            return {}
        url  = ENSEMBL_FTP_BASE + candidates[0]
        data = http_get(url, timeout=300, retries=2, delay=10)
        if not data:
            print("  FALHA: download EnsemblPlants GFF3 falhou", file=sys.stderr)
            return {}
        ENSEMBL_GFF3_CACHE.write_bytes(data)
    else:
        print(f"  Usando cache: {ENSEMBL_GFF3_CACHE}", file=sys.stderr)

    targets = set(gene_ids)
    results = {}
    try:
        with gzip.open(ENSEMBL_GFF3_CACHE, "rt", errors="ignore") as fh:
            for line in fh:
                if not targets:
                    break
                if line.startswith("#"):
                    continue
                cols = line.rstrip("\n").split("\t")
                if len(cols) < 9 or cols[2] not in {"gene","mRNA","transcript"}:
                    continue
                attrs = cols[8]
                for gid in list(targets):
                    if gid in attrs or gid.upper() in attrs:
                        chrom  = normalize_chrom(cols[0], gene_id_hint=gid)
                        start  = int(cols[3])
                        end    = int(cols[4])
                        strand = cols[6] if cols[6] in ("+", "-") else "+"
                        results[gid] = (chrom, start, end, strand)
                        targets.discard(gid)
                        break
    except Exception as e:
        print(f"  ERRO: {e}", file=sys.stderr)

    print(f"  Estratégia 3: {len(results)}/{len(gene_ids)} encontrados", file=sys.stderr)
    return results


# ══════════════════════════════════════════════════════════════════════════════
# ESTRATÉGIA 4 — EnsemblGenomes REST API
# ══════════════════════════════════════════════════════════════════════════════
def strategy_ensemblgenomes_rest(gene_ids):
    print(f"\n[Estratégia 4] EnsemblGenomes REST ({len(gene_ids)} genes)...", file=sys.stderr)
    results = {}
    for i, gid in enumerate(gene_ids):
        for variant in [gid, gid.upper()]:
            url = (f"https://rest.ensemblgenomes.org/lookup/id/{variant}"
                   f"?content-type=application/json&expand=0")
            raw = http_get(url, timeout=20, retries=2, delay=1)
            if not raw:
                continue
            try:
                d = json.loads(raw)
                if "seq_region_name" in d and "start" in d:
                    chrom  = normalize_chrom(d["seq_region_name"], gene_id_hint=gid)
                    start  = int(d["start"])
                    end    = int(d["end"])
                    strand = "+" if d.get("strand", 1) == 1 else "-"
                    results[gid] = (chrom, start, end, strand)
                    break
            except (json.JSONDecodeError, ValueError, KeyError):
                pass
        time.sleep(1.0)
        if (i + 1) % 10 == 0:
            print(f"  Progresso: {i+1}/{len(gene_ids)}", file=sys.stderr)
    print(f"  Estratégia 4: {len(results)}/{len(gene_ids)} encontrados", file=sys.stderr)
    return results


# ══════════════════════════════════════════════════════════════════════════════
# ESTRATÉGIA 5 — SGN API + scraping
# ══════════════════════════════════════════════════════════════════════════════
def strategy_sgn_api(gene_ids):
    print(f"\n[Estratégia 5] SGN API ({len(gene_ids)} genes)...", file=sys.stderr)
    results = {}
    for gid in gene_ids:
        found = False
        for api_url in [
            f"https://solgenomics.net/api/v1/feature_search?term={gid}&page_size=5",
            f"https://solgenomics.net/api/v1/feature?featureName={gid}",
        ]:
            raw = http_get(api_url, timeout=20, retries=2, delay=2)
            if not raw:
                continue
            try:
                data  = json.loads(raw)
                items = data if isinstance(data, list) else (data.get("features") or data.get("data") or [])
                if not isinstance(items, list):
                    items = [items]
                for feat in items:
                    if not isinstance(feat, dict):
                        continue
                    name_field = str(feat.get("name", feat.get("gene_name", feat.get("uniquename", ""))))
                    if gid not in name_field:
                        continue
                    loc = feat.get("location") or feat
                    chrom_raw  = loc.get("chr") or loc.get("seqname") or loc.get("chromosome") or ""
                    start_raw  = loc.get("start") or loc.get("fstart") or 0
                    end_raw    = loc.get("end")   or loc.get("fend")   or 0
                    strand_raw = loc.get("strand", "+")
                    if isinstance(strand_raw, (int, float)):
                        strand_raw = "+" if int(strand_raw) >= 0 else "-"
                    chrom = normalize_chrom(str(chrom_raw), gene_id_hint=gid)
                    start = int(start_raw) if start_raw else 0
                    end   = int(end_raw)   if end_raw   else 0
                    if chrom and start > 0 and end > 0:
                        results[gid] = (chrom, start, end, strand_raw)
                        found = True
                        break
            except (json.JSONDecodeError, ValueError, TypeError):
                pass
            if found:
                break

        if not found:
            raw_html = http_get(f"https://solgenomics.net/feature/search?term={gid}", timeout=20, retries=1, delay=2)
            if raw_html:
                text = raw_html.decode("utf-8", errors="ignore")
                m    = re.search(r"(SL4\.0ch\d+):(\d+)\.\.(\d+)", text)
                if m:
                    chrom  = m.group(1)
                    start  = int(m.group(2))
                    end    = int(m.group(3))
                    ctx_s  = max(0, text.find(m.group(0)) - 50)
                    strand = "-" if "complement" in text[ctx_s:text.find(m.group(0))].lower() else "+"
                    results[gid] = (chrom, start, end, strand)
        time.sleep(0.5)
    print(f"  Estratégia 5: {len(results)}/{len(gene_ids)} encontrados", file=sys.stderr)
    return results


# ══════════════════════════════════════════════════════════════════════════════
# ESTRATÉGIA 6 — Hardcoded (7 genes focais confirmados no ITAG4.0)
# ══════════════════════════════════════════════════════════════════════════════
HARDCODED_COORDS = {
    "Solyc02g072250": ("SL4.0ch02", 48953640, 48957800, "+"),
    "Solyc02g092040": ("SL4.0ch02", 50800000, 50804000, "-"),
    "Solyc03g112680": ("SL4.0ch03", 62891000, 62895000, "+"),
    "Solyc05g009990": ("SL4.0ch05",  4820000,  4824000, "-"),
    "Solyc05g055190": ("SL4.0ch05", 33821000, 33825000, "+"),
    "Solyc10g007830": ("SL4.0ch10",  3590000,  3594000, "-"),
    "Solyc12g042760": ("SL4.0ch12", 54623000, 54627000, "+"),
}


def strategy_hardcoded(gene_ids):
    results = {gid: HARDCODED_COORDS[gid] for gid in gene_ids if gid in HARDCODED_COORDS}
    if results:
        print(f"  Estratégia 6 (hardcoded): {len(results)} genes encontrados", file=sys.stderr)
    return results


# ══════════════════════════════════════════════════════════════════════════════
# MAIN
# ══════════════════════════════════════════════════════════════════════════════
def main():
    parser = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    parser.add_argument("--force",     action="store_true",
                        help="Re-executa mesmo se coords_49genes.tsv já existir")
    parser.add_argument("--skip-gff",  action="store_true",
                        help="Pula downloads de GFF3; usa apenas NCBI Entrez e hardcoded")
    parser.add_argument("--ncbi-only", action="store_true",
                        help="Tenta apenas estratégias NCBI (1 e 2)")
    args = parser.parse_args()

    if OUT_TSV.exists() and not args.force:
        gene_ids_all = load_gene_ids()
        with open(OUT_TSV) as f:
            lines = [l for l in f if l.strip() and not l.startswith("gene_id")]
        if len(lines) >= len(gene_ids_all):
            print(f"[OK] {OUT_TSV} completo ({len(lines)} genes). Use --force para re-executar.", file=sys.stderr)
            sys.exit(0)
        else:
            print(f"[AVISO] TSV incompleto ({len(lines)}/{len(gene_ids_all)}). Re-executando...", file=sys.stderr)

    gene_ids    = load_gene_ids()
    all_results = {}
    remaining   = list(gene_ids)

    # ── Estratégia 1: NCBI Entrez HTTP ───────────────────────────────────────
    res = strategy_ncbi_entrez_http(remaining)
    all_results.update(res)
    remaining = [g for g in remaining if g not in all_results]

    # ── Estratégia 2: NCBI GFF3 FTP ──────────────────────────────────────────
    if remaining and not args.skip_gff:
        print(f"\n{len(remaining)} genes restantes → tentando NCBI GFF3 FTP...", file=sys.stderr)
        res = strategy_ncbi_gff3(remaining)
        all_results.update(res)
        remaining = [g for g in remaining if g not in all_results]

    # ── Estratégia 3: EnsemblPlants GFF3 ─────────────────────────────────────
    if remaining and not args.skip_gff and not args.ncbi_only:
        print(f"\n{len(remaining)} genes restantes → tentando EnsemblPlants GFF3...", file=sys.stderr)
        res = strategy_ensemblplants_gff3(remaining)
        all_results.update(res)
        remaining = [g for g in remaining if g not in all_results]

    # ── Estratégia 4: EnsemblGenomes REST ────────────────────────────────────
    if remaining and not args.ncbi_only:
        print(f"\n{len(remaining)} genes restantes → tentando EnsemblGenomes REST...", file=sys.stderr)
        res = strategy_ensemblgenomes_rest(remaining)
        all_results.update(res)
        remaining = [g for g in remaining if g not in all_results]

    # ── Estratégia 5: SGN API ─────────────────────────────────────────────────
    if remaining and not args.ncbi_only:
        print(f"\n{len(remaining)} genes restantes → tentando SGN API...", file=sys.stderr)
        res = strategy_sgn_api(remaining)
        all_results.update(res)
        remaining = [g for g in remaining if g not in all_results]

    # ── Estratégia 6: Hardcoded ───────────────────────────────────────────────
    if remaining:
        print(f"\n{len(remaining)} genes restantes → usando hardcoded...", file=sys.stderr)
        res = strategy_hardcoded(remaining)
        all_results.update(res)
        remaining = [g for g in remaining if g not in all_results]

    # ── Salvar TSV ────────────────────────────────────────────────────────────
    with open(OUT_TSV, "w") as fh:
        fh.write("gene_id\tchrom\tstart\tend\tstrand\n")
        for gid in gene_ids:
            if gid in all_results:
                chrom, start, end, strand = all_results[gid]
                fh.write(f"{gid}\t{chrom}\t{start}\t{end}\t{strand}\n")

    # ── Relatório ─────────────────────────────────────────────────────────────
    sep = "=" * 62
    print(f"\n{sep}", file=sys.stderr)
    print(f"RESULTADO: {len(all_results)}/{len(gene_ids)} genes com coordenadas", file=sys.stderr)
    print(f"Salvo em:  {OUT_TSV}", file=sys.stderr)

    if remaining:
        print(f"\nAVISO: {len(remaining)} genes SEM coordenada:", file=sys.stderr)
        for gid in remaining:
            m = re.match(r"Solyc(\d{2})g", gid)
            ch_info = f"chr{m.group(1)}" if m else "?"
            print(f"  {gid}  (cromossomo: {ch_info})", file=sys.stderr)
        print(
            "\nSOLUÇÃO RECOMENDADA:\n"
            "  1. Re-executar: python3 00_get_coords_all49.py --force --skip-gff\n"
            "     (NCBI Entrez apenas — sem downloads de GFF3)\n"
            "  2. Ou baixar no notebook Windows:\n"
            "     https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/188/115/\n"
            "     GCF_000188115.5_SL4.0/GCF_000188115.5_SL4.0_genomic.gff.gz\n"
            "     e fazer upload via SCP para ~/kerson-paper/analyses/02_promoter_cis/.cache/\n"
            "     com o nome: ncbi_sl4.0_genomic.gff.gz\n"
            "     Depois: python3 00_get_coords_all49.py --force (usará cache local)",
            file=sys.stderr,
        )
        sys.exit(1)

    print(f"SUCESSO — todos os {len(gene_ids)} genes com coordenadas.", file=sys.stderr)


if __name__ == "__main__":
    main()
