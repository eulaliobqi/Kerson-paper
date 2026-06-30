#!/usr/bin/env python3
"""
00_get_coords_all49.py
Obtém coordenadas genômicas ITAG4.0 dos 49 LRR-RLPs de Solanum lycopersicum.

Estratégias em cascata (para quando a primeira falha):
  1. EnsemblPlants GFF3 via EBI FTP  — download único, resolve todos de uma vez
  2. EnsemblGenomes REST API          — por gene, com bypass SSL e retry
  3. SGN API / scraping mínimo        — solgenomics.net, por gene
  4. NCBI Entrez via Biopython        — busca por Gene Name + organismo
  5. Hardcoded                        — 7 posições confirmadas do ITAG4.0

Saída: coords_49genes.tsv
  Colunas: gene_id  chrom  start  end  strand
  chrom no formato SL4.0chNN (compatível com S_lycopersicum_chromosomes.4.00.fa)

Uso:
  python3 00_get_coords_all49.py            # pula se TSV já existe
  python3 00_get_coords_all49.py --force    # re-executa sempre
  python3 00_get_coords_all49.py --skip-gff # pula download do GFF3 (apenas REST/SGN/NCBI)
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

# ── SSL context sem verificação de certificado ────────────────────────────────
CTX = ssl.create_default_context()
CTX.check_hostname = False
CTX.verify_mode    = ssl.CERT_NONE


def http_get(url, timeout=30, retries=3, delay=2):
    """GET com retry, bypass SSL e User-Agent definido. Retorna bytes ou None."""
    headers = {"User-Agent": "kerson-paper-bioinf/1.0 (eulalio.santos@ufv.br)"}
    for attempt in range(retries):
        try:
            req = urllib.request.Request(url, headers=headers)
            with urllib.request.urlopen(req, timeout=timeout, context=CTX) as r:
                return r.read()
        except urllib.error.HTTPError as e:
            print(f"    HTTP {e.code} — {url}", file=sys.stderr)
            break  # HTTP 4xx/5xx: não retenta
        except Exception as e:
            print(f"    [tentativa {attempt + 1}/{retries}] {type(e).__name__}: {e}", file=sys.stderr)
            if attempt < retries - 1:
                time.sleep(delay)
    return None


# ── Normalização de nomes de cromossomo → SL4.0chNN ──────────────────────────
def normalize_chrom(raw):
    """Converte qualquer nome de cromossomo para SL4.0chNN."""
    c = str(raw).strip()
    if re.match(r"^SL4\.0ch\d+$", c):
        return c
    # Ensembl Plants usa apenas o número: "1", "2", ...
    if re.match(r"^\d{1,2}$", c):
        return f"SL4.0ch{int(c):02d}"
    # "ch01", "Ch1", "Chr01", "chr1" etc.
    m = re.match(r"^[Cc][Hh][Rr]?0*(\d+)$", c)
    if m:
        return f"SL4.0ch{int(m.group(1)):02d}"
    # Formas como "SL2.40ch01" (versão antiga) — mapear para SL4.0
    m2 = re.match(r"^SL\d+\.\d+ch0*(\d+)$", c)
    if m2:
        return f"SL4.0ch{int(m2.group(1)):02d}"
    return c  # mantém como está; será validado pelo usuário


# ── Carga dos 49 IDs ─────────────────────────────────────────────────────────
def load_gene_ids():
    with open(IDS_FILE) as f:
        ids = [l.strip() for l in f if l.strip()]
    print(f"IDs carregados de {IDS_FILE}: {len(ids)}", file=sys.stderr)
    return ids


# ══════════════════════════════════════════════════════════════════════════════
# ESTRATÉGIA 1 — EnsemblPlants GFF3 (EBI FTP)
# ══════════════════════════════════════════════════════════════════════════════
ENSEMBL_FTP_BASE = "https://ftp.ensemblgenomes.ebi.ac.uk/pub/plants/current/gff3/solanum_lycopersicum/"
ENSEMBL_GFF3_CACHE = CACHE_DIR / "ensemblplants_slycopersicum.gff3.gz"


def strategy_ensemblplants_gff3(gene_ids):
    """Baixa GFF3 completo do EnsemblPlants e extrai coordenadas em uma passagem."""
    print("\n[Estratégia 1] EnsemblPlants GFF3 (EBI FTP)...", file=sys.stderr)

    if not ENSEMBL_GFF3_CACHE.exists():
        print("  Obtendo índice do diretório EBI FTP...", file=sys.stderr)
        listing_raw = http_get(ENSEMBL_FTP_BASE, timeout=25)
        if not listing_raw:
            print("  FALHA: EBI FTP inacessível", file=sys.stderr)
            return {}

        listing = listing_raw.decode("utf-8", errors="ignore")

        # Procurar o arquivo gene-level GFF3 (sem toplevel, chromosome, abinitio)
        candidates = re.findall(
            r'href="(Solanum_lycopersicum[^"]+\.gff3\.gz)"', listing
        )
        candidates = [
            c for c in candidates
            if not any(x in c for x in ("toplevel", "chromosome", "abinitio", "cdna", "cds", "ncrna", "pep"))
        ]

        if not candidates:
            print("  FALHA: nenhum arquivo GFF3 encontrado no diretório", file=sys.stderr)
            return {}

        gff3_file = candidates[0]
        gff3_url  = ENSEMBL_FTP_BASE + gff3_file
        print(f"  Baixando {gff3_url} (pode demorar ~2 min)...", file=sys.stderr)

        data = http_get(gff3_url, timeout=300, retries=2, delay=10)
        if not data:
            print("  FALHA: download do GFF3 falhou", file=sys.stderr)
            return {}

        ENSEMBL_GFF3_CACHE.write_bytes(data)
        print(f"  GFF3 cacheado: {ENSEMBL_GFF3_CACHE}", file=sys.stderr)
    else:
        print(f"  Usando cache: {ENSEMBL_GFF3_CACHE}", file=sys.stderr)

    # Parsear GFF3: busca genes pelos IDs alvo
    targets  = set(gene_ids)
    results  = {}
    feat_types = {"gene", "mRNA", "transcript"}

    try:
        with gzip.open(ENSEMBL_GFF3_CACHE, "rt", errors="ignore") as fh:
            for line in fh:
                if not targets:
                    break
                if line.startswith("#"):
                    continue
                cols = line.rstrip("\n").split("\t")
                if len(cols) < 9 or cols[2] not in feat_types:
                    continue
                attrs = cols[8]
                for gid in list(targets):
                    # Solyc IDs podem aparecer em maiúsculas no GFF3 do EnsemblPlants
                    if gid in attrs or gid.upper() in attrs:
                        chrom  = normalize_chrom(cols[0])
                        start  = int(cols[3])
                        end    = int(cols[4])
                        strand = cols[6] if cols[6] in ("+", "-") else "+"
                        results[gid] = (chrom, start, end, strand)
                        targets.discard(gid)
                        print(f"    {gid}: {chrom}:{start}-{end} ({strand})", file=sys.stderr)
                        break
    except Exception as e:
        print(f"  ERRO ao parsear GFF3: {e}", file=sys.stderr)

    print(f"  Estratégia 1: {len(results)}/{len(gene_ids)} encontrados", file=sys.stderr)
    return results


# ══════════════════════════════════════════════════════════════════════════════
# ESTRATÉGIA 2 — EnsemblGenomes REST API (por gene)
# ══════════════════════════════════════════════════════════════════════════════
def strategy_ensemblgenomes_rest(gene_ids):
    """EnsemblGenomes REST API com bypass SSL, retry e rate-limit de 1 req/s."""
    print(f"\n[Estratégia 2] EnsemblGenomes REST API ({len(gene_ids)} genes)...", file=sys.stderr)
    results = {}

    for i, gid in enumerate(gene_ids):
        # Tentar: ID original, ID em maiúsculas, sem versão
        id_variants = [gid, gid.upper()]
        for variant in id_variants:
            url = (
                f"https://rest.ensemblgenomes.org/lookup/id/{variant}"
                f"?content-type=application/json&expand=0"
            )
            raw = http_get(url, timeout=25, retries=2, delay=1)
            if not raw:
                continue
            try:
                d = json.loads(raw)
                if "seq_region_name" in d and "start" in d:
                    chrom  = normalize_chrom(d["seq_region_name"])
                    start  = int(d["start"])
                    end    = int(d["end"])
                    strand = "+" if d.get("strand", 1) == 1 else "-"
                    results[gid] = (chrom, start, end, strand)
                    print(f"    {gid}: {chrom}:{start}-{end} ({strand})", file=sys.stderr)
                    break
                if "error" in d:
                    print(f"    REST erro {gid}: {d['error']}", file=sys.stderr)
            except (json.JSONDecodeError, ValueError, KeyError):
                pass

        time.sleep(1.0)  # respeitar rate-limit EnsemblGenomes: ≤1 req/s sem token

        if (i + 1) % 10 == 0:
            print(f"  Progresso: {i + 1}/{len(gene_ids)}", file=sys.stderr)

    print(f"  Estratégia 2: {len(results)}/{len(gene_ids)} encontrados", file=sys.stderr)
    return results


# ══════════════════════════════════════════════════════════════════════════════
# ESTRATÉGIA 3 — SGN API + scraping mínimo
# ══════════════════════════════════════════════════════════════════════════════
def strategy_sgn_api(gene_ids):
    """SGN feature search API e scraping mínimo de página HTML."""
    print(f"\n[Estratégia 3] SGN API ({len(gene_ids)} genes)...", file=sys.stderr)
    results = {}

    for gid in gene_ids:
        found = False

        # 3a. API JSON
        for api_url in [
            f"https://solgenomics.net/api/v1/feature_search?term={gid}&page_size=5",
            f"https://solgenomics.net/api/v1/feature?featureName={gid}",
        ]:
            raw = http_get(api_url, timeout=20, retries=2, delay=2)
            if not raw:
                continue
            try:
                data = json.loads(raw)
                # SGN pode retornar lista direta ou objeto com "features"/"data"
                items = data if isinstance(data, list) else (
                    data.get("features") or data.get("data") or []
                )
                if not isinstance(items, list):
                    items = [items]

                for feat in items:
                    if not isinstance(feat, dict):
                        continue
                    name_field = str(feat.get("name", feat.get("gene_name", feat.get("uniquename", ""))))
                    if gid not in name_field:
                        continue
                    # Tentar extrair localização de diferentes layouts JSON do SGN
                    loc = feat.get("location") or feat
                    chrom_raw = (
                        loc.get("chr") or loc.get("seqname") or
                        loc.get("chromosome") or loc.get("srcfeature_id") or ""
                    )
                    start_raw = loc.get("start") or loc.get("fstart") or loc.get("pos_start") or 0
                    end_raw   = loc.get("end")   or loc.get("fend")   or loc.get("pos_end")   or 0
                    strand_raw = loc.get("strand", loc.get("fstrand", "+"))
                    if isinstance(strand_raw, (int, float)):
                        strand_raw = "+" if int(strand_raw) >= 0 else "-"
                    chrom = normalize_chrom(str(chrom_raw))
                    start = int(start_raw) if start_raw else 0
                    end   = int(end_raw)   if end_raw   else 0
                    if chrom and start > 0 and end > 0:
                        results[gid] = (chrom, start, end, strand_raw)
                        print(f"    {gid} (API): {chrom}:{start}-{end} ({strand_raw})", file=sys.stderr)
                        found = True
                        break
            except (json.JSONDecodeError, ValueError, TypeError):
                pass
            if found:
                break

        # 3b. Scraping HTML mínimo da página de busca do SGN
        if not found:
            html_url = f"https://solgenomics.net/feature/search?term={gid}"
            raw_html = http_get(html_url, timeout=20, retries=1, delay=2)
            if raw_html:
                text = raw_html.decode("utf-8", errors="ignore")
                # Padrão SL4.0chXX:nnnnn..nnnnn no HTML
                m = re.search(r"(SL4\.0ch\d+):(\d+)\.\.(\d+)", text)
                if m:
                    chrom  = m.group(1)
                    start  = int(m.group(2))
                    end    = int(m.group(3))
                    # "complement" antes da posição indica strand -
                    ctx_start = max(0, text.find(m.group(0)) - 50)
                    strand = "-" if "complement" in text[ctx_start:text.find(m.group(0))].lower() else "+"
                    results[gid] = (chrom, start, end, strand)
                    print(f"    {gid} (scraping): {chrom}:{start}-{end} ({strand})", file=sys.stderr)
                    found = True

        time.sleep(0.5)

    print(f"  Estratégia 3: {len(results)}/{len(gene_ids)} encontrados", file=sys.stderr)
    return results


# ══════════════════════════════════════════════════════════════════════════════
# ESTRATÉGIA 4 — NCBI Entrez via Biopython
# ══════════════════════════════════════════════════════════════════════════════
def strategy_ncbi_entrez(gene_ids):
    """Busca via NCBI Gene por nome do gene + Solanum lycopersicum."""
    print(f"\n[Estratégia 4] NCBI Entrez ({len(gene_ids)} genes)...", file=sys.stderr)

    try:
        from Bio import Entrez
    except ImportError:
        print("  Biopython não disponível — pulando NCBI Entrez", file=sys.stderr)
        return {}

    Entrez.email = "eulalio.santos@ufv.br"
    Entrez.tool  = "kerson-paper-rlp-analysis"
    results = {}

    for gid in gene_ids:
        try:
            # Busca pelo símbolo do gene
            for query in [
                f'"{gid}"[Gene Symbol] AND "Solanum lycopersicum"[Organism]',
                f'"{gid}"[Gene Name] AND "Solanum lycopersicum"[Organism]',
                f"{gid} AND txid4081[Organism]",
            ]:
                handle = Entrez.esearch(db="gene", term=query, retmax=3)
                record = Entrez.read(handle)
                handle.close()
                if record["IdList"]:
                    break

            if not record["IdList"]:
                time.sleep(0.35)
                continue

            ncbi_gene_id = record["IdList"][0]
            handle = Entrez.efetch(
                db="gene", id=ncbi_gene_id,
                rettype="gene_table", retmode="text"
            )
            gene_table = handle.read()
            handle.close()

            # gene_table é texto tabulado; parsear pela linha de coordenadas genômicas
            for line in str(gene_table).split("\n"):
                if not line.strip() or line.startswith("chr") or "start" in line.lower():
                    continue
                cols = line.strip().split("\t")
                if len(cols) >= 3:
                    try:
                        chrom  = normalize_chrom(cols[0])
                        start  = int(cols[1])
                        end    = int(cols[2])
                        strand = cols[3] if len(cols) > 3 and cols[3] in ("+", "-") else "+"
                        if start > 0 and end > 0:
                            results[gid] = (chrom, start, end, strand)
                            print(f"    {gid}: {chrom}:{start}-{end} ({strand})", file=sys.stderr)
                            break
                    except (ValueError, IndexError):
                        pass

        except Exception as e:
            print(f"    NCBI falhou {gid}: {e}", file=sys.stderr)

        time.sleep(0.35)  # NCBI rate-limit: 3 req/s sem API key

    print(f"  Estratégia 4: {len(results)}/{len(gene_ids)} encontrados", file=sys.stderr)
    return results


# ══════════════════════════════════════════════════════════════════════════════
# ESTRATÉGIA 5 — Hardcoded (posições confirmadas do ITAG4.0)
# ══════════════════════════════════════════════════════════════════════════════
# Fonte das 7 posições: scripts legados do pipeline (coordenadas ITAG4.0 confirmadas).
# Para adicionar os 42 restantes: consultar https://solgenomics.net/feature/search?term=GENE_ID
# ou https://ensemblgenomes.org/id/GENE_ID e preencher abaixo.
HARDCODED_COORDS = {
    "Solyc02g072250": ("SL4.0ch02", 48953640, 48957800, "+"),
    "Solyc02g092040": ("SL4.0ch02", 50800000, 50804000, "-"),
    "Solyc03g112680": ("SL4.0ch03", 62891000, 62895000, "+"),
    "Solyc05g009990": ("SL4.0ch05",  4820000,  4824000, "-"),
    "Solyc05g055190": ("SL4.0ch05", 33821000, 33825000, "+"),
    "Solyc10g007830": ("SL4.0ch10",  3590000,  3594000, "-"),
    "Solyc12g042760": ("SL4.0ch12", 54623000, 54627000, "+"),
    # Adicionar aqui os 42 restantes se as APIs acima falharem:
    # "SolycXXgYYYYYY": ("SL4.0chXX", start, end, "+/-"),
}


def strategy_hardcoded(gene_ids):
    results = {}
    for gid in gene_ids:
        if gid in HARDCODED_COORDS:
            results[gid] = HARDCODED_COORDS[gid]
    if results:
        print(f"  Estratégia 5 (hardcoded): {len(results)} genes encontrados", file=sys.stderr)
    return results


# ══════════════════════════════════════════════════════════════════════════════
# MAIN
# ══════════════════════════════════════════════════════════════════════════════
def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--force",    action="store_true", help="Re-executa mesmo se coords_49genes.tsv já existir")
    parser.add_argument("--skip-gff", action="store_true", help="Pula download do GFF3 do EnsemblPlants (mais rápido, apenas REST/SGN/NCBI)")
    args = parser.parse_args()

    if OUT_TSV.exists() and not args.force:
        print(f"[OK] {OUT_TSV} já existe. Use --force para re-executar.", file=sys.stderr)
        # Verificar quantas linhas tem
        with open(OUT_TSV) as f:
            lines = [l for l in f if l.strip() and not l.startswith("gene_id")]
        print(f"    Contém {len(lines)} genes com coordenadas.", file=sys.stderr)
        sys.exit(0)

    gene_ids   = load_gene_ids()
    all_results = {}
    remaining  = list(gene_ids)

    # ── Estratégia 1: EnsemblPlants GFF3 (download único, mais confiável) ───
    if not args.skip_gff:
        res1 = strategy_ensemblplants_gff3(remaining)
        all_results.update(res1)
        remaining = [g for g in remaining if g not in all_results]
    else:
        print("\n[--skip-gff] Pulando download de GFF3", file=sys.stderr)

    # ── Estratégia 2: EnsemblGenomes REST API ───────────────────────────────
    if remaining:
        print(f"\n{len(remaining)} genes sem coordenada → tentando REST EnsemblGenomes...", file=sys.stderr)
        res2 = strategy_ensemblgenomes_rest(remaining)
        all_results.update(res2)
        remaining = [g for g in remaining if g not in all_results]

    # ── Estratégia 3: SGN API ────────────────────────────────────────────────
    if remaining:
        print(f"\n{len(remaining)} genes sem coordenada → tentando SGN API...", file=sys.stderr)
        res3 = strategy_sgn_api(remaining)
        all_results.update(res3)
        remaining = [g for g in remaining if g not in all_results]

    # ── Estratégia 4: NCBI Entrez ────────────────────────────────────────────
    if remaining:
        print(f"\n{len(remaining)} genes sem coordenada → tentando NCBI Entrez...", file=sys.stderr)
        res4 = strategy_ncbi_entrez(remaining)
        all_results.update(res4)
        remaining = [g for g in remaining if g not in all_results]

    # ── Estratégia 5: Hardcoded ──────────────────────────────────────────────
    if remaining:
        print(f"\n{len(remaining)} genes sem coordenada → usando hardcoded...", file=sys.stderr)
        res5 = strategy_hardcoded(remaining)
        all_results.update(res5)
        remaining = [g for g in remaining if g not in all_results]

    # ── Salvar TSV ───────────────────────────────────────────────────────────
    with open(OUT_TSV, "w") as fh:
        fh.write("gene_id\tchrom\tstart\tend\tstrand\n")
        for gid in gene_ids:
            if gid in all_results:
                chrom, start, end, strand = all_results[gid]
                fh.write(f"{gid}\t{chrom}\t{start}\t{end}\t{strand}\n")

    # ── Relatório final ──────────────────────────────────────────────────────
    sep = "=" * 62
    print(f"\n{sep}", file=sys.stderr)
    print(f"RESULTADO: {len(all_results)}/{len(gene_ids)} genes com coordenadas", file=sys.stderr)
    print(f"Salvo em:  {OUT_TSV}", file=sys.stderr)

    if remaining:
        print(f"\nAVISO: {len(remaining)} genes SEM coordenada:", file=sys.stderr)
        for gid in remaining:
            print(f"  {gid}  ->  https://solgenomics.net/feature/search?term={gid}", file=sys.stderr)
        print(
            "\nPara corrigir:\n"
            "  1. Acesse os links SGN acima e copie as coordenadas\n"
            "  2. Adicione-as ao dicionario HARDCODED_COORDS neste script\n"
            "  3. Re-execute: python3 00_get_coords_all49.py --force --skip-gff\n"
            "  Ou: baixe o GFF3 SGN manualmente e extraia as coords via grep",
            file=sys.stderr,
        )
        sys.exit(1)

    print(f"SUCESSO — todos os {len(gene_ids)} genes com coordenadas.", file=sys.stderr)


if __name__ == "__main__":
    main()
