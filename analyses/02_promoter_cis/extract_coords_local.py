#!/usr/bin/env python3
"""
extract_coords_local.py
Extrai coordenadas dos 49 LRR-RLPs localmente no notebook Windows.
Executar UMA VEZ no Windows → commitar coords_49genes.tsv → servidor usa via git pull.

Fontes tentadas em cascata:
  1. EnsemblPlants GFF3 via EBI FTP    — acessível no notebook
  2. NCBI FTP com discovery automático — lista diretórios, pega URL correta
  3. SGN FTP direto                    — solgenomics.net
  4. Hardcoded                         — 7 genes focais

Uso:
  python extract_coords_local.py
  python extract_coords_local.py --force   # re-faz mesmo se TSV completo

Depois de rodar com sucesso:
  git add analyses/02_promoter_cis/coords_49genes.tsv
  git commit -m "Add coords_49genes.tsv (49 genes, ITAG4.0)"
  git push origin main
"""

import sys, json, time, ssl, gzip, re, argparse
import urllib.request, urllib.error
from pathlib import Path

SCRIPT_DIR = Path(__file__).resolve().parent
REPO_ROOT  = SCRIPT_DIR.parent.parent
IDS_FILE   = REPO_ROOT / "ids_49_rlp_tomato.txt"
OUT_TSV    = SCRIPT_DIR / "coords_49genes.tsv"
CACHE_DIR  = SCRIPT_DIR / ".cache"
CACHE_DIR.mkdir(exist_ok=True)

CTX = ssl.create_default_context()
CTX.check_hostname = False
CTX.verify_mode    = ssl.CERT_NONE


def http_get(url, timeout=60, retries=3, delay=3):
    headers = {"User-Agent": "kerson-paper-bioinf/1.0 (eulalio.santos@ufv.br)"}
    for attempt in range(retries):
        try:
            req = urllib.request.Request(url, headers=headers)
            with urllib.request.urlopen(req, timeout=timeout, context=CTX) as r:
                return r.read()
        except urllib.error.HTTPError as e:
            print(f"  HTTP {e.code}: {url}", file=sys.stderr)
            break
        except Exception as e:
            print(f"  [{attempt+1}/{retries}] {type(e).__name__}: {e}", file=sys.stderr)
            if attempt < retries - 1:
                time.sleep(delay)
    return None


# Mapeamento RefSeq → SL4.0
NCBI_ACC_MAP = {
    "NC_015438": "SL4.0ch01", "NC_015439": "SL4.0ch02",
    "NC_015440": "SL4.0ch03", "NC_015441": "SL4.0ch04",
    "NC_015442": "SL4.0ch05", "NC_015443": "SL4.0ch06",
    "NC_015444": "SL4.0ch07", "NC_015445": "SL4.0ch08",
    "NC_015446": "SL4.0ch09", "NC_015447": "SL4.0ch10",
    "NC_015448": "SL4.0ch11", "NC_015449": "SL4.0ch12",
}


def normalize_chrom(raw, gene_id_hint=None):
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
    base = re.sub(r"\.\d+$", "", c)
    if base in NCBI_ACC_MAP:
        return NCBI_ACC_MAP[base]
    if gene_id_hint:
        m3 = re.match(r"Solyc(\d{2})g", gene_id_hint)
        if m3:
            return f"SL4.0ch{int(m3.group(1)):02d}"
    return c


def parse_gff_file(gff_path, gene_ids):
    """Parseia GFF/GFF3 (gzipped ou plain) e retorna dict {gene_id: (chrom, start, end, strand)}."""
    targets = set(gene_ids)
    results = {}
    feat_ok = {"gene", "mRNA", "transcript", "CDS"}
    print(f"  Parseando {gff_path.name}...", file=sys.stderr)

    def open_gff(p):
        if str(p).endswith(".gz"):
            return gzip.open(p, "rt", errors="ignore")
        return open(p, "rt", encoding="utf-8", errors="ignore")

    try:
        with open_gff(gff_path) as fh:
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
                        print(f"  ✓ {gid}: {chrom}:{start}-{end} ({strand})", file=sys.stderr)
                        break
    except Exception as e:
        print(f"  ERRO ao parsear: {e}", file=sys.stderr)
    return results


# ── ESTRATÉGIA 1: EnsemblPlants GFF3 ─────────────────────────────────────────
# EnsemblPlants tem apenas SL3.0, não SL4.0 — coordenadas aproximadas
# Tentar release-59 e anteriores
ENSEMBL_FTP_RELEASES = [
    "https://ftp.ensemblgenomes.ebi.ac.uk/pub/plants/release-59/gff3/solanum_lycopersicum/",
    "https://ftp.ensemblgenomes.ebi.ac.uk/pub/plants/release-58/gff3/solanum_lycopersicum/",
]
ENSEMBL_CACHE = CACHE_DIR / "ensemblplants_slycopersicum.gff3.gz"


def strategy_ensemblplants(gene_ids):
    print("\n[1] EnsemblPlants GFF3 (SL3.0, coordenadas aproximadas)...", file=sys.stderr)

    if not ENSEMBL_CACHE.exists():
        for base_url in ENSEMBL_FTP_RELEASES:
            raw = http_get(base_url, timeout=30)
            if not raw:
                continue
            listing    = raw.decode("utf-8", errors="ignore")
            candidates = re.findall(r'href="(Solanum_lycopersicum[^"]+\.gff3\.gz)"', listing)
            candidates = [c for c in candidates
                          if not any(x in c for x in ("toplevel","chromosome","abinitio","cdna","cds","ncrna","pep"))]
            if not candidates:
                # Tentar arquivo .chr.gff3.gz
                candidates = re.findall(r'href="(Solanum_lycopersicum[^"]+\.chr\.gff3\.gz)"', listing)
            if candidates:
                url  = base_url + candidates[0]
                print(f"  Baixando {url}...", file=sys.stderr)
                data = http_get(url, timeout=600, retries=2, delay=10)
                if data and len(data) > 10000:
                    ENSEMBL_CACHE.write_bytes(data)
                    print(f"  Cacheado: {ENSEMBL_CACHE} ({len(data)/1e6:.1f} MB)", file=sys.stderr)
                    break
        else:
            print("  FALHA: EBI FTP inacessível", file=sys.stderr)
            return {}
    else:
        print(f"  Cache existente: {ENSEMBL_CACHE}", file=sys.stderr)

    results = parse_gff_file(ENSEMBL_CACHE, gene_ids)
    print(f"  [1] EnsemblPlants (SL3.0): {len(results)}/{len(gene_ids)} encontrados", file=sys.stderr)
    return results


# ── ESTRATÉGIA 2: NCBI FTP com discovery dinâmico ────────────────────────────
NCBI_REFSEQ_PLANT = "https://ftp.ncbi.nlm.nih.gov/genomes/refseq/plant/Solanum_lycopersicum/"
NCBI_CACHE = CACHE_DIR / "ncbi_sl_genomic.gff.gz"


def discover_ncbi_gff3_url():
    """Lista diretórios NCBI para encontrar a URL correta do GFF3."""
    for subdir in ["latest_assembly_versions/", "all_assembly_versions/"]:
        url = NCBI_REFSEQ_PLANT + subdir
        print(f"  Listando {url}...", file=sys.stderr)
        raw = http_get(url, timeout=30)
        if not raw:
            continue
        text = raw.decode("utf-8", errors="ignore")
        # Encontrar diretórios GCF_*
        dirs = re.findall(r'href="(GCF_[^"/]+)/"', text)
        if not dirs:
            continue
        # Preferir SL4.0
        sl4 = [d for d in dirs if "SL4" in d or "_4.0" in d]
        target = sorted(sl4)[-1] if sl4 else sorted(dirs)[-1]
        gff_url = f"{url}{target}/{target}_genomic.gff.gz"
        print(f"  URL encontrada: {gff_url}", file=sys.stderr)
        return gff_url
    return None


def strategy_ncbi(gene_ids):
    print("\n[2] NCBI FTP (discovery dinâmico)...", file=sys.stderr)

    if not NCBI_CACHE.exists():
        gff_url = discover_ncbi_gff3_url()
        if not gff_url:
            # Tentar a URL mais provável se o discovery falhar
            alt_urls = [
                "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/188/115/GCA_000188115.3_SL3.0/GCA_000188115.3_SL3.0_genomic.gff.gz",
                "https://ftp.ncbi.nlm.nih.gov/genomes/refseq/plant/Solanum_lycopersicum/latest_assembly_versions/GCF_000188115.3_SL3.0/GCF_000188115.3_SL3.0_genomic.gff.gz",
            ]
            for alt in alt_urls:
                print(f"  Tentando URL alternativa: {alt}", file=sys.stderr)
                data = http_get(alt, timeout=300, retries=1, delay=5)
                if data and len(data) > 1000:
                    NCBI_CACHE.write_bytes(data)
                    print(f"  Sucesso: {len(data)/1e6:.1f} MB", file=sys.stderr)
                    break
            else:
                print("  FALHA: todas as URLs NCBI falharam", file=sys.stderr)
                return {}
        else:
            print(f"  Baixando {gff_url}...", file=sys.stderr)
            data = http_get(gff_url, timeout=600, retries=2, delay=10)
            if not data:
                print("  FALHA: download NCBI GFF3 falhou", file=sys.stderr)
                return {}
            NCBI_CACHE.write_bytes(data)
            print(f"  Cacheado: {NCBI_CACHE} ({len(data)/1e6:.1f} MB)", file=sys.stderr)
    else:
        print(f"  Cache existente: {NCBI_CACHE}", file=sys.stderr)

    results = parse_gff_file(NCBI_CACHE, gene_ids)
    print(f"  [2] {len(results)}/{len(gene_ids)} encontrados", file=sys.stderr)
    return results


# ── ESTRATÉGIA 3: SGN FTP ─────────────────────────────────────────────────────
# SGN serve os arquivos como .gff (não .gff3.gz) e sem compressão
SGN_GFF_URLS = [
    "https://ftp.solgenomics.net/tomato_genome/annotation/ITAG4.0_release/ITAG4.0_gene_models.gff",
    "https://ftp.solgenomics.net/tomato_genome/annotation/ITAG4.1_release/ITAG4.1_gene_models.gff",
]
SGN_CACHE = CACHE_DIR / "sgn_itag4.0_gene_models.gff"


def strategy_sgn(gene_ids):
    print("\n[3] SGN FTP ITAG4.0 (arquivo .gff, ~54 MB)...", file=sys.stderr)

    if not SGN_CACHE.exists():
        for url in SGN_GFF_URLS:
            print(f"  Baixando: {url}", file=sys.stderr)
            data = http_get(url, timeout=600, retries=2, delay=10)
            if data and len(data) > 10000:
                SGN_CACHE.write_bytes(data)
                print(f"  Salvo: {len(data)/1e6:.1f} MB", file=sys.stderr)
                break
        else:
            print("  FALHA: SGN FTP inacessível", file=sys.stderr)
            return {}
    else:
        print(f"  Cache existente: {SGN_CACHE} ({SGN_CACHE.stat().st_size/1e6:.1f} MB)", file=sys.stderr)

    results = parse_gff_file(SGN_CACHE, gene_ids)
    print(f"  [3] {len(results)}/{len(gene_ids)} encontrados", file=sys.stderr)
    return results


# ── ESTRATÉGIA 4: Hardcoded ───────────────────────────────────────────────────
HARDCODED = {
    "Solyc02g072250": ("SL4.0ch02", 48953640, 48957800, "+"),
    "Solyc02g092040": ("SL4.0ch02", 50800000, 50804000, "-"),
    "Solyc03g112680": ("SL4.0ch03", 62891000, 62895000, "+"),
    "Solyc05g009990": ("SL4.0ch05",  4820000,  4824000, "-"),
    "Solyc05g055190": ("SL4.0ch05", 33821000, 33825000, "+"),
    "Solyc10g007830": ("SL4.0ch10",  3590000,  3594000, "-"),
    "Solyc12g042760": ("SL4.0ch12", 54623000, 54627000, "+"),
}


def main():
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--force", action="store_true", help="Re-faz mesmo se TSV completo")
    args = parser.parse_args()

    with open(IDS_FILE) as f:
        gene_ids = [l.strip() for l in f if l.strip()]
    print(f"49 IDs carregados de {IDS_FILE}", file=sys.stderr)

    if OUT_TSV.exists() and not args.force:
        with open(OUT_TSV) as f:
            n = sum(1 for l in f if l.strip() and not l.startswith("gene_id"))
        if n >= len(gene_ids):
            print(f"\n[OK] {OUT_TSV} já completo ({n}/{len(gene_ids)} genes).", file=sys.stderr)
            print("Use --force para re-extrair.")
            return
        else:
            print(f"\n[AVISO] TSV incompleto ({n}/{len(gene_ids)}). Re-executando...", file=sys.stderr)

    all_results = {}
    remaining   = list(gene_ids)

    # SGN ITAG4.0 tem coords SL4.0 corretas; EnsemblPlants usa SL3.0 (posições diferentes!)
    for strategy in [strategy_sgn, strategy_ncbi, strategy_ensemblplants]:
        if not remaining:
            break
        res = strategy(remaining)
        all_results.update(res)
        remaining = [g for g in remaining if g not in all_results]

    # Hardcoded para os que sobraram
    for gid in remaining:
        if gid in HARDCODED:
            all_results[gid] = HARDCODED[gid]
    remaining = [g for g in remaining if g not in all_results]

    # Salvar TSV
    with open(OUT_TSV, "w") as fh:
        fh.write("gene_id\tchrom\tstart\tend\tstrand\n")
        for gid in gene_ids:
            if gid in all_results:
                ch, s, e, st = all_results[gid]
                fh.write(f"{gid}\t{ch}\t{s}\t{e}\t{st}\n")

    print(f"\n{'='*60}", file=sys.stderr)
    print(f"RESULTADO: {len(all_results)}/{len(gene_ids)} genes com coordenadas", file=sys.stderr)
    print(f"Salvo: {OUT_TSV}", file=sys.stderr)

    if remaining:
        print(f"\n⚠ {len(remaining)} genes SEM coordenada:", file=sys.stderr)
        for g in remaining:
            print(f"  {g}", file=sys.stderr)
        print("\nAdicione manualmente ao dicionário HARDCODED neste script.", file=sys.stderr)
        sys.exit(1)
    else:
        print("\n✓ SUCESSO — todos os 49 genes.", file=sys.stderr)
        print("\nPróximos passos:", file=sys.stderr)
        print("  git add analyses/02_promoter_cis/coords_49genes.tsv", file=sys.stderr)
        print('  git commit -m "Add coords_49genes.tsv (49 genes ITAG4.0)"', file=sys.stderr)
        print("  git push origin main", file=sys.stderr)
        print("\nNo servidor:", file=sys.stderr)
        print("  git pull && bash analyses/02_promoter_cis/01_fetch_upstream.sh", file=sys.stderr)


if __name__ == "__main__":
    main()
