#!/usr/bin/env python3
"""
fetch_upstream_local.py
Busca regiões 2 kb upstream dos 49 LRR-RLPs via Ensembl Plants REST API.
Executar no Windows (com internet).

Uso:
    python3 fetch_upstream_local.py
    python3 fetch_upstream_local.py --upstream 2000 --out rlp_upstream_2kb.fa

Fonte: https://rest.ensemblgenomes.org (Ensembl Plants r61, SL4.0)
"""

import sys, time, argparse, csv, ssl
from pathlib import Path
from urllib.request import urlopen, Request
from urllib.error import HTTPError, URLError

NCBI_EFETCH = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi"
NCBI_EMAIL  = "eulalio.santos@ufv.br"

# Mapeamento SL4.0chXX → accession NCBI do assembly SL4.0 (GCA_000188115.3)
# Tomato Heinz 1706 ITAG4.0 — cromossomos no GenBank
CHROM_MAP = {
    "SL4.0ch01": "CM001064.3", "SL4.0ch02": "CM001065.3",
    "SL4.0ch03": "CM001066.3", "SL4.0ch04": "CM001067.3",
    "SL4.0ch05": "CM001068.3", "SL4.0ch06": "CM001069.3",
    "SL4.0ch07": "CM001070.3", "SL4.0ch08": "CM001071.3",
    "SL4.0ch09": "CM001072.3", "SL4.0ch10": "CM001073.3",
    "SL4.0ch11": "CM001074.3", "SL4.0ch12": "CM001075.3",
}

# Contexto SSL sem verificação (necessário no Windows para Ensembl Plants;
# NCBI geralmente OK, mas mantemos como fallback)
SSL_CTX = ssl.create_default_context()
SSL_CTX.check_hostname = False
SSL_CTX.verify_mode    = ssl.CERT_NONE


def fetch_seq(chrom_accession, start, end, strand, upstream=2000):
    """
    Busca região upstream via NCBI efetch (nuccore).
    start/end: coordenadas 1-based do gene.
    strand: '+' ou '-'
    Retorna (sequência_str, None) ou (None, mensagem_erro).
    """
    if strand == "+":
        reg_start = max(1, start - upstream)
        reg_end   = start - 1
        strand_id = 1          # NCBI: 1=plus, 2=minus
    else:
        reg_start = end + 1
        reg_end   = end + upstream
        strand_id = 2

    if reg_start >= reg_end:
        return None, f"região upstream inválida ({reg_start}–{reg_end})"

    url = (
        f"{NCBI_EFETCH}?db=nuccore"
        f"&id={chrom_accession}"
        f"&seq_start={reg_start}&seq_stop={reg_end}"
        f"&strand={strand_id}"
        f"&rettype=fasta&retmode=text"
        f"&tool=kerson-paper&email={NCBI_EMAIL}"
    )

    req = Request(url, headers={"User-Agent": "kerson-paper/1.0"})
    for attempt in range(3):
        try:
            with urlopen(req, timeout=30, context=SSL_CTX) as resp:
                raw = resp.read().decode("utf-8").strip()
            # Remover header FASTA e juntar sequência
            lines = raw.splitlines()
            seq = "".join(l.strip() for l in lines if not l.startswith(">")).upper()
            if not seq:
                return None, "resposta FASTA vazia"
            return seq, None
        except HTTPError as e:
            if e.code == 429:
                time.sleep(2 ** (attempt + 1))
                continue
            return None, f"HTTP {e.code}"
        except URLError as e:
            return None, f"URLError: {e.reason}"
    return None, "3 tentativas falharam"


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--upstream", type=int, default=2000)
    parser.add_argument("--coords",   default=None)
    parser.add_argument("--out",      default=None)
    args = parser.parse_args()

    script_dir  = Path(__file__).parent
    coords_file = Path(args.coords) if args.coords else script_dir / "coords_49genes.tsv"
    out_file    = Path(args.out)    if args.out    else script_dir / "rlp_upstream_2kb.fa"

    if not coords_file.exists():
        sys.exit(f"ERRO: {coords_file} não encontrado")

    genes = []
    with open(coords_file, newline='') as f:
        for row in csv.DictReader(f, delimiter='\t'):
            genes.append(row)

    print(f"Genes a processar: {len(genes)}")
    print(f"Upstream: {args.upstream} bp | Fonte: NCBI efetch (SL4.0/ITAG4.0)")
    print(f"Saída: {out_file}\n")

    ok, fail = 0, 0
    with open(out_file, "w") as out:
        for i, g in enumerate(genes, 1):
            gid    = g["gene_id"]
            chrom  = g["chrom"]
            start  = int(g["start"])
            end    = int(g["end"])
            strand = g["strand"]

            chrom_e = CHROM_MAP.get(chrom)
            if chrom_e is None:
                print(f"[{i:02d}/{len(genes)}] {gid}  AVISO: cromossomo '{chrom}' sem mapeamento")
                fail += 1
                continue

            seq, err = fetch_seq(chrom_e, start, end, strand, args.upstream)
            if err:
                print(f"[{i:02d}/{len(genes)}] {gid}  ERRO: {err}")
                fail += 1
            else:
                header = f">{gid}::{chrom}:{start}-{end}({strand}) upstream_{args.upstream}bp"
                out.write(f"{header}\n{seq}\n")
                print(f"[{i:02d}/{len(genes)}] {gid}  OK  ({len(seq)} bp)")
                ok += 1

            time.sleep(0.34)   # respeitar rate limit (~3 req/s Ensembl)

    print(f"\nConcluído: {ok} OK | {fail} falhas")
    if ok > 0:
        print(f"FASTA salvo: {out_file}")
        print(f"\nPróximo passo:")
        print(f"  Acesse https://bioinformatics.psb.ugent.be/webtools/plantcare/html/")
        print(f"  Submeta o arquivo: {out_file}")


if __name__ == "__main__":
    main()
