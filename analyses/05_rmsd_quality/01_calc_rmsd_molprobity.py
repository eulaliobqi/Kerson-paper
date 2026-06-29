#!/usr/bin/env python3
"""
Extrai RMSD numérico em Å entre os 7 modelos 3D de LRR-RLPs e verifica qualidade.

Modo 1 (PyMOL batch): calcula RMSD pairwise entre todas as estruturas
Modo 2 (MolProbity): prepara sequência de submissão ao MolProbity web

Uso:
  python3 01_calc_rmsd_molprobity.py --pdb-dir ./pdb_models/
  python3 01_calc_rmsd_molprobity.py --pdb-dir ./pdb_models/ --mode molprobity

Dependência para modo PyMOL:
  mamba install -c conda-forge pymol-open-source
"""

import os, sys, glob, argparse
from pathlib import Path

GENES = {
    "Solyc05g055190": "SlRLP1_CLV2",
    "Solyc03g112680":  "SlRLP2",
    "Solyc05g009990":  "SlRLP3_RIC7",
    "Solyc12g042760":  "SlRLP4_TMM",
    "Solyc02g072250":  "SlRLP5_SNC23",
    "Solyc02g092040":  "SlRLP6",
    "Solyc10g007830":  "SlRLP7",
}

# ── Modo 1: PyMOL RMSD ────────────────────────────────────────────────────────
PYMOL_SCRIPT = """#!/usr/bin/env python3
# PyMOL script — calcular RMSD pairwise entre 7 modelos LRR-RLP
# Executar: pymol -c rmsd_calc.py
# Saída: rmsd_matrix.csv

import pymol
from pymol import cmd
import csv, os, sys

pdb_dir = sys.argv[1] if len(sys.argv) > 1 else "."
genes = {GENE_DICT}

# Carregar todos os modelos
for gene_id, label in genes.items():
    pdb_file = os.path.join(pdb_dir, f"{{gene_id}}.pdb")
    if not os.path.exists(pdb_file):
        # Tentar nomes alternativos
        alternatives = glob.glob(os.path.join(pdb_dir, f"{{gene_id}}*.pdb"))
        if alternatives:
            pdb_file = alternatives[0]
        else:
            print(f"AVISO: PDB não encontrado para {{gene_id}}")
            continue
    cmd.load(pdb_file, label)
    cmd.select(f"{{label}}_CA", f"{{label}} and name CA")
    print(f"Carregado: {{label}}")

labels = list(genes.values())
cmd.align("all")  # alinhamento inicial de referência

# Calcular RMSD pairwise
rows = []
print("\\nMatrix de RMSD (Å) — alinhamento Ca:")
header = [""] + labels
print("\\t".join(header))

for i, lbl1 in enumerate(labels):
    row = [lbl1]
    for j, lbl2 in enumerate(labels):
        if i == j:
            row.append("0.000")
        elif i < j:
            try:
                result = cmd.align(f"{{lbl1}} and name CA",
                                   f"{{lbl2}} and name CA",
                                   cycles=0)
                rmsd = round(result[0], 3)
            except Exception as e:
                rmsd = "NA"
            row.append(str(rmsd))
        else:
            # Usar valor simétrico já calculado
            row.append(rows[j][i+1])
    rows.append(row)
    print("\\t".join(row))

# Salvar CSV
with open("rmsd_matrix.csv", "w", newline="") as f:
    writer = csv.writer(f)
    writer.writerow(header)
    writer.writerows(rows)

print("\\nMatriz RMSD salva: rmsd_matrix.csv")

# Obter pLDDT médio de modelos AlphaFold (b-factor)
print("\\n=== pLDDT médio por modelo ===")
for label in labels:
    try:
        stored = []
        cmd.iterate(f"{{label}} and name CA", "stored.append(b)", space={{"stored": stored}})
        if stored:
            mean_plddt = sum(stored)/len(stored)
            print(f"  {{label}}: {{mean_plddt:.1f}}")
    except:
        pass

cmd.quit()
"""

# ── Modo 2: MolProbity ───────────────────────────────────────────────────────
def generate_molprobity_instructions(pdb_dir: str, outfile: str):
    """Gera instruções + script curl para submissão ao MolProbity."""
    pdbs = sorted(glob.glob(os.path.join(pdb_dir, "*.pdb")))

    with open(outfile, "w") as f:
        f.write("# Instruções para MolProbity — Qualidade dos Modelos 3D\n")
        f.write("# MolProbity: https://molprobity.biochem.duke.edu/\n\n")
        f.write("# Para cada PDB abaixo:\n")
        f.write("#   1. Acesse: https://molprobity.biochem.duke.edu/index.php\n")
        f.write("#   2. Upload do arquivo PDB\n")
        f.write("#   3. Run 'Multi-criterion chart'\n")
        f.write("#   4. Anotar: Ramachandran favored (%), rotamer outliers (%),\n")
        f.write("#              Clashscore, MolProbity score\n\n")
        f.write("# Métricas alvo para publicação (Plant Journal / Frontiers):\n")
        f.write("#   Ramachandran favored: > 95%\n")
        f.write("#   Rotamer outliers: < 2%\n")
        f.write("#   Clashscore: < 20 (idealmente < 10)\n\n")
        f.write("gene_id\tfile\trama_favored_%\trotamer_outliers_%\tclashscore\tmolprobity_score\n")
        for pdb in pdbs:
            basename = os.path.basename(pdb)
            gene = basename.replace(".pdb","").split("_")[0]
            f.write(f"{gene}\t{basename}\t\t\t\t\n")

    print(f"Template MolProbity salvo: {outfile}")
    print("Preencher manualmente após submissão ao site.")


# ── Main ─────────────────────────────────────────────────────────────────────
def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--pdb-dir", default="./pdb_models",
                        help="Diretório com os 7 arquivos PDB dos modelos 3D")
    parser.add_argument("--mode", choices=["pymol","molprobity","both"],
                        default="both")
    args = parser.parse_args()

    pdb_dir  = Path(args.pdb_dir)
    out_dir  = Path(__file__).parent

    pdbs = list(pdb_dir.glob("*.pdb"))
    print(f"PDBs encontrados em {pdb_dir}: {len(pdbs)}")
    for p in pdbs:
        print(f"  {p.name}")

    if len(pdbs) == 0:
        print(f"\nAVISO: Nenhum PDB em {pdb_dir}")
        print("Coloque os 7 arquivos PDB (Swiss-Model ou AlphaFold3) em:")
        print(f"  {pdb_dir.resolve()}/")
        print("  Nomeie como: Solyc05g055190.pdb, Solyc03g112680.pdb, ...")

    # Modo PyMOL
    if args.mode in ("pymol","both"):
        gene_dict_str = repr(GENES)
        script = PYMOL_SCRIPT.replace("{GENE_DICT}", gene_dict_str)
        pymol_script_path = out_dir / "rmsd_calc.py"
        pymol_script_path.write_text(script)
        print(f"\nScript PyMOL gerado: {pymol_script_path}")
        print(f"Executar: pymol -c {pymol_script_path} {pdb_dir.resolve()}")

    # Modo MolProbity
    if args.mode in ("molprobity","both"):
        mp_out = out_dir / "molprobity_template.tsv"
        generate_molprobity_instructions(str(pdb_dir), str(mp_out))

    print("\n=== Onde obter os PDBs ===")
    print("Swiss-Model: https://swissmodel.expasy.org/")
    print("AlphaFold3:  https://alphafoldserver.com/")
    print("AlphaFold DB (pré-computado): https://alphafold.ebi.ac.uk/")
    print("  Buscar por: Solyc05g055190 (tomate ITAG4.0)")


if __name__ == "__main__":
    main()
