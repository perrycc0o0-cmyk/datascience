"""Build the excitatory IEG-high gene-category membership network."""

import argparse
import csv
import json
import math
from collections import defaultdict
from pathlib import Path


def read_table(path):
    with path.open(encoding="utf-8-sig", newline="") as handle:
        return list(csv.DictReader(handle))


def write_table(path, rows, columns):
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=columns)
        writer.writeheader()
        writer.writerows(rows)


def number(value):
    return math.nan if value in ("", "NA") else float(value)


def direction(row):
    return "Upregulated" if number(row["log2FoldChange"]) > 0 else "Downregulated"


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("deseq2_output_dir", type=Path)
    parser.add_argument("gsea_output_dir", type=Path)
    parser.add_argument("output_dir", type=Path)
    args = parser.parse_args()
    key = "Excitatory__RetExt_vs_Ext"
    de = read_table(args.deseq2_output_dir / "datasets" / key / "tables" / "ALL_genes.csv")
    ranks = read_table(args.gsea_output_dir / "source_data" / f"{key}_ranked_genes.csv")
    tests = read_table(args.gsea_output_dir / "tables" / f"{key}_P005_displayed_GSEA.csv")
    membership = read_table(args.gsea_output_dir / "source_data" / "union_gene_set_membership.csv")
    config = read_table(Path(__file__).resolve().parent / "config" / "network_categories.csv")

    pathways = {row["ID"]: row for row in tests}
    if set(pathways) != {row["pathway_ID"] for row in config} or not all(
        0 <= number(row["pvalue"]) < 0.05 for row in tests
    ):
        raise ValueError("The GSEA results do not match the frozen nominal-P pathway panel.")
    de_by_gene = {row["gene"]: row for row in de}
    rank_by_gene = {row["gene"]: row for row in ranks}
    if len(de_by_gene) != len(de) or len(rank_by_gene) != len(ranks):
        raise ValueError("Feature identifiers must be unique.")
    for row in ranks:
        native = de_by_gene[row["gene"]]
        if not all(math.isclose(number(row[r]), number(native[d]), rel_tol=1e-12, abs_tol=tol)
                   for r, d, tol in [("DE_pvalue", "pvalue", 1e-300),
                                     ("log2FoldChange", "log2FoldChange", 1e-15)]):
            raise ValueError("Rank and DE records differ for feature " + row["gene"])

    gene_sets = defaultdict(set)
    for row in membership:
        gene_sets[row["ID"]].add(row["Entrez"])
    categories = {}
    dictionary = []
    for row in config:
        cid, label, pid = (row[k] for k in ("category_id", "category_name", "pathway_ID"))
        categories.setdefault(cid, (label, []))[1].append(pid)
        test = pathways[pid]
        dictionary.append(dict(category_id=cid, category_name=label, pathway_ID=pid,
                               pathway_name=test["Description"], GSEA_NES=test["NES"],
                               GSEA_P=test["pvalue"], GSEA_FDR=test["BH"]))
    selected_de = [row for row in de if 0 <= number(row["pvalue"]) < 0.05]
    links = defaultdict(list)
    edges = []
    for cid, (label, ids) in categories.items():
        for row in selected_de:
            gene = row["gene"]
            mapped = rank_by_gene.get(gene)
            if not mapped:
                continue
            entrez = mapped["ENTREZID"]
            matches = [pid for pid in ids if entrez in gene_sets[pid]]
            if not matches:
                continue
            links[gene].append(cid)
            edges.append(dict(
                direction=direction(row), gene=gene, symbol=mapped["symbol"],
                ENTREZID=entrez, category_id=cid, category_name=label,
                member_pathway_IDs=";".join(matches),
                member_pathway_names=";".join(pathways[pid]["Description"] for pid in matches),
                log2FoldChange=row["log2FoldChange"], P=row["pvalue"], BH=row["padj"],
                edge_type="gene-category membership"))

    node_columns = ["direction", "node_id", "label", "node_type", "log2FC", "count", "P", "BH", "categories"]
    edge_columns = ["direction", "gene", "symbol", "ENTREZID", "category_id", "category_name",
                    "member_pathway_IDs", "member_pathway_names", "log2FoldChange", "P", "BH", "edge_type"]
    nodes = []
    for side in ("Downregulated", "Upregulated"):
        subset = [row for row in edges if row["direction"] == side]
        for gene in sorted({row["gene"] for row in subset}):
            row = de_by_gene[gene]
            nodes.append(dict(direction=side, node_id=gene, label=gene, node_type="gene",
                              log2FC=row["log2FoldChange"], count=1, P=row["pvalue"],
                              BH=row["padj"], categories=";".join(links[gene])))
        for cid, (label, _) in categories.items():
            values = [number(row["log2FoldChange"]) for row in subset if row["category_id"] == cid]
            if values:
                nodes.append(dict(direction=side, node_id=cid, label=label, node_type="category",
                                  log2FC=sum(values) / len(values), count=len(values), P="", BH="", categories=cid))

    # The isolated two-gene upregulated component is omitted from the display only.
    # Do not remove it if a different input no longer has this exact isolated component.
    c6 = [row for row in edges if row["direction"] == "Upregulated" and row["category_id"] == "C6"]
    removed = {"Inhba", "Hmgcs1"}
    omit = {row["gene"] for row in c6} == removed and all(links[gene] == ["C6"] for gene in removed)
    display_nodes = [row for row in nodes if not (
        omit and row["direction"] == "Upregulated" and row["node_id"] in removed | {"C6"})]
    display_edges = [row for row in edges if not (
        omit and row["direction"] == "Upregulated" and row["category_id"] == "C6")]
    audits = []
    for row in selected_de:
        gene = row["gene"]
        mapped = rank_by_gene.get(gene)
        memberships = links.get(gene, [])
        status = ("included" if memberships else "no_membership_in_selected_pathways" if mapped
                  else "not_in_frozen_GSEA_mapped_rank_universe")
        audits.append(dict(gene=gene, direction=direction(row), log2FoldChange=row["log2FoldChange"],
                           P=row["pvalue"], BH=row["padj"], ENTREZID=mapped["ENTREZID"] if mapped else "",
                           categories=";".join(memberships), included=bool(memberships), status=status))

    out = args.output_dir
    if out.exists() and any(out.iterdir()):
        raise ValueError(f"Choose an empty output directory: {out}")
    out.mkdir(parents=True, exist_ok=True)
    write_table(out / "network_nodes.csv", nodes, node_columns)
    write_table(out / "network_edges.csv", edges, edge_columns)
    write_table(out / "display_nodes.csv", display_nodes, node_columns)
    write_table(out / "display_edges.csv", display_edges, edge_columns)
    write_table(out / "category_pathway_dictionary.csv", dictionary,
                ["category_id", "category_name", "pathway_ID", "pathway_name", "GSEA_NES", "GSEA_P", "GSEA_FDR"])
    write_table(out / "nominal_DE_inclusion_audit.csv", audits,
                ["gene", "direction", "log2FoldChange", "P", "BH", "ENTREZID", "categories", "included", "status"])
    summary = dict(nominal_DE_genes=len(selected_de),
                   mapped_nominal_DE_genes=sum(row["gene"] in rank_by_gene for row in selected_de),
                   included_genes=len(links), membership_edges=len(edges), display_edges=len(display_edges),
                   isolated_component_display_omitted=omit)
    (out / "network_summary.json").write_text(json.dumps(summary, indent=2) + "\n", encoding="utf-8")


if __name__ == "__main__":
    main()
