"""Run low-resolution module GO enrichment and annotate the module browser."""

import csv
import json
import pathlib
import re
import time

import requests


ROOT = pathlib.Path(__file__).resolve().parents[1]
OUTPUTS = ROOT / "analysis_outputs"
TABLES = OUTPUTS / "tables"
HTML = OUTPUTS / "html"

MODULE_SUMMARY = TABLES / "low_resolution_module_summary.tsv"
MODULE_MEMBERS = TABLES / "low_resolution_module_members.tsv"
MODULE_BROWSER = HTML / "low_resolution_module_browser.html"
GO_SUMMARY = TABLES / "low_resolution_module_go_summary.tsv"
GO_ENRICHMENT = TABLES / "low_resolution_module_go_enrichment.tsv"
GO_BROWSER = HTML / "low_resolution_go_module_browser.html"

GPROFILER_URL = "https://biit.cs.ut.ee/gprofiler/api/gost/profile/"
GO_SOURCES = ["GO:BP", "GO:MF", "GO:CC"]


def read_tsv(path):
    with path.open(newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def write_tsv(path, rows, columns):
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=columns, delimiter="\t", extrasaction="ignore")
        writer.writeheader()
        writer.writerows({column: row.get(column, "") for column in columns} for row in rows)


def run_gprofiler(query_genes, background_genes, retries=3):
    if len(query_genes) < 2:
        return []
    payload = {
        "organism": "mmusculus",
        "query": query_genes,
        "sources": GO_SOURCES,
        "user_threshold": 1.0,
        "all_results": True,
        "domain_scope": "custom",
        "background": background_genes,
        "no_evidences": False,
        "significance_threshold_method": "fdr",
    }
    error = None
    for attempt in range(retries):
        try:
            response = requests.post(GPROFILER_URL, json=payload, timeout=60)
            response.raise_for_status()
            return response.json().get("result", [])
        except Exception as exception:
            error = exception
            time.sleep(attempt + 1)
    raise RuntimeError(f"g:Profiler failed for module genes {query_genes[:5]}") from error


def intersection_genes(result, query_genes):
    evidence = result.get("intersections") or []
    return ",".join(gene for gene, hit in zip(query_genes, evidence) if hit)


def add_go_to_browser(labels):
    text = MODULE_BROWSER.read_text()
    match = re.search(
        r'<script id="heatmap-data" type="application/json">(.*?)</script>',
        text,
        re.S,
    )
    if not match:
        raise ValueError(f"Embedded browser data not found in {MODULE_BROWSER}")
    data = json.loads(match.group(1))
    for module in data["modules"]:
        go = labels.get(str(module["module_id"]), {})
        module["auto_theme"] = go.get("go_auto_theme", "")
        module["theme_hits"] = go.get("go_term_id", "")
        module["go_source"] = go.get("go_source", "")
        module["go_p_value"] = go.get("go_p_value", "")
        module["go_significant"] = go.get("go_significant", "")
    encoded = json.dumps(data, separators=(",", ":")).replace("</", "<\\/")
    output = text[: match.start(1)] + encoded + text[match.end(1) :]
    output = re.sub(r"<title>.*?</title>", "<title>Low-Resolution GO Module Browser</title>", output, count=1)
    output = re.sub(r"<h1>.*?</h1>", "<h1>Low-Resolution GO Module Browser</h1>", output, count=1)
    output = output.replace("<th>auto theme</th>", "<th>top GO term</th>")
    output = output.replace("Theme: ${esc(m.auto_theme)}", "top GO term: ${esc(m.auto_theme)}")
    GO_BROWSER.write_text(output)


def main():
    summaries = read_tsv(MODULE_SUMMARY)
    members = read_tsv(MODULE_MEMBERS)
    genes_by_module = {}
    for row in members:
        if row["gene"]:
            genes_by_module.setdefault(row["module_id"], set()).add(row["gene"])
    background = sorted({row["gene"] for row in members if row["gene"]})

    labels = {}
    enrichment_rows = []
    for summary in summaries:
        module_id = summary["module_id"]
        genes = sorted(genes_by_module.get(module_id, set()))
        results = sorted(
            run_gprofiler(genes, background),
            key=lambda result: (result.get("p_value", 1), result.get("source", "")),
        )
        top = results[0] if results else {}
        labels[module_id] = {
            "go_auto_theme": top.get("name", ""),
            "go_term_id": top.get("native", ""),
            "go_source": top.get("source", ""),
            "go_p_value": f"{top.get('p_value'):.6g}" if isinstance(top.get("p_value"), (int, float)) else "",
            "go_intersection_size": top.get("intersection_size", ""),
            "go_term_size": top.get("term_size", ""),
            "go_significant": top.get("significant", ""),
        }
        for rank, result in enumerate(results[:20], start=1):
            enrichment_rows.append({
                "module_id": module_id,
                "rank": rank,
                "source": result.get("source", ""),
                "term_id": result.get("native", ""),
                "term_name": result.get("name", ""),
                "p_value": result.get("p_value", ""),
                "significant": result.get("significant", ""),
                "query_size": result.get("query_size", ""),
                "term_size": result.get("term_size", ""),
                "intersection_size": result.get("intersection_size", ""),
                "precision": result.get("precision", ""),
                "recall": result.get("recall", ""),
                "intersection_genes": intersection_genes(result, genes),
                "description": result.get("description", ""),
            })

    go_summaries = [dict(row, **labels[row["module_id"]]) for row in summaries]
    go_columns = list(summaries[0]) + [
        "go_auto_theme", "go_term_id", "go_source", "go_p_value",
        "go_intersection_size", "go_term_size", "go_significant",
    ]
    enrichment_columns = [
        "module_id", "rank", "source", "term_id", "term_name", "p_value",
        "significant", "query_size", "term_size", "intersection_size",
        "precision", "recall", "intersection_genes", "description",
    ]
    write_tsv(GO_SUMMARY, go_summaries, go_columns)
    write_tsv(GO_ENRICHMENT, enrichment_rows, enrichment_columns)
    add_go_to_browser(labels)
    print(f"GO enrichment completed for {len(summaries)} modules")
    print(GO_SUMMARY)
    print(GO_ENRICHMENT)
    print(GO_BROWSER)


if __name__ == "__main__":
    main()
