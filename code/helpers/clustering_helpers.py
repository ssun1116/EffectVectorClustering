import base64
import collections
import json
import math

import numpy as np


def load_plotly_heatmap(path):
    text = path.read_text()
    call_start = text.rfind("Plotly.newPlot(")
    if call_start < 0:
        raise ValueError(f"Could not find Plotly.newPlot in {path}")
    data_start = text.find("[", call_start)
    data, _ = json.JSONDecoder().raw_decode(text[data_start:])
    heatmap = next(trace for trace in data if trace.get("type") == "heatmap")
    x_labels = heatmap["x"]
    y_labels = heatmap["y"]
    z = heatmap["z"]
    if isinstance(z, dict) and "bdata" in z:
        matrix = np.frombuffer(base64.b64decode(z["bdata"]), dtype=np.dtype(z["dtype"]))
        matrix = matrix.reshape(tuple(map(int, z["shape"].split(","))))
    else:
        matrix = np.asarray(z, dtype=np.float64)

    # Plotly heatmap y-axis is reversed relative to x. Reorder rows so A[i, j]
    # is similarity between x_labels[i] and x_labels[j].
    y_pos = {label: i for i, label in enumerate(y_labels)}
    aligned = np.array([matrix[y_pos[label], :] for label in x_labels])
    return x_labels, aligned


def parse_label(label):
    cell, gene = label.rsplit(" - ", 1)
    return cell, gene


def prefix_sums(matrix):
    return np.pad(matrix, ((1, 0), (1, 0))).cumsum(axis=0).cumsum(axis=1)


def block_sum(prefix, a, b, c, d):
    return prefix[b, d] - prefix[a, d] - prefix[b, c] + prefix[a, c]


def block_mean(prefix, a, b, exclude_diag=False):
    size = b - a
    denom = size * size
    total = block_sum(prefix, a, b, a, b)
    if exclude_diag:
        total -= size
        denom -= size
    return float(total / denom) if denom > 0 else float("nan")


def cross_mean(prefix, a, b, c, d):
    denom = (b - a) * (d - c)
    return float(block_sum(prefix, a, b, c, d) / denom) if denom > 0 else float("nan")


def recursive_segments(matrix, min_size=6, max_size=80, min_gain=0.045, min_child_within=0.02):
    prefix = prefix_sums(matrix)
    segments = []

    def split(a, b):
        n = b - a
        if n <= max_size:
            segments.append((a, b))
            return

        parent = block_mean(prefix, a, b, exclude_diag=True)
        best = None
        for cut in range(a + min_size, b - min_size + 1):
            l_n = cut - a
            r_n = b - cut
            left = block_mean(prefix, a, cut, exclude_diag=True)
            right = block_mean(prefix, cut, b, exclude_diag=True)
            cross = cross_mean(prefix, a, cut, cut, b)
            child = (left * l_n + right * r_n) / n
            gain = child - cross
            score = gain + 0.25 * (child - parent)
            if best is None or score > best[0]:
                best = (score, cut, left, right, cross, child)

        if best is None:
            segments.append((a, b))
            return

        score, cut, left, right, cross, child = best
        do_split = score >= min_gain and max(left, right, child) >= min_child_within
        if do_split:
            split(a, cut)
            split(cut, b)
        else:
            segments.append((a, b))

    split(0, matrix.shape[0])
    return segments, prefix


def local_contrast(prefix, a, b, n_total):
    size = b - a
    left_a = max(0, a - size)
    right_b = min(n_total, b + size)
    vals = []
    if left_a < a:
        vals.append(cross_mean(prefix, a, b, left_a, a))
    if b < right_b:
        vals.append(cross_mean(prefix, a, b, b, right_b))
    vals = [v for v in vals if not math.isnan(v)]
    return float(np.mean(vals)) if vals else float("nan")


def entropy_from_counts(counts):
    total = sum(counts.values())
    if total == 0:
        return float("nan")
    probs = np.array([n / total for n in counts.values() if n > 0])
    return float(-(probs * np.log2(probs)).sum())


def summarize_segments(labels, matrix, segments, prefix):
    summaries = []
    member_rows = []
    for module_id, (a, b) in enumerate(segments, start=1):
        parsed = [parse_label(label) for label in labels[a:b]]
        cells = [c for c, _ in parsed]
        genes = [g for _, g in parsed]
        gene_counts = collections.Counter(genes)
        cell_counts = collections.Counter(cells)
        within = block_mean(prefix, a, b, exclude_diag=True)
        neighbor = local_contrast(prefix, a, b, len(labels))
        top_genes = gene_counts.most_common(12)
        top_cells = cell_counts.most_common(10)
        top_gene_count = top_genes[0][1] if top_genes else 0
        summaries.append(
            {
                "module_id": module_id,
                "manual_label": "",
                "start": a,
                "end": b,
                "size": b - a,
                "within_mean": within,
                "neighbor_mean": neighbor,
                "contrast": within - neighbor if not math.isnan(neighbor) else float("nan"),
                "n_unique_genes": len(gene_counts),
                "top_gene_fraction": top_gene_count / max(1, b - a),
                "gene_entropy": entropy_from_counts(gene_counts),
                "n_unique_cell_types": len(cell_counts),
                "top_genes": "; ".join(f"{g}:{n}" for g, n in top_genes),
                "top_cells": "; ".join(f"{c}:{n}" for c, n in top_cells),
            }
        )
        for pos, label in enumerate(labels[a:b], start=a):
            cell, gene = parse_label(label)
            member_rows.append(
                {
                    "module_id": module_id,
                    "position": pos,
                    "label": label,
                    "cell_type": cell,
                    "gene": gene,
                }
            )
    return summaries, member_rows
