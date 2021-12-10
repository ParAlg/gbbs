"""
This is a script that takes all the outputs in ~/scan-experiment-results and
organizes the data into several files in ~/scan-experiments-results/summary
representing the data for each plot in the paper "Parallel Index-Based
Structural Graph Clustering and Its Approximation".
"""
#!/usr/bin/env python3.6

import collections
import csv
import itertools
import os

INPUT_DIRECTORY = os.path.expanduser("~/scan-experiment-results")
OUTPUT_DIRECTORY = INPUT_DIRECTORY + "/summary"
if not os.path.exists(OUTPUT_DIRECTORY):
    os.makedirs(OUTPUT_DIRECTORY)

def get_all_data():
    all_data = []
    for filename in os.listdir(INPUT_DIRECTORY):
        if not filename.endswith(".txt"):
            continue
        split_filename = filename[:-4].split("-")
        if len(split_filename) < 2:
            print("ignoring file " + filename)
            continue

        path = INPUT_DIRECTORY + "/" + filename
        graph = split_filename[0]
        is_mkl = len(split_filename) >= 3 and split_filename[2] == "mkl"
        is_serial = split_filename[-1] == "serial"
        if split_filename[1] == "gbbs":
            algorithm = "GBBSIndexSCAN"
            if is_mkl:
                algorithm += "_matrix-multiply"
            if is_serial:
                algorithm += "_1-thread"

            with open(path, newline="") as csvfile:
                truncated_file = itertools.dropwhile(
                    lambda line: not line.startswith("BEGIN GBBS EXPERIMENTS OUTPUT"),
                    csvfile,
                )
                next(truncated_file)
                reader = csv.DictReader(truncated_file)
                for row in reader:
                    row["graph"] = graph
                    row["algorithm"] = algorithm
                    all_data.append(dict(row))
        elif split_filename[1] == "pp":
            with open(path, newline="") as csvfile:
                reader = csv.DictReader(csvfile)
                for row in reader:
                    row["graph"] = graph
                    row["algorithm"] = "ppSCAN"
                    row["similarity measure"] = "cosine similarity"
                    row["operation"] = "cluster"
                    all_data.append(row)
            pass
        else:
            print("ignoring file " + filename)
    return all_data


FIELD_NAMES = [
    "graph",
    "algorithm",
    "similarity measure",
    "approximation samples",
    "operation",
    "mu",
    "epsilon",
    "median time",
    "quality",
]


def write_summary(filename, rows):
    with open(OUTPUT_DIRECTORY + "/" + filename, "w", newline="") as csvfile:
        writer = csv.DictWriter(csvfile, fieldnames=FIELD_NAMES)
        writer.writeheader()
        for row in rows:
            writer.writerow(row)


# Figure 5: Exact index construction times
def fig5(data):
    rows = []
    for row in data:
        if row["operation"] == "index construction" and not row[
            "similarity measure"
        ].startswith("approx"):
            rows.append(row)
    rows.sort(key=lambda r: (r["graph"], r["algorithm"]))
    write_summary("figure-5.txt", rows)


# Figure 6: Clustering time with mu=5 and varying epsilon
def fig6(data):
    rows = []
    for row in data:
        if row["operation"] == "cluster" and int(row["mu"]) == 5:
            rows.append(row)
    rows.sort(key=lambda r: (r["graph"], r["algorithm"], float(r["epsilon"])))
    write_summary("figure-6.txt", rows)


# Figure 7: Clustering time with mu=5 and varying epsilon
def fig7(data):
    rows = []
    for row in data:
        if (
            row["operation"] == "cluster"
            and int(row["mu"]) != 5
            and abs(float(row["epsilon"]) - 0.6) < 1e-3
        ):
            rows.append(row)
    rows.sort(key=lambda r: (r["graph"], r["algorithm"], int(r["mu"])))
    write_summary("figure-7.txt", rows)


# Figure 8: Approximate index construction time
def fig8(data):
    rows = []
    for row in data:
        if (
            row["operation"] == "index construction"
            and row["algorithm"] == "GBBSIndexSCAN"
        ):
            rows.append(row)
    rows.sort(
        key=lambda r: (
            r["graph"],
            r["algorithm"],
            r["similarity measure"],
            int(r["approximation samples"] or 0),
        )
    )
    write_summary("figure-8.txt", rows)


# Figure 9: Approximate index construction time vs. modularity
def fig9(data):
    Key = collections.namedtuple("Key", ["graph", "measure", "samples"])
    points = {}
    for row in data:
        if row["algorithm"] != "GBBSIndexSCAN":
            continue
        key = Key(row["graph"], row["similarity measure"], row["approximation samples"])
        val = None
        if row["operation"] == "index construction":
            val = {"median time": row["median time"]}
        elif row["operation"] == "best exact clustering":
            val = {"quality": row["quality"]}
        elif row["operation"] == "mean modularity at best approximate clustering":
            val = {"quality": row["quality"]}
        if not val:
            continue
        points[key] = {**points.get(key, {}), **val}

    rows = []
    for key, val in points.items():
        if key.measure.startswith("jaccard"):
            cosine_key = key._replace(measure="cosine similarity")
            val["median time"] = points[cosine_key]["median time"]

        rows.append(
            {
                "graph": key.graph,
                "similarity measure": key.measure,
                "approximation samples": key.samples,
                **val,
            }
        )
    rows.sort(
        key=lambda r: (
            r["graph"],
            r["similarity measure"],
            int(r["approximation samples"] or 0),
        )
    )
    write_summary("figure-9.txt", rows)


# Figure 10: Approximate index construction time vs. ARI
def fig10(data):
    Key = collections.namedtuple("Key", ["graph", "measure", "samples"])
    points = {}
    for row in data:
        if row["algorithm"] != "GBBSIndexSCAN":
            continue
        key = Key(row["graph"], row["similarity measure"], row["approximation samples"])
        val = None
        if row["operation"] == "index construction":
            val = {"median time": row["median time"]}
            if not row["similarity measure"].startswith("approx"):
                val["quality"] = 1
        elif row["operation"] == "mean ari at best exact clustering":
            val = {"quality": row["quality"]}
        if not val:
            continue
        points[key] = {**points.get(key, {}), **val}

    rows = []
    for key, val in points.items():
        rows.append(
            {
                "graph": key.graph,
                "similarity measure": key.measure,
                "approximation samples": key.samples,
                **val,
            }
        )
    rows.sort(
        key=lambda r: (
            r["graph"],
            r["similarity measure"],
            int(r["approximation samples"] or 0),
        )
    )
    write_summary("figure-10.txt", rows)


if __name__ == "__main__":
    data = get_all_data()
    fig5(data)
    fig6(data)
    fig7(data)
    fig8(data)
    fig9(data)
    fig10(data)
