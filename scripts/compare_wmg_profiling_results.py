import argparse
import ast
import datetime
import json
import os
import sys

import matplotlib.pyplot as plt

# Add the root directory to the Python module search path so you can reference backend
# without needing to move this script to the root directory to run it.
scripts_dir = os.path.dirname(os.path.abspath(__file__))
root_dir = os.path.dirname(scripts_dir)
sys.path.append(root_dir)


def generate_scatter_plots(json1, json2):
    # Load the profiling results
    with open(json1, "r") as f1, open(json2, "r") as f2:
        profiling_results1 = json.load(f1)["profiles"]
        profiling_results2 = json.load(f2)["profiles"]

    # Create a new figure for each metric
    metrics = ["query-tiledb", "build-response", "total_time_backend"]
    compare_options = {d["params"]["compare"] for d in profiling_results1}
    sex_criteria = [
        ast.literal_eval(i)
        for i in {str(d["params"]["payload"]["filter"]["sex_ontology_term_ids"]) for d in profiling_results1}
    ]
    num_conditions = len(compare_options) * len(sex_criteria)

    fig, axs = plt.subplots(len(metrics), num_conditions)
    fig.set_size_inches(4 * num_conditions, 12)

    # Generate scatter plots for the collected metrics
    for i, metric in enumerate(metrics):
        for j, compare_option in enumerate(compare_options):
            for k, sex in enumerate(sex_criteria):
                metric_values1 = [
                    d[metric]
                    for d in profiling_results1
                    if d["params"]["compare"] == compare_option
                    and d["params"]["payload"]["filter"]["sex_ontology_term_ids"] == sex
                ]
                metric_values2 = [
                    d[metric]
                    for d in profiling_results2
                    if d["params"]["compare"] == compare_option
                    and d["params"]["payload"]["filter"]["sex_ontology_term_ids"] == sex
                ]

                max_value = max(metric_values1 + metric_values2)
                row, col = i, j * len(sex_criteria) + k
                axs[row, col].plot([0, max_value], [0, max_value], "k:")
                axs[row, col].scatter(metric_values1, metric_values2)
                axs[row, col].set_title(f"compare: {compare_option}, sex: {sex}")
                axs[row, col].set_xlabel("Run 1")
                axs[row, col].set_ylabel("Run 2")
                axs[row, col].set_xlim(0, max_value + 0.5)
                axs[row, col].set_ylim(0, max_value + 0.5)

    # Manually add super-titles for each row
    for i, metric in enumerate(metrics):
        # Calculate y position for the super-title of each row
        y_position = 1 - (i - 0.01) / len(metrics)
        fig.text(0.5, y_position, metric, ha="center", va="center", fontsize=14, weight="bold")

    fig.tight_layout(pad=3.0)
    plt.rcParams.update({"font.size": 10})

    # Save the scatter plots to a PDF
    timestamp = datetime.datetime.now().strftime("%Y%m%d_%H%M%S")
    fig.savefig(f"scatter_plots_{timestamp}.pdf")
    print(f"Scatter plots output to scatter_plots_{timestamp}.pdf")


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("json1", help="The first profiling results JSON file")
    parser.add_argument("json2", help="The second profiling results JSON file")
    args = parser.parse_args()
    generate_scatter_plots(args.json1, args.json2)
