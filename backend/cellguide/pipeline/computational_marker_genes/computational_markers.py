import itertools
import pathlib

import numpy as np
import pandas as pd

from backend.cellguide.pipeline.computational_marker_genes.utils import (
    post_process_stats,
    run_ttest,
)
from backend.common.utils.rollup import rollup_across_cell_type_descendants, rollup_across_cell_type_descendants_array
from backend.wmg.data.snapshot import WmgSnapshot

ENSEMBL_GENE_ID_TO_DESCRIPTION_FILENAME = ""


class MarkerGenesCalculator:
    def __init__(self, snapshot: WmgSnapshot, all_cell_type_ids_in_corpus: list[str], groupby_terms: list[str]):
        self.all_cell_type_ids_in_corpus = all_cell_type_ids_in_corpus
        self.organism_id_to_name = {
            k: v for d in snapshot.primary_filter_dimensions["organism_terms"] for k, v in d.items()
        }
        self.tissue_id_to_name = {
            k: v
            for organism in snapshot.primary_filter_dimensions["tissue_terms"]
            for i in snapshot.primary_filter_dimensions["tissue_terms"][organism]
            for k, v in i.items()
        }

        file_path = self._get_symbol_to_gene_descriptions_file_path()
        gene_metadata = pd.read_csv(file_path, sep="\t")
        self.gene_id_to_name = gene_metadata.set_index("Ensembl GeneIDs")["Description"].to_dict()
        self.gene_id_to_symbol = gene_metadata.set_index("Ensembl GeneIDs")["Symbols"].to_dict()

        if "cell_type_ontology_term_id" in groupby_terms:
            groupby_terms.remove("cell_type_ontology_term_id")
        self.groupby_terms = groupby_terms
        self.group_by_terms_with_celltype = groupby_terms + ["cell_type_ontology_term_id"]
        self.group_by_terms_with_celltype_and_gene = self.group_by_terms_with_celltype + ["cell_type_ontology_term_id"]

        cell_counts_df = snapshot.cell_counts_cube.df[:]
        expressions_df = snapshot.expression_summary_fmg_cube.df[:]
        (self.universe_cell_counts_df, self.expressions_df) = self._prepare_cell_counts_and_gene_expression_dfs(
            cell_counts_df, expressions_df
        )

    def _get_symbol_to_gene_descriptions_file_path(self) -> pathlib.Path:
        """
        Returns the file path for the table mapping gene symbols to gene descriptions.
        Returns
        -------
        pathlib.Path: The file path for the ensembl to gene descriptions file.
        """

        return (
            pathlib.Path(__file__)
            .parent.absolute()
            .parent.joinpath("fixtures", ENSEMBL_GENE_ID_TO_DESCRIPTION_FILENAME)
        )

    def _prepare_cell_counts_and_gene_expression_dfs(self, cell_counts_df, expressions_df):
        # group by cell counts
        cell_counts_df = cell_counts_df.groupby(self.groupby_terms).sum(numeric_only=True)
        cell_counts_df = cell_counts_df[
            cell_counts_df.index.get_level_values("cell_type_ontology_term_id").isin(self.all_cell_type_ids_in_corpus)
        ]

        # group by gene expressions
        expressions_df = expressions_df.groupby(self.groupby_terms + ["gene_ontology_term_id"]).sum(numeric_only=True)
        expressions_df = expressions_df.reset_index()
        expressions_df = expressions_df[
            expressions_df["cell_type_ontology_term_id"].isin(self.all_cell_type_ids_in_corpus)
        ]

        # rollup cell counts
        index = pd.Index(
            list(
                itertools.product(
                    *[cell_counts_df.index.get_level_values(term).unique() for term in self.groupby_terms]
                    + [self.all_cell_type_ids_in_corpus]
                )
            )
        )
        index = index.set_names(self.group_by_terms_with_celltype)
        universe_cell_counts_df = pd.DataFrame(index=index)
        universe_cell_counts_df["n_cells"] = 0
        universe_cell_counts_df["n_cells"][cell_counts_df.index] = cell_counts_df["n_cells"]

        universe_cell_counts_df = rollup_across_cell_type_descendants(universe_cell_counts_df.reset_index())
        universe_cell_counts_df = universe_cell_counts_df[universe_cell_counts_df["n_cells"] > 0]
        universe_cell_counts_df = universe_cell_counts_df.groupby(self.groupby_terms).sum()
        return universe_cell_counts_df, expressions_df

    def get_computational_marker_genes(self):
        universe_cell_counts_df = self.universe_cell_counts_df

        x = list(zip(*self.expressions_df[self.group_by_terms_with_celltype].values.T))
        y = list(self.expressions_df["gene_ontology_term_id"])

        xu = list(set(x))
        yu = list(set(y))

        x_index = pd.Series(index=pd.Index(xu), data=np.arange(len(xu)))
        y_index = pd.Series(index=pd.Index(yu), data=np.arange(len(yu)))

        e_nnz = np.zeros((len(xu), len(yu)))
        e_sum = np.zeros((len(xu), len(yu)))
        e_sqsum = np.zeros((len(xu), len(yu)))

        e_nnz[x_index[x].values, y_index[y].values] = self.expressions_df["nnz"]
        e_sum[x_index[x].values, y_index[y].values] = self.expressions_df["sum"]
        e_sqsum[x_index[x].values, y_index[y].values] = self.expressions_df["sqsum"]

        available_combinations = set(universe_cell_counts_df.index.values)
        missing_combinations = available_combinations.difference(xu)

        xu = xu + list(missing_combinations)
        x_index = pd.Series(index=pd.Index(xu), data=np.arange(len(xu)))

        e_nnz = np.vstack((e_nnz, np.zeros((len(missing_combinations), e_nnz.shape[1]))))
        e_sum = np.vstack((e_sum, np.zeros((len(missing_combinations), e_sum.shape[1]))))
        e_sqsum = np.vstack((e_sqsum, np.zeros((len(missing_combinations), e_sqsum.shape[1]))))

        groupby_term_to_values = [x_index.index.get_level_values(i) for i in range(len(self.groupby_terms))]
        groupby_term_to_unique_values = [list(set(i)) for i in groupby_term_to_values]
        cell_types = x_index.index.get_level_values(-1)

        e_nnz_rollup = np.zeros_like(e_nnz)
        e_sum_rollup = np.zeros_like(e_sum)
        e_sqsum_rollup = np.zeros_like(e_sqsum)

        n_cells = universe_cell_counts_df["n_cells"][x_index.index]
        n_cells = np.tile(n_cells.values[:, None], (1, e_nnz_rollup.shape[1]))

        all_results = []
        for combination in itertools.product(*groupby_term_to_unique_values.values()):
            filt = groupby_term_to_values[0] == combination[0]
            for _i in range(1, len(combination)):
                filt = np.logical_and(filt, groupby_term_to_values[1] == combination[1])

            if filt.sum() == 0:
                continue

            cell_types_o = cell_types[filt]

            e_nnz_o = e_nnz[filt]
            e_sum_o = e_sum[filt]
            e_sqsum_o = e_sqsum[filt]
            n_cells_o = n_cells[filt]

            e_nnz_o = rollup_across_cell_type_descendants_array(e_nnz_o, cell_types_o)
            e_sum_o = rollup_across_cell_type_descendants_array(e_sum_o, cell_types_o)
            e_sqsum_o = rollup_across_cell_type_descendants_array(e_sqsum_o, cell_types_o)

            e_nnz_rollup[filt] = e_nnz_o
            e_sum_rollup[filt] = e_sum_o
            e_sqsum_rollup[filt] = e_sqsum_o

            for i in range(e_sum_o.shape[0]):
                sum1 = e_sum_o[i][None, :].copy()
                sumsq1 = e_sqsum_o[i][None, :].copy()
                n1 = n_cells_o[i][None, :].copy()

                pvals, effects = run_ttest(sum1, sumsq1, n1, e_sum_o, e_sqsum_o, n_cells_o)

                pvals[i] = np.nan
                effects[i] = np.nan

                res = post_process_stats(
                    cell_types_o[i], cell_types_o, y_index.index.values, pvals, effects, percentile=0.05
                )

                res2 = pd.DataFrame()
                res2.index = pd.Index(list(res))
                res2["p_value"] = [res[k]["p_value"] for k in res]
                res2["effect_size"] = [res[k]["effect_size"] for k in res]
                res = res2

                res["cell_type_ontology_term_id"] = cell_types_o[i]
                for j, term in enumerate(self.groupby_terms):
                    res[term] = combination[j]
                res["gene_ontology_term_id"] = res.index
                res = res.reset_index(drop=True)
                res = res[res["effect_size"].notnull()]
                res = res[res["effect_size"] > 0]
                all_results.append(res)

        x_new, y_new = (e_nnz_rollup + e_sum_rollup).nonzero()

        r_x_index = pd.Series(index=x_index.values, data=x_index.index.values)
        r_y_index = pd.Series(index=y_index.values, data=y_index.index.values)

        y_r_new = r_y_index[y_new].values

        x_r_new = r_x_index[x_new].values

        new_index = pd.Index([i + (j,) for i, j in zip(x_r_new, y_r_new)])

        nnz_flat = e_nnz_rollup[x_new, y_new]
        sum_flat = e_sum_rollup[x_new, y_new]

        new_index = new_index.set_names(self.group_by_terms_with_celltype_and_gene)
        new_expression_rollup = pd.DataFrame(index=new_index)

        new_expression_rollup["nnz"] = nnz_flat
        new_expression_rollup["sum"] = sum_flat

        new_expression_rollup = new_expression_rollup.reset_index()

        markers_df = pd.concat(all_results, axis=0)
        markers_df = markers_df[markers_df["p_value"] < 1e-5]

        markers_df = markers_df.groupby(self.group_by_terms_with_celltype_and_gene).first()
        new_expression_rollup = new_expression_rollup.groupby(self.group_by_terms_with_celltype_and_gene).first()

        new_expression_rollup = new_expression_rollup.join(markers_df, on=self.group_by_terms_with_celltype_and_gene)

        new_expression_rollup = new_expression_rollup.reset_index()
        markers_df = markers_df.reset_index()

        top_per_group = markers_df.groupby(self.group_by_terms_with_celltype).apply(
            lambda x: x.nlargest(100, "effect_size")
        )

        marker_gene_groups = list(zip(*top_per_group[self.group_by_terms_with_celltype_and_gene].values.T))

        # TODO: why skip organism in jupyter?
        filt = (
            np.array(
                [
                    new_expression_rollup[term].isin(top_per_group[term].unique())
                    for term in self.group_by_terms_with_celltype_and_gene
                ]
            ).prod(0)
            > 0
        )

        new_expression_rollup = new_expression_rollup[filt]

        new_expression_rollup.index = pd.Index(
            list(zip(*new_expression_rollup[self.group_by_terms_with_celltype_and_gene].values.T))
        )

        universe_cell_counts_df = universe_cell_counts_df.groupby(self.group_by_terms_with_celltype).sum()["n_cells"]

        data = {}
        for group in marker_gene_groups:
            gene = group[-1]
            celltype = group[-2]
            group_excluding_gene = group[:-1]

            nnz = new_expression_rollup["nnz"][i]
            s = new_expression_rollup["sum"][i]
            es = new_expression_rollup["effect_size"][i]

            n_cells = universe_cell_counts_df[group_excluding_gene]

            a = data.get(celltype, [])
            entry = {
                "me": s / nnz if nnz > 0 else 0,
                "pc": nnz / n_cells,
                "marker_score": es,
                "symbol": self.gene_id_to_symbol[gene],
                "name": self.gene_id_to_name.get(gene, self.gene_id_to_symbol[gene]),
            }
            for i, term in enumerate(self.groupby_terms):
                entry[term] = group[i]
            a.append(entry)
            data[celltype] = a

            return data
