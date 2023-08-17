import itertools
import logging
import pathlib

import numpy as np
import pandas as pd

from backend.cellguide.pipeline.computational_marker_genes.utils import (
    post_process_stats,
    run_ttest,
)
from backend.common.utils.rollup import rollup_across_cell_type_descendants, rollup_across_cell_type_descendants_array
from backend.wmg.data.snapshot import WmgSnapshot

logger = logging.getLogger(__name__)

ENSEMBL_GENE_ID_TO_DESCRIPTION_FILENAME = "ensembl_gene_ids_to_descriptions.tsv.gz"

"""
This module contains the MarkerGenesCalculator class which is used to calculate marker genes for cell types.
The class takes in a snapshot of the data, a list of all cell type IDs in the corpus, and a list of terms to group by.
It then prepares the cell counts and gene expression dataframes, and compiles the marker genes.

The groupby terms provided by the user determine the dimensions across which the marker gene computation will be stratified. 
For instance, users can stratify marker gene calculation across just organisms, each combination of organism and tissue, 
or any arbitrary combinations of metadata dimensions.
"""


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

        # Groupby variables used to group the data various operations
        # cell types are removed as they are treated differently
        # all other metadata (like tissue and organism) are dimensions across
        # which marker gene computation for cell types will be stratified
        if "cell_type_ontology_term_id" in groupby_terms:
            groupby_terms.remove("cell_type_ontology_term_id")
        else:
            raise ValueError("cell_type_ontology_term_id must be one of the groupby terms")

        self.groupby_terms = groupby_terms
        self.groupby_terms_with_celltype = groupby_terms + ["cell_type_ontology_term_id"]
        self.groupby_terms_with_celltype_and_gene = self.groupby_terms_with_celltype + ["cell_type_ontology_term_id"]

        # load the cell counts and expression summary cubes fully in memory
        cell_counts_df = snapshot.cell_counts_cube.df[:]
        expressions_df = snapshot.expression_summary_fmg_cube.df[:]

        # prep the cell counts and expressions dataframes
        (self.universe_cell_counts_df, self.expressions_df) = self._prepare_cell_counts_and_gene_expression_dfs(
            cell_counts_df, expressions_df
        )

    def _get_symbol_to_gene_descriptions_file_path(self) -> pathlib.Path:
        """
        TODO: Refactor this as a utility function as it will be shared with canonical marker genes

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
        cell_counts_df = cell_counts_df.groupby(self.groupby_terms_with_celltype).sum(numeric_only=True)
        cell_counts_df = cell_counts_df[
            cell_counts_df.index.get_level_values("cell_type_ontology_term_id").isin(self.all_cell_type_ids_in_corpus)
        ]

        # group by gene expressions
        expressions_df = expressions_df.groupby(self.groupby_terms_with_celltype + ["gene_ontology_term_id"]).sum(
            numeric_only=True
        )
        expressions_df = expressions_df.reset_index()
        expressions_df = expressions_df[
            expressions_df["cell_type_ontology_term_id"].isin(self.all_cell_type_ids_in_corpus)
        ]

        # get the cartesian product of the groupby metadata and all the cell types in the corpus
        index = pd.Index(
            list(
                itertools.product(
                    *[cell_counts_df.index.get_level_values(term).unique() for term in self.groupby_terms]
                    + [self.all_cell_type_ids_in_corpus]
                )
            )
        )
        index = index.set_names(self.groupby_terms_with_celltype)
        # instantiate an empty dataframe with the groups from the cartesian product as the index
        universe_cell_counts_df = pd.DataFrame(index=index)
        universe_cell_counts_df["n_cells"] = 0
        # populate the dataframe with the number of cells from the groups in the cell counts df
        universe_cell_counts_df["n_cells"][cell_counts_df.index] = cell_counts_df["n_cells"]
        # roll up the cell counts
        universe_cell_counts_df = rollup_across_cell_type_descendants(universe_cell_counts_df.reset_index())
        # remove groups that still have 0 cells after the rollup operation
        universe_cell_counts_df = universe_cell_counts_df[universe_cell_counts_df["n_cells"] > 0]
        # remake the multi-index
        universe_cell_counts_df = universe_cell_counts_df.groupby(self.groupby_terms_with_celltype).sum()
        return universe_cell_counts_df, expressions_df

    def get_computational_marker_genes(self):
        universe_cell_counts_df = self.universe_cell_counts_df

        # the metadata groups (incl cell type) will be treated as row coordinates
        groupby_coords = list(zip(*self.expressions_df[self.groupby_terms_with_celltype].values.T))
        groupby_coords_unique = list(set(groupby_coords))
        groupby_index = pd.Series(index=pd.Index(groupby_coords_unique), data=np.arange(len(groupby_coords_unique)))

        # the genes will be treated as column coordinates
        gene_coords = list(self.expressions_df["gene_ontology_term_id"])
        gene_coords_unique = list(set(gene_coords))
        gene_index = pd.Series(index=pd.Index(gene_coords_unique), data=np.arange(len(gene_coords_unique)))

        # instantiate empty numpy arrays with number of rows equal to the number of groups
        # and the number of cols equal to the number of genes
        e_nnz = np.zeros((len(groupby_coords_unique), len(gene_coords_unique)))
        e_sum = np.zeros((len(groupby_coords_unique), len(gene_coords_unique)))
        e_sqsum = np.zeros((len(groupby_coords_unique), len(gene_coords_unique)))

        # populate the arrays with the numeric data from the expressions dataframe
        e_nnz[groupby_index[groupby_coords].values, gene_index[gene_coords].values] = self.expressions_df["nnz"]
        e_sum[groupby_index[groupby_coords].values, gene_index[gene_coords].values] = self.expressions_df["sum"]
        e_sqsum[groupby_index[groupby_coords].values, gene_index[gene_coords].values] = self.expressions_df["sqsum"]

        # get all available combinations from the augmented cell counts dataframe
        available_combinations = set(universe_cell_counts_df.index.values)
        # get all the combinations that are missing from the expressions dataframe
        missing_combinations = available_combinations.difference(groupby_coords_unique)

        # augment the row coordinates with the missing groups
        groupby_coords_unique = groupby_coords_unique + list(missing_combinations)
        groupby_index = pd.Series(index=pd.Index(groupby_coords_unique), data=np.arange(len(groupby_coords_unique)))

        # augment the expression arrays with empty values for the missing groups
        e_nnz = np.vstack((e_nnz, np.zeros((len(missing_combinations), e_nnz.shape[1]))))
        e_sum = np.vstack((e_sum, np.zeros((len(missing_combinations), e_sum.shape[1]))))
        e_sqsum = np.vstack((e_sqsum, np.zeros((len(missing_combinations), e_sqsum.shape[1]))))

        # for each groupby term, get the corresponding values from the multiindex
        groupby_term_to_values = [groupby_index.index.get_level_values(i) for i in range(len(self.groupby_terms))]
        # get the unique values for each groupby term
        groupby_term_to_unique_values = [list(set(i)) for i in groupby_term_to_values]
        # cell types will always be the last level in the multiindex by convention
        cell_types = groupby_index.index.get_level_values(-1)

        # instantiate empty rolled up expression arrays
        e_nnz_rollup = np.zeros_like(e_nnz)
        e_sum_rollup = np.zeros_like(e_sum)
        e_sqsum_rollup = np.zeros_like(e_sqsum)

        # get the number of cells for each group and tile it into an array that is
        # the same shape as the expression arrays above
        n_cells = universe_cell_counts_df["n_cells"][groupby_index.index]
        n_cells = np.tile(n_cells.values[:, None], (1, e_nnz_rollup.shape[1]))

        all_results = []
        # for example, if self.groupby_terms contains organism and tissue, then this loop
        # iterates through each organism and tissue combination.
        logger.info(f"Iterating through all combinations of groupby dimensions {self.groupby_terms}")
        for combination in itertools.product(*groupby_term_to_unique_values):
            logger.info(f"Getting marker genes for cell types in combination {combination}")
            # get the rows corresponding to groups that match the current "combination"
            filt = groupby_term_to_values[0] == combination[0]
            for _i in range(1, len(combination)):
                filt = np.logical_and(filt, groupby_term_to_values[1] == combination[1])

            if filt.sum() == 0:
                continue

            # filter down to the specified subset
            cell_types_o = cell_types[filt]
            e_nnz_o = e_nnz[filt]
            e_sum_o = e_sum[filt]
            e_sqsum_o = e_sqsum[filt]
            n_cells_o = n_cells[filt]

            # roll up the arrays across the rows given the cell types corresponding to each row
            e_nnz_o = rollup_across_cell_type_descendants_array(e_nnz_o, cell_types_o)
            e_sum_o = rollup_across_cell_type_descendants_array(e_sum_o, cell_types_o)
            e_sqsum_o = rollup_across_cell_type_descendants_array(e_sqsum_o, cell_types_o)

            # populate the original rollup arrays with the rolled up values for this combination
            e_nnz_rollup[filt] = e_nnz_o
            e_sum_rollup[filt] = e_sum_o
            e_sqsum_rollup[filt] = e_sqsum_o

            # for each cell type in the current subset
            for i in range(len(cell_types_o)):
                logger.info(f"Getting marker genes for cell type {cell_types_o[i]} ({i+1} / {len(cell_types_o)})...")
                # get the expressions and cell counts corresponding to the cell type
                sum1 = e_sum_o[i][None, :].copy()
                sumsq1 = e_sqsum_o[i][None, :].copy()
                n1 = n_cells_o[i][None, :].copy()

                # run the t-test elementwise against all other rows
                pvals, effects = run_ttest(sum1, sumsq1, n1, e_sum_o, e_sqsum_o, n_cells_o)

                # set the results for the current cell type's row to nan
                # (cell type should not be compared to itself)
                pvals[i] = np.nan
                effects[i] = np.nan

                # process the stats to get the ranked list of differentially expressed genes
                ranked_genes = post_process_stats(
                    cell_types_o[i], cell_types_o, gene_index.index.values, pvals, effects, percentile=0.05
                )

                # create a dataframe containing the ranked genes, p-values, effect sizes, current combination
                # and current cell type
                ranked_genes_df = pd.DataFrame()
                ranked_genes_df.index = pd.Index(list(ranked_genes))
                ranked_genes_df["p_value"] = [ranked_genes[k]["p_value"] for k in ranked_genes]
                ranked_genes_df["effect_size"] = [ranked_genes[k]["effect_size"] for k in ranked_genes]

                ranked_genes_df["cell_type_ontology_term_id"] = cell_types_o[i]
                for j, term in enumerate(self.groupby_terms):
                    ranked_genes_df[term] = combination[j]
                ranked_genes_df["gene_ontology_term_id"] = ranked_genes_df.index
                ranked_genes_df = ranked_genes_df.reset_index(drop=True)
                ranked_genes_df = ranked_genes_df[ranked_genes_df["effect_size"].notnull()]
                ranked_genes_df = ranked_genes_df[ranked_genes_df["effect_size"] > 0]
                all_results.append(ranked_genes_df)

        # concatenate all the results into one marker gene dataframe
        markers_df = pd.concat(all_results, axis=0)
        # filter out rows with p-values >= 1e-5 (arbitrary heuristic)
        markers_df = markers_df[markers_df["p_value"] < 1e-5]
        # use the groupby operation to convert the groupby columns into a MultiIndex
        markers_df = markers_df.groupby(self.groupby_terms_with_celltype_and_gene).first()

        # get the row and col indices corresponding to nonzero expression values
        groupby_i_coords_new, gene_i_coords_new = (e_nnz_rollup + e_sum_rollup).nonzero()

        # build the reverse indexers (going from integer indices to the original metadata values)
        reverse_groupby_index = pd.Series(index=groupby_index.values, data=groupby_index.index.values)
        reverse_gene_index = pd.Series(index=gene_index.values, data=gene_index.index.values)

        # get the metadata groups and genes corresponding to each integer index
        reverse_gene_coords_new = reverse_gene_index[gene_i_coords_new].values
        reverse_groupby_coords_new = reverse_groupby_index[groupby_i_coords_new].values

        # combine the metadata groups and genes into one MultiIndex
        new_index = pd.Index([i + (j,) for i, j in zip(reverse_groupby_coords_new, reverse_gene_coords_new)])
        new_index = new_index.set_names(self.groupby_terms_with_celltype_and_gene)

        # instantiate the rolled up expression dataframe
        new_expression_rollup = pd.DataFrame(index=new_index)

        # extract the nonzero values from the expression data
        nnz_flat = e_nnz_rollup[groupby_i_coords_new, gene_i_coords_new]
        sum_flat = e_sum_rollup[groupby_i_coords_new, gene_i_coords_new]

        # populate the rolled up expression dataframe
        new_expression_rollup["nnz"] = nnz_flat
        new_expression_rollup["sum"] = sum_flat
        self.markers_df = markers_df
        self.groupby_i_coords_new = groupby_i_coords_new
        self.gene_i_coords_new = gene_i_coords_new
        self.reverse_groupby_index = reverse_groupby_index
        self.reverse_gene_index = reverse_gene_index
        self.reverse_gene_coords_new = reverse_gene_coords_new
        self.reverse_groupby_coords_new = reverse_groupby_coords_new
        self.new_index = new_index
        self.new_expression_rollup = new_expression_rollup
        self.nnz_flat = nnz_flat
        self.sum_flat = sum_flat
        # Add the numpy arrays
        self.e_nnz = e_nnz
        self.e_sum = e_sum
        self.e_sqsum = e_sqsum
        self.e_nnz_rollup = e_nnz_rollup
        self.e_sum_rollup = e_sum_rollup
        self.e_sqsum_rollup = e_sqsum_rollup

        return self

        # join the rolled up expression dataframe to the marker genes dataframe along the index
        # to combine the expression information with the marker gene scores
        new_expression_rollup = new_expression_rollup.join(markers_df, on=self.groupby_terms_with_celltype_and_gene)

        # reset the index to convert MultiIndex back into columns
        markers_df = markers_df.reset_index()

        # get the top 100 genes per metadata group
        top_per_group = markers_df.groupby(self.groupby_terms_with_celltype).apply(
            lambda x: x.nlargest(100, "effect_size")
        )
        # get the marker gene groups
        marker_gene_groups = list(zip(*top_per_group[self.groupby_terms_with_celltype_and_gene].values.T))

        # filter the rollup expression df down to the rows that are among the top marker gene groups
        filt = (
            np.array(
                [
                    new_expression_rollup.index.get_level_values(term).isin(top_per_group[term].unique())
                    for term in self.groupby_terms_with_celltype_and_gene
                ]
            ).prod(0)
            > 0
        )
        new_expression_rollup = new_expression_rollup[filt]
        # copy the groupby+gene columns into a multi-index
        new_expression_rollup.index = pd.Index(
            list(zip(*new_expression_rollup[self.groupby_terms_with_celltype_and_gene].values.T))
        )

        universe_cell_counts_df = universe_cell_counts_df.groupby(self.groupby_terms_with_celltype).sum()["n_cells"]

        data = {}
        for group in marker_gene_groups:
            # each group is by convention going to be (...user-specified groupby dimensions..., cell type, gene)
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
