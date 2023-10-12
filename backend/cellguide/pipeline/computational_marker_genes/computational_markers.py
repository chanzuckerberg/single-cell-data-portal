import concurrent.futures
import itertools
import logging
import warnings
from typing import Tuple

import numpy as np
import pandas as pd
from dask import compute, delayed
from dask.diagnostics import ProgressBar
from tqdm import tqdm

from backend.cellguide.pipeline.computational_marker_genes.constants import (
    MARKER_SCORE_THRESHOLD,
)
from backend.cellguide.pipeline.computational_marker_genes.types import ComputationalMarkerGenes
from backend.cellguide.pipeline.computational_marker_genes.utils import (
    bootstrap_rows_percentiles,
    calculate_cohens_d,
    calculate_pvalue_excluding_nans,
    query_gene_info_for_gene_description,
)
from backend.cellguide.pipeline.constants import CELLGUIDE_PIPELINE_NUM_CPUS
from backend.cellguide.pipeline.utils import get_gene_id_to_name_and_symbol
from backend.common.utils.rollup import (
    are_cell_types_colinear,
    get_overlapping_cell_type_descendants,
    rollup_across_cell_type_descendants,
    rollup_across_cell_type_descendants_array,
)
from backend.wmg.data.snapshot import WmgSnapshot

logger = logging.getLogger(__name__)


"""
This module contains the MarkerGenesCalculator class which is used to calculate marker genes for cell types.
The class takes in a snapshot of the data, a list of all cell type IDs in the corpus, and a list of terms to group by.
It then prepares the cell counts and gene expression dataframes, and compiles the marker genes.

The groupby terms provided by the user determine the dimensions across which the marker gene computation will be stratified. 
For instance, users can stratify marker gene calculation across just organisms, each combination of organism and tissue, 
or any arbitrary combinations of metadata dimensions.
"""


class MarkerGenesCalculator:
    def __init__(self, *, snapshot: WmgSnapshot, all_cell_type_ids_in_corpus: list[str], groupby_terms: list[str]):
        self.all_cell_type_ids_in_corpus = all_cell_type_ids_in_corpus

        gene_metadata = get_gene_id_to_name_and_symbol()
        self.gene_id_to_name = gene_metadata.gene_id_to_name
        self.gene_id_to_symbol = gene_metadata.gene_id_to_symbol

        # add the primary filter gene metadata to the gene_id_to_symbol dictionary read from Ensembl
        # this is done in case the gene metadata file is missing some gene IDs, which could occur
        # due to mismatch in gene reference versions.
        primary_filters__gene_id_to_symbol = {
            k: v
            for organism in snapshot.primary_filter_dimensions["gene_terms"]
            for i in snapshot.primary_filter_dimensions["gene_terms"][organism]
            for k, v in i.items()
        }
        self.gene_id_to_symbol.update(primary_filters__gene_id_to_symbol)

        # Groupby variables used to group the data various operations
        # cell types are removed as they are treated differently
        # all other metadata (like tissue and organism) are dimensions across
        # which marker gene computation for cell types will be stratified
        if "cell_type_ontology_term_id" in groupby_terms:
            groupby_terms.remove("cell_type_ontology_term_id")

        self.groupby_terms = groupby_terms
        self.groupby_terms_with_celltype = groupby_terms + ["cell_type_ontology_term_id"]
        self.groupby_terms_with_celltype_and_gene = self.groupby_terms_with_celltype + ["gene_ontology_term_id"]

        # load the cell counts and expression summary cubes fully in memory
        cell_counts_df = snapshot.cell_counts_cube.df[:]
        expressions_df = snapshot.expression_summary_default_cube.df[:]

        # prep the cell counts and expressions dataframes
        (
            self.cell_counts_df,
            self.cell_counts_df_orig,
            self.expressions_df,
        ) = self._prepare_cell_counts_and_gene_expression_dfs(cell_counts_df, expressions_df)

    def _get_gene_symbol_from_id(self, gene_id: str) -> str:
        return self.gene_id_to_symbol.get(gene_id, gene_id)

    def _get_gene_name_from_id(self, gene_id: str) -> str:
        # If gene_id is already in the dictionary, return the name
        if gene_id in self.gene_id_to_name:
            return self.gene_id_to_name[gene_id]

        # If gene_id_to_name_memory doesn't exist, initialize it
        if not hasattr(self, "gene_id_to_name_memory"):
            self.gene_id_to_name_memory = {}

        # If gene_id is not in the memory, query the gene info and store it
        if gene_id not in self.gene_id_to_name_memory:
            gene_name = query_gene_info_for_gene_description(gene_id)
            # If gene_name is in the gene_id_to_symbol dictionary, replace it
            if gene_name in self.gene_id_to_symbol:
                gene_name = self.gene_id_to_symbol[gene_name]
            self.gene_id_to_name_memory[gene_id] = gene_name

        return self.gene_id_to_name_memory[gene_id]

    def _prepare_cell_counts_and_gene_expression_dfs(
        self, cell_counts_df: pd.DataFrame, expressions_df: pd.DataFrame
    ) -> Tuple[pd.DataFrame, pd.DataFrame]:
        """
        This method prepares the cell counts and gene expression dataframes for further processing.
        It groups the cell counts and expressions dataframes by the terms specified in the groupby_terms_with_celltype attribute.

        Arguments
        ---------
        cell_counts_df - pd.DataFrame
            The dataframe containing cell counts data.
        expressions_df - pd.DataFrame
            The dataframe containing gene expression data.

        Returns
        -------
        Tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]
            A tuple containing the prepared cell counts, unrolled cell counts, and gene expression dataframes.
        """
        logger.info("Preparing cell counts and gene expression dataframes")

        # group by cell counts
        cell_counts_df = cell_counts_df.groupby(self.groupby_terms_with_celltype).sum(numeric_only=True)
        cell_counts_df = cell_counts_df[
            cell_counts_df.index.get_level_values("cell_type_ontology_term_id").isin(self.all_cell_type_ids_in_corpus)
        ]

        # group by gene expressions
        expressions_df = expressions_df.groupby(self.groupby_terms_with_celltype_and_gene).sum(numeric_only=True)
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

        # create an unrolled copy of the cell counts dataframe
        cell_counts_df_orig = universe_cell_counts_df.copy()
        cell_counts_df_orig["n_cells"] = 0
        cell_counts_df_orig["n_cells"][cell_counts_df.index] = cell_counts_df["n_cells"]
        return universe_cell_counts_df, cell_counts_df_orig, expressions_df

    def _process_cell_type__parallel(
        self,
        *,
        i: int,
        cell_types_o: np.ndarray,
        e_sum_o: np.ndarray,
        e_sqsum_o: np.ndarray,
        n_cells_o: np.ndarray,
        e_sum_o_orig: np.ndarray,
        e_sqsum_o_orig: np.ndarray,
        n_cells_o_orig: np.ndarray,
        filter_genes: np.ndarray,
    ) -> tuple[str, np.ndarray, np.ndarray, np.ndarray]:
        """
        This function is used in dask scheduler to process each cell type in parallel.

        Parameters
        ----------
        i : int
            The index of the current cell type being processed.
        cell_types_o : np.ndarray
            Array of cell types.
        e_sum_o : np.ndarray
            Array of gene expression sums for each cell type.
        e_sqsum_o : np.ndarray
            Array of squared gene expression sums for each cell type.
        n_cells_o : np.ndarray
            Array of cell counts for each cell type.
        e_sum_o_orig : np.ndarray
            Array of raw gene expression sums for each cell type.
        e_sqsum_o_orig : np.ndarray
            Array of raw squared gene expression sums for each cell type.
        n_cells_o_orig : np.ndarray
            Array of cell counts for each cell type.
        filter_genes: np.ndarray
            A boolean array indicating which genes to filter out.

        Returns
        -------
        tuple[np.ndarray, np.ndarray]
            A tuple containing the filtered effect sizes and the column indices of the remaining genes.
        """
        e_sum_o = e_sum_o.copy()
        e_sqsum_o = e_sqsum_o.copy()
        n_cells_o = n_cells_o.copy()
        cell_type_target = cell_types_o[i]

        indexer = pd.Series(index=cell_types_o, data=np.arange(cell_types_o.size))

        for j, cell_type in enumerate(cell_types_o):
            if cell_type_target == cell_type or are_cell_types_colinear(cell_type, cell_type_target):
                continue

            overlapping_descendants = get_overlapping_cell_type_descendants(cell_type, cell_type_target)
            overlapping_descendants = list(set(overlapping_descendants).intersection(cell_types_o))
            if len(overlapping_descendants) > 0:
                e_sum_o[j] -= e_sum_o_orig[indexer[overlapping_descendants].values].sum(0)
                e_sqsum_o[j] -= e_sqsum_o_orig[indexer[overlapping_descendants].values].sum(0)
                n_cells_o[j] -= n_cells_o_orig[indexer[overlapping_descendants].values].sum(0)

        # get the expressions and cell counts corresponding to the cell type
        sum1 = e_sum_o[i][None, :]
        sumsq1 = e_sqsum_o[i][None, :]
        n1 = n_cells_o[i][None, :]

        # calculate cohens d effect size against all other rows
        effects = calculate_cohens_d(sum1=sum1, sumsq1=sumsq1, n1=n1, sum2=e_sum_o, sumsq2=e_sqsum_o, n2=n_cells_o)

        # zero out nans
        effects[np.isnan(effects)] = 0

        # get valid genes
        unique_cols = np.where(~filter_genes)[0]

        # filter out rows that are colinear with the cell type target and filter out invalid genes
        is_colinear = np.array([are_cell_types_colinear(cell_type, cell_type_target) for cell_type in cell_types_o])
        effects_sub = effects[~is_colinear][:, unique_cols]

        return effects_sub, unique_cols

    def get_computational_marker_genes(
        self,
        num_marker_genes=100,
        minimum_nnz=25,
        num_replicates=100,
        percentile=10,
    ) -> dict[str, list[ComputationalMarkerGenes]]:
        """
        Calculate the computational marker genes for each cell type.

        Arguments
        ---------
        num_marker_genes : int, optional
            The maximum number of marker genes to return, by default 100
        minimum_nnz : int, optional
            The minimum number of non-zero expressions for a gene to be considered, by default 25
        num_replicates : int, optional
            The number of bootstrap replicates to generate for calcualting percentiles, by default 100
        percentile : int, optional
            The percentile for the bootstrap, by default 10

        Returns
        -------
        dict[str, list[ComputationalMarkerGenes]]
            A dictionary where the keys are cell types and the values are lists of ComputationalMarkerGenes objects.
        """
        logger.info("Getting computational marker genes")
        cell_counts_df = self.cell_counts_df
        cell_counts_df_orig = self.cell_counts_df_orig

        # the metadata groups (incl cell type) will be treated as row coordinates
        groupby_coords = list(zip(*self.expressions_df[self.groupby_terms_with_celltype].values.T))
        groupby_coords_unique = sorted(set(groupby_coords))
        groupby_index = pd.Series(index=pd.Index(groupby_coords_unique), data=np.arange(len(groupby_coords_unique)))

        # the genes will be treated as column coordinates
        gene_coords = list(self.expressions_df["gene_ontology_term_id"])
        gene_coords_unique = sorted(set(gene_coords))
        gene_index = pd.Series(index=pd.Index(gene_coords_unique), data=np.arange(len(gene_coords_unique)))

        # instantiate empty numpy arrays with number of rows equal to the number of groups
        # and the number of cols equal to the number of genes
        e_nnz = np.zeros((len(groupby_coords_unique), len(gene_coords_unique)))
        e_sum = np.zeros((len(groupby_coords_unique), len(gene_coords_unique)))
        e_sqsum = np.zeros((len(groupby_coords_unique), len(gene_coords_unique)))

        logger.info("Populating arrays with numeric data from the expressions dataframe")
        # populate the arrays with the numeric data from the expressions dataframe
        e_nnz[groupby_index[groupby_coords].values, gene_index[gene_coords].values] = self.expressions_df["nnz"].values
        e_sum[groupby_index[groupby_coords].values, gene_index[gene_coords].values] = self.expressions_df["sum"].values
        e_sqsum[groupby_index[groupby_coords].values, gene_index[gene_coords].values] = self.expressions_df[
            "sqsum"
        ].values

        # get all available combinations from the augmented cell counts dataframe
        available_combinations = set(cell_counts_df.index.values)
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
        groupby_term_to_unique_values = [sorted(set(i)) for i in groupby_term_to_values]
        # cell types will always be the last level in the multiindex by convention
        cell_types = groupby_index.index.get_level_values(-1)

        # instantiate empty rolled up expression arrays
        e_nnz_rollup = np.zeros_like(e_nnz)
        e_sum_rollup = np.zeros_like(e_sum)
        e_sqsum_rollup = np.zeros_like(e_sqsum)

        # get the number of cells for each group and tile it into an array that is
        # the same shape as the expression arrays above
        n_cells = cell_counts_df["n_cells"][groupby_index.index]
        n_cells = np.tile(n_cells.values[:, None], (1, e_nnz_rollup.shape[1]))

        n_cells_orig = cell_counts_df_orig["n_cells"][groupby_index.index]
        n_cells_orig = np.tile(n_cells_orig.values[:, None], (1, e_nnz_rollup.shape[1]))

        all_results = []
        # for example, if self.groupby_terms contains organism and tissue, then this loop
        # iterates through each organism and tissue combination.
        logger.info(f"Iterating through all combinations of groupby dimensions {self.groupby_terms}")
        for combination in itertools.product(*groupby_term_to_unique_values):
            # get the rows corresponding to groups that match the current "combination"
            filt = groupby_term_to_values[0] == combination[0]
            for _ in range(1, len(combination)):
                filt = np.logical_and(filt, groupby_term_to_values[1] == combination[1])

            if filt.sum() == 0:
                continue

            # filter down to the specified subset
            cell_types_o = cell_types[filt]
            e_nnz_o = e_nnz[filt]
            e_sum_o = e_sum[filt]
            e_sqsum_o = e_sqsum[filt]
            n_cells_o = n_cells[filt]
            n_cells_orig_o = n_cells_orig[filt]

            # roll up the arrays across the rows given the cell types corresponding to each row
            e_nnz_o_rollup = rollup_across_cell_type_descendants_array(e_nnz_o, cell_types_o)
            e_sum_o_rollup = rollup_across_cell_type_descendants_array(e_sum_o, cell_types_o)
            e_sqsum_o_rollup = rollup_across_cell_type_descendants_array(e_sqsum_o, cell_types_o)

            # populate the original rollup arrays with the rolled up values for this combination
            e_nnz_rollup[filt] = e_nnz_o_rollup
            e_sum_rollup[filt] = e_sum_o_rollup
            e_sqsum_rollup[filt] = e_sqsum_o_rollup

            delayed_results = [
                delayed(self._process_cell_type__parallel)(
                    i=i,
                    cell_types_o=cell_types_o,
                    e_sum_o=e_sum_o_rollup,
                    e_sqsum_o=e_sqsum_o_rollup,
                    n_cells_o=n_cells_o,
                    e_sum_o_orig=e_sum_o,
                    e_sqsum_o_orig=e_sqsum_o,
                    n_cells_o_orig=n_cells_orig_o,
                    filter_genes=e_nnz_o[i] < minimum_nnz,
                )
                for i in range(len(cell_types_o))
            ]
            logger.info(
                f"Getting marker genes for {len(cell_types_o)} cell types in combination {combination} using {CELLGUIDE_PIPELINE_NUM_CPUS} CPUs..."
            )
            with warnings.catch_warnings():
                warnings.simplefilter("ignore")

                with ProgressBar():
                    results = compute(*delayed_results, num_workers=CELLGUIDE_PIPELINE_NUM_CPUS)

                effect_sizes = []
                for iteration in tqdm(
                    range(len(results)),
                    desc="Bootstrapping marker scores",
                ):
                    effect_sizes_chunk, col_idx = results[iteration]
                    effects = np.full(gene_index.size, fill_value=np.nan)
                    bootstrapped_percentiles = np.full((num_replicates, gene_index.size), fill_value=np.nan)

                    if effect_sizes_chunk.shape[0] > 0:
                        bootstrapped_percentiles[:, col_idx] = bootstrap_rows_percentiles(
                            effect_sizes_chunk,
                            num_replicates=num_replicates,
                            num_samples=effect_sizes_chunk.shape[0],
                            percentile=percentile,
                        )

                        effects[col_idx] = bootstrapped_percentiles[:, col_idx].mean(0)

                    effect_sizes.append(effects)
                effect_sizes = np.vstack(effect_sizes)

                for iteration in tqdm(range(len(effect_sizes)), desc="Processing cell type marker genes"):
                    cell_type = cell_types_o[iteration]

                    effect_size = effect_sizes[iteration]
                    not_overlapping = np.array([not are_cell_types_colinear(ct, cell_type) for ct in cell_types_o])

                    specificity = calculate_pvalue_excluding_nans(effect_size, effect_sizes[not_overlapping])
                    ranked_genes_df = pd.DataFrame()
                    ranked_genes_df["gene_ontology_term_id"] = gene_index.index.values
                    ranked_genes_df["specificity"] = specificity
                    ranked_genes_df["effect_size"] = effect_size
                    ranked_genes_df["cell_type_ontology_term_id"] = cell_type

                    for j, term in enumerate(self.groupby_terms):
                        ranked_genes_df[term] = combination[j]

                    ranked_genes_df = ranked_genes_df[ranked_genes_df["effect_size"].notnull()]
                    ranked_genes_df = ranked_genes_df[ranked_genes_df["effect_size"] > MARKER_SCORE_THRESHOLD]
                    all_results.append(ranked_genes_df)

        # concatenate all the results into one marker gene dataframe
        markers_df = pd.concat(all_results, axis=0)

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

        # join the rolled up expression dataframe to the marker genes dataframe along the index
        # to combine the expression information with the marker gene scores
        new_expression_rollup = new_expression_rollup.join(markers_df, on=self.groupby_terms_with_celltype_and_gene)

        # reset the index to convert MultiIndex back into columns
        markers_df = markers_df.reset_index()

        # get the top `num_marker_genes` genes per metadata group
        top_per_group = (
            markers_df.groupby(self.groupby_terms_with_celltype)
            .apply(lambda x: x.nlargest(num_marker_genes, "effect_size"))
            .reset_index(drop=True)
        )
        # get the marker gene groups
        marker_gene_groups = list(zip(*top_per_group[self.groupby_terms_with_celltype_and_gene].values.T))

        # convert columns to MultiIndex
        top_per_group.set_index(self.groupby_terms_with_celltype_and_gene, inplace=True)
        # filter the rollup expression df down to the rows that are among the top marker gene groups
        filt = new_expression_rollup.index.isin(top_per_group.index)
        new_expression_rollup = new_expression_rollup[filt].reset_index()
        # set the groupby+gene columns as a multi-index
        new_expression_rollup.set_index(self.groupby_terms_with_celltype_and_gene, inplace=True)

        # filter the expressions down to the top marker gene groups
        new_expression_rollup = new_expression_rollup.loc[marker_gene_groups]

        # join the cell counts dataframe to the expressions dataframe along its index
        new_expression_rollup = new_expression_rollup.join(cell_counts_df, on=cell_counts_df.index.names)

        # calculate the mean expression and percent cells for each marker gene group
        new_expression_rollup["me"] = new_expression_rollup["sum"] / new_expression_rollup["nnz"]
        new_expression_rollup["pc"] = new_expression_rollup["nnz"] / new_expression_rollup["n_cells"]
        # ensure that the percent cells is between 0 and 1
        assert new_expression_rollup["pc"].max() <= 1.0

        # get all the records from the expressions dataframe
        records = new_expression_rollup.reset_index().to_dict(orient="records")

        # get gene names from IDs
        gene_names_to_ids = {}
        with concurrent.futures.ThreadPoolExecutor() as executor:
            futures = {
                executor.submit(self._get_gene_name_from_id, datum["gene_ontology_term_id"]): datum[
                    "gene_ontology_term_id"
                ]
                for datum in records
            }

        for future in concurrent.futures.as_completed(futures):
            gene_name = future.result()
            gene_id = futures[future]
            gene_names_to_ids[gene_id] = gene_name

        # format the records into a dictionary mapping cell type IDs to lists of marker genes
        formatted_data = {}
        for datum in records:
            marker_gene_list = formatted_data.get(datum["cell_type_ontology_term_id"], [])
            entry = {
                "me": datum["me"],
                "pc": datum["pc"],
                "marker_score": datum["effect_size"],
                "specificity": datum["specificity"],
                "gene_ontology_term_id": datum["gene_ontology_term_id"],
                "symbol": self._get_gene_symbol_from_id(datum["gene_ontology_term_id"]),
                "name": gene_names_to_ids[datum["gene_ontology_term_id"]],
            }
            entry["groupby_dims"] = {term: datum[term] for term in self.groupby_terms}

            marker_gene_list.append(ComputationalMarkerGenes(**entry))
            formatted_data[datum["cell_type_ontology_term_id"]] = marker_gene_list

        return formatted_data
