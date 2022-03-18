from typing import Iterable

from pandas import DataFrame


def apply_cell_type_ordering(cell_type_ontology_term_ids: Iterable, cell_type_ordering: DataFrame) -> DataFrame:
    unordered = DataFrame(cell_type_ontology_term_ids, columns=['cell_type_ontology_term_id'])
    expr_summary_with_ordering = unordered.join(cell_type_ordering,
                                                on=['tissue_ontology_term_id',
                                                    'cell_type_ontology_term_id'])
    return expr_summary_with_ordering. \
        sort_values(by=['tissue_ontology_term_id', 'cell_type_ontology_term_order'],
                    ignore_index=True).drop(columns=['cell_type_ontology_term_order'])
