from backend.cellguide.pipeline.canonical_marker_genes.canonical_markers import CanonicalMarkerGenesCompiler
from backend.cellguide.pipeline.constants import CANONICAL_MARKER_GENES_FILENAME
from backend.cellguide.pipeline.utils import output_json
from backend.wmg.data.snapshot import WmgSnapshot, load_snapshot


def run(output_directory):
    snapshot: WmgSnapshot = load_snapshot()
    cell_counts_df = snapshot.cell_counts_cube.df[:]

    marker_gene_compiler = CanonicalMarkerGenesCompiler(cell_counts_df)
    parsed_asctb_table_entries = marker_gene_compiler.get_processed_asctb_table_entries()

    output_json(parsed_asctb_table_entries, f"{output_directory}/{CANONICAL_MARKER_GENES_FILENAME}")
