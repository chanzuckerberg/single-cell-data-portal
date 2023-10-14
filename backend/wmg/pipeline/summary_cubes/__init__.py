import time

from backend.wmg.pipeline.summary_cubes.builder import SummaryCubesBuilder

ORGANISM_INFO = [
    {"label": "homo_sapiens", "id": "NCBITaxon:9606"},
    {"label": "mus_musculus", "id": "NCBITaxon:10090"},
]


def run():
    snapshot_id = str(int(time.time()))
    for organism in ORGANISM_INFO:
        builder = SummaryCubesBuilder(corpus_path=snapshot_id, organismInfo=organism)
        builder.run_pipeline()


# from backend.common.utils.exceptions import CubeValidationException
# from backend.common.utils.result_notification import (
#     format_failed_batch_issue_slack_alert,
#     gen_wmg_pipeline_failure_message,
#     notify_slack,
# )
# from backend.wmg.data.validation.validation import Validation
# from backend.wmg.pipeline.summary_cubes.cell_count import create_cell_count_cube
# from backend.wmg.pipeline.summary_cubes.expression_summary.job import create_expression_summary_cube
# from backend.wmg.pipeline.summary_cubes.expression_summary_fmg.job import create_expression_summary_fmg_cube
# from backend.wmg.pipeline.summary_cubes.marker_genes import create_marker_genes_cube


# def run(corpus_path: str, validate_cube: bool) -> dict:
#     """
#     Build expression summary cube and cell count cube based
#     on cell data stored in integrated corpus
#     validate expression summary cube based on biological expectations
#     if indicated by param
#     """
#     create_expression_summary_cube(corpus_path)
#     create_expression_summary_cube(corpus_path, default=True)
#     create_expression_summary_fmg_cube(corpus_path)
#     create_cell_count_cube(corpus_path)
#     create_marker_genes_cube(corpus_path)

#     if validate_cube:
#         try:
#             is_valid = Validation(corpus_path).validate_cube()
#         except Exception as e:
#             raise CubeValidationException from e
#         if is_valid is False:
#             pipeline_failure_message = gen_wmg_pipeline_failure_message(
#                 "Issue with cube validation, see logs for more detail"
#             )
#             data = format_failed_batch_issue_slack_alert(pipeline_failure_message)
#             notify_slack(data)
#             raise CubeValidationException
