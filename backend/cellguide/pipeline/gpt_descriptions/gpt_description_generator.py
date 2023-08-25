import logging
import time

from backend.cellguide.pipeline.config import CellGuideConfig
from backend.cellguide.pipeline.gpt_descriptions.constants import (
    GPT_CELLTYPE_DESCRIPTION_SYSTEM_ROLE,
    GPT_CELLTYPE_DESCRIPTION_USER_ROLE,
)
from backend.cellguide.pipeline.providers.openai_provider import OpenAIProvider
from backend.common.providers.s3_provider import S3Provider

logger = logging.getLogger(__name__)


def generate_new_gpt_descriptions(all_cell_type_ids_to_labels_in_corpus: dict[str, str]) -> dict[str, str]:
    """
    Generates GPT descriptions for all cell type ID's in the corpus if they haven't
    already been generated. Cell type ids should be stored in the cellguide-data-public-<env>
    bucket in the gpt_descriptions folder. For each cell type ID, we'll have a JSON file
    that contains the GPT description

    📁 gpt_descriptions
    ﹂CL_0000000.json
    ﹂CL_0000001.json
    ﹂...
    """
    bucket_name = CellGuideConfig().bucket
    s3_provider = S3Provider()
    openai_provider = OpenAIProvider()

    new_gpt_descriptions: dict[str, str] = {}
    for cell_type_id in all_cell_type_ids_to_labels_in_corpus:
        cell_type_id_filename = cell_type_id.replace(":", "_")
        object_key = f"gpt_descriptions/{cell_type_id_filename}.json"
        if not s3_provider.does_object_exist(bucket_name=bucket_name, object_key=object_key):
            cell_type_label = all_cell_type_ids_to_labels_in_corpus[cell_type_id]
            logger.info(f"GPT description does not exist yet for {cell_type_id} ({cell_type_label}), generating...")
            description_str = openai_provider.generate_gpt_output(
                user_role=GPT_CELLTYPE_DESCRIPTION_USER_ROLE(cell_type_label),
                system_role=GPT_CELLTYPE_DESCRIPTION_SYSTEM_ROLE,
            )
            new_gpt_descriptions[cell_type_id] = description_str
            logger.info(f"Generated GPT description for {cell_type_id}: {description_str}")

            # extra safety to avoid rate limiting
            time.sleep(1)

    return new_gpt_descriptions
