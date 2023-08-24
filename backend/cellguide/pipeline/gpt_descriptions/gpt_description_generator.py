import logging
import os
from tempfile import TemporaryDirectory

from backend.cellguide.pipeline.config import CellGuideConfig
from backend.cellguide.pipeline.providers.openai_provider import OpenAIProvider
from backend.common.providers.s3_provider import S3Provider

logger = logging.getLogger(__name__)


def generate_new_gpt_descriptions(all_cell_type_ids_in_corpus: list[str]) -> dict[str, str]:
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
    new_gpt_descriptions: dict[str, str] = {}

    for cell_type_id in all_cell_type_ids_in_corpus:
        object_key = f"gpt_descriptions/{cell_type_id}.json"
        if not S3Provider.does_object_exist(bucket_name=bucket_name, object_key=object_key):
            logger.info(f"GPT description does not exist yet for {cell_type_id}, generating...")
            description_str = OpenAIProvider.generate_chatgpt4_output(cell_type_id=cell_type_id)

            new_gpt_descriptions[cell_type_id] = description_str
            logger.info(f"Generated GPT description for {cell_type_id}: {description_str}")


def upload_descriptions_to_s3(new_gpt_descriptions: dict[str, str]) -> None:
    bucket_name = CellGuideConfig().bucket

    with TemporaryDirectory() as temp_dir:
        # write to a local temp directory gpt_descriptions/
        for cell_type_id, description_str in new_gpt_descriptions.items():
            local_file_path = os.path.join(temp_dir, f"{cell_type_id}.json")
            with open(local_file_path, "w") as file:
                file.write(description_str)

        # upload to s3
        S3Provider.sync_directory(src_dir=temp_dir, s3_uri=f"s3://{bucket_name}/gpt_descriptions/")
        logger.info(
            f"Uploaded {len(new_gpt_descriptions.keys())} GPT descriptions to s3://{bucket_name}/gpt_descriptions/"
        )
