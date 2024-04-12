import logging
import time

from backend.geneguide.pipeline.config import GeneGuideConfig
from backend.geneguide.pipeline.constants import GPT_OUTPUT_DIRECTORY_FOLDERNAME, GPT_SEO_OUTPUT_DIRECTORY_FOLDERNAME
from backend.geneguide.pipeline.gpt_descriptions.constants import (
    GPT_GENEGUIDE_DESCRIPTION_SYSTEM_ROLE,
    GPT_GOTERM_DESCRIPTION_USER_ROLE,
    GPT_GOTERM_SEO_DESCRIPTION_USER_ROLE,
)
from backend.geneguide.pipeline.providers.openai_provider import OpenAIProvider
from backend.geneguide.pipeline.providers.s3_provider import S3Provider
from backend.geneguide.pipeline.utils import get_object_key

logger = logging.getLogger(__name__)


def generate_new_gpt_descriptions(all_go_term_ids: dict[str, str], generate_all: bool = False) -> dict[str, str]:
    """
    Generates GPT descriptions for all GO terms in the corpus if they haven't
    already been generated. GO terms should be stored in the cellguide-data-public-<env>
    bucket in the gpt_descriptions folder. For each cell type ID, we'll have a JSON file
    that contains the GPT description

    📁 gpt_descriptions
    ﹂GO_0000000.json
    ﹂GO_0000001.json
    ﹂...
    """
    bucket_name = GeneGuideConfig().bucket

    openai_provider = OpenAIProvider()

    new_gpt_descriptions: dict[str, str] = {}

    s3_provider = S3Provider() if not generate_all else None

    for go_term in all_go_term_ids:
        go_term_filename = go_term.replace(":", "_")
        object_key = get_object_key(object=f"{GPT_OUTPUT_DIRECTORY_FOLDERNAME}/{go_term_filename}.json")
        if not s3_provider or not s3_provider.does_object_exist(bucket_name=bucket_name, object_key=object_key):
            go_term_label = all_go_term_ids[go_term]
            logger.info(f"GPT description does not exist yet for {go_term} ({go_term_label}), generating...")
            description_str = openai_provider.generate_gpt_output(
                user_role=GPT_GOTERM_DESCRIPTION_USER_ROLE(go_term_label),
                system_role=GPT_GENEGUIDE_DESCRIPTION_SYSTEM_ROLE,
            )
            new_gpt_descriptions[go_term] = description_str
            logger.info(f"Generated GPT description for {go_term}: {description_str}")

            # extra safety to avoid rate limiting
            time.sleep(1)

    return new_gpt_descriptions


def generate_new_seo_gpt_descriptions(
    new_gpt_descriptions: dict[str, str], generate_all: bool = False
) -> dict[str, str]:
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
    bucket_name = GeneGuideConfig().bucket
    s3_provider = S3Provider() if not generate_all else None
    openai_provider = OpenAIProvider()

    new_gpt_seo_descriptions: dict[str, str] = {}
    for go_term in new_gpt_descriptions:
        go_term_filename = go_term.replace(":", "_")
        object_key = get_object_key(object=f"{GPT_SEO_OUTPUT_DIRECTORY_FOLDERNAME}/{go_term_filename}.json")
        if not s3_provider or not s3_provider.does_object_exist(bucket_name=bucket_name, object_key=object_key):
            logger.info(f"GPT SEO description does not exist yet for {go_term}, generating...")
            seo_description_str = openai_provider.generate_gpt_output(
                user_role=GPT_GOTERM_SEO_DESCRIPTION_USER_ROLE(new_gpt_descriptions[go_term]),
                system_role=GPT_GENEGUIDE_DESCRIPTION_SYSTEM_ROLE,
            )
            new_gpt_seo_descriptions[go_term] = seo_description_str
            logger.info(f"Generated GPT description for {go_term}: {seo_description_str}")

            # extra safety to avoid rate limiting
            time.sleep(1)

    return new_gpt_seo_descriptions
