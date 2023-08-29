import logging
import time

import openai

from backend.cellguide.pipeline.config import CellGuideConfig

logger: logging.Logger = logging.getLogger(__name__)


class OpenAIProvider:
    def __init__(self, model_name="gpt-4", fallback_model_name="gpt-3.5-turbo") -> None:
        self.organization = CellGuideConfig().openai_organization_id
        self.api_key = CellGuideConfig().openai_api_key

        valid_gpt_models = [
            i.id
            for i in openai.Model.list(api_key=self.api_key, organization=self.organization)["data"]
            if i.id.startswith("gpt")
        ]
        if model_name not in valid_gpt_models and fallback_model_name not in valid_gpt_models:
            raise ValueError(
                f"Model {model_name} and {fallback_model_name} not available. Valid models are {valid_gpt_models}"
            )

        elif model_name not in valid_gpt_models:
            logger.error(f"Model {model_name} not available, falling back to {fallback_model_name}")
            self.model_name = fallback_model_name
        else:
            self.model_name = model_name

    def generate_gpt_output(self, *, system_role: str, user_role: str, max_retries: int = 10) -> str:
        for i in range(max_retries):
            try:
                result = openai.ChatCompletion.create(
                    model=self.model_name,
                    messages=[
                        {"role": "system", "content": system_role},
                        {"role": "user", "content": user_role},
                    ],
                    api_key=self.api_key,
                    organization=self.organization,
                )
                return result["choices"][0]["message"]["content"]
            except Exception as e:
                logger.error(f"Error on attempt {i+1}: {e}")
                if i < max_retries - 1:  # Don't sleep after the last attempt
                    time.sleep(10)

        raise Exception("Failed to generate GPT output after 5 attempts")
