from typing import Optional

from pydantic import BaseModel, HttpUrl


class IngestionManifest(BaseModel):
    """
    # Deserialize JSON to Pydantic model
    data_obj = IngestionManifest.model_validate_json(json_data)

    # Convert Pydantic object to JSON
    json_output = data_obj.model_dump_json(indent=2)
    """

    anndata: HttpUrl
    atac_seq_fragment: Optional[HttpUrl] = None  # Optional field
