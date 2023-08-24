import re

USERNAME_REGEX = r"(?P<username>[\w\-\|]+)"
ID_REGEX = r"[0-9a-fA-F]{8}-[0-9a-fA-F]{4}-[0-9a-fA-F]{4}-[0-9a-fA-F]{4}-[0-9a-fA-F]{12}"
DATASET_ID_REGEX = f"(?P<dataset_id>{ID_REGEX})"
COLLECTION_ID_REGEX = f"(?P<collection_id>{ID_REGEX})"
CONTROL_CHARS = r"[\x00-\x1f\x7f\x81\x8d\x8f\x90\x9d\xa0]"

DOI_REGEX_COMPILED = re.compile(r"^10.\d{4,9}/[-._;()/:A-Z0-9]+$", flags=re.I)
CURIE_REFERENCE_REGEX = r"^\d{2}\.\d{4}.*$"

EMAIL_REGEX = re.compile(r"(.+)@(.+)\.(.+)")
S3_URI_REGEX = re.compile(r"^s3://(?P<bucket>[^/]+)/((?P<prefix>.*/$)|(?P<key>.*))")
