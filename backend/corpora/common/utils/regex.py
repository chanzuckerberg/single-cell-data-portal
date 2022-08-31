import re

USERNAME_REGEX = r"(?P<username>[\w\-\|]+)"
ID_REGEX = r"[0-9a-fA-F]{8}-[0-9a-fA-F]{4}-[0-9a-fA-F]{4}-[0-9a-fA-F]{4}-[0-9a-fA-F]{12}"
DATASET_ID_REGEX = f"(?P<dataset_id>{ID_REGEX})"
COLLECTION_ID_REGEX = f"(?P<collection_id>{ID_REGEX})"
CONTROL_CHARS = r"[\x00-\x1f\x7f-\xa0]"
CURATOR_TAG_REGEX = r"(?P<tag>.+)"


def validate_curator_tag(curator_tag: str) -> bool:
    """
    Verify the correct curator tag format is obeyed (i.e., it is not a UUID)

    :param curator_tag: the tag name to validate.
    :return: True if CURATOR_TAG_REGEX is matched.
    """
    regex = f"^({DATASET_ID_REGEX}|{CURATOR_TAG_REGEX})$"
    matched = re.match(regex, curator_tag)
    if matched:
        matches = matched.groupdict()
        if matches.get("tag"):
            return True
        elif matches.get("dataset_id"):
            raise ValueError("Curator tag cannot assume UUID format.")
    raise ValueError("Curator tag cannot be empty.")
