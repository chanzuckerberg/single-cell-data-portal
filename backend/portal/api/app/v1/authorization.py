from typing import Optional

from backend.common.utils.authorization_checks import (
    is_user_owner_or_allowed as is_user_owner_or_allowed_common,
    owner_or_allowed as owner_or_allowed_common,
    is_super_curator as is_super_curator_common,
)


def is_user_owner_or_allowed(token_info: dict, owner: str) -> bool:
    return is_user_owner_or_allowed_common(token_info.get("sub"), token_info.get("scope", ""), owner)


def owner_or_allowed(token_info: dict) -> Optional[str]:
    return owner_or_allowed_common(token_info.get("sub"), token_info.get("scope", ""))


def is_super_curator(token_info: dict) -> bool:
    return is_super_curator_common(token_info.get("scope", ""))
