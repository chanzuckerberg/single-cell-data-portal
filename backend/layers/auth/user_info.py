from typing import Optional

from backend.common.utils import authorization_checks as auth


class UserInfo:
    def __init__(self, token_info: Optional[dict]) -> None:
        self.token_info = token_info

    def is_none(self):
        return not self.token_info

    def is_super_curator(self):
        if self.token_info is None:
            return False
        else:
            return auth.is_super_curator(self.token_info.get("scope", ""))

    def is_cxg_admin(self):
        if self.token_info is None:
            return False
        else:
            return auth.is_cxg_admin(self.token_info.get("scope", ""))

    @property
    def user_id(self):
        if self.token_info is None:
            raise Exception("Bad")  # TODO: improve
        return self.token_info.get("sub")

    def is_user_owner_or_allowed(self, owner: str) -> bool:
        if self.is_none():
            return False
        return self.user_id == owner or self.is_super_curator()
