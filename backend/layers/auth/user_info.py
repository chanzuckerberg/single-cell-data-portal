from typing import Optional

from backend.common.utils import authorization_checks as auth


class UserInfo:
    def __init__(self, token_info: Optional[dict]) -> None:
        self.token_info = token_info

    # def is_user_owner_or_allowed(token_info: dict, owner: str) -> bool:
    # return is_user_owner_or_allowed_common(token_info.get("sub"), token_info.get("scope", ""), owner)

    # def owner_or_allowed(token_info: dict) -> Optional[str]:
    #     return owner_or_allowed_common(token_info.get("sub"), token_info.get("scope", ""))

    # def is_super_curator(token_info: dict) -> bool:
    #     return is_super_curator_common(token_info.get("scope", ""))

    def is_none(self):
        return not self.token_info

    def is_super_curator(self):
        if self.token_info is None:
            return False
        else:
            return auth.is_super_curator(self.token_info.get("scope", ""))

    def user_id(self):
        if self.token_info is None:
            raise Exception("Bad")  # TODO: improve
        return self.token_info.get("sub")

    def is_user_owner_or_allowed(self, owner: str) -> bool:
        if self.is_none():
            return False
        return self.user_id() == owner or self.is_super_curator()

    # def get_collection_query_filter_for_owner(self, is_published: Optional[bool] = None):
    #     if self.is_super_curator():
    #         return CollectionQueryFilter(is_published, None)
    #     elif :
    #         return CollectionQueryFilter()
