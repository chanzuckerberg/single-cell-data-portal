class PortalApiInterface:

    def get_collections_list(from_date: int = None, to_date: int = None, token_info: Optional[dict] = None):
        pass

    def get_collection_details(collection_id: str, token_info: dict):
        pass

    def get_collections_index():
        pass

    def post_collection_revision(collection_id: str, token_info: dict):
        pass

    def create_collection(body: dict, user: str):
        pass

    def delete_collection(collection_id: str, token_info: dict):
        pass

    def update_collection(collection_id: str, body: dict, token_info: dict):
        pass


    def post(collection_id: str, body: object, token_info: dict): # publish
        pass

    def link(collection_id: str, body: dict, token_info: dict):
        pass

    def relink(collection_id: str, body: dict, token_info: dict):
        pass

    def upload_from_link(collection_id: str, token_info: dict, url: str, dataset_id: str = None):
        pass



    def post_dataset_asset(dataset_id: str, asset_id: str):
        pass

    def get_dataset_assets(dataset_id: str):
        pass

    def get_status(dataset_id: str, token_info: dict):
        pass

    def get_datasets_index():
        pass

    def delete_dataset(dataset_id: str, token_info: dict):
        pass

    def get_dataset_identifiers(url: str):
        pass

    def post_dataset_gene_sets(dataset_id: str, body: object, token_info: dict):
        pass
