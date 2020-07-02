def get_submissions_list(query_user_uuid: str):
    raise NotImplementedError


def create_new_submission(request_body: dict):
    raise NotImplementedError


def get_submission_details(path_project_uuid: str):
    raise NotImplementedError


def delete_submission(path_project_uuid: str):
    raise NotImplementedError


def add_file_to_submission(path_project_uuid: str, request_body: dict):
    raise NotImplementedError


def delete_dataset_from_submission(path_project_uuid: str, path_dataset_uuid: str):
    raise NotImplementedError


def validate_submission(path_project_uuid: str):
    raise NotImplementedError


def save_submission(path_project_uuid: str, request_body: dict):
    raise NotImplementedError


def publish_submission(path_project_uuid: str):
    raise NotImplementedError
