from typing import Iterable, List


class S3ProviderInterface:
    def get_file_size(self, path: str) -> int:
        pass

    def generate_presigned_url(self, path: str) -> str:
        pass

    def upload_file(self, src_file: str, bucket_name: str, dst_file: str, extra_args: dict):
        pass

    def delete_files(self, bucket_name: str, object_keys: List[str]):
        pass

    def delete_prefix(self, bucket_name: str, prefix: str) -> None:
        pass

    def download_file(self, bucket_name: str, object_key: str, local_filename: str):
        pass

    def upload_directory(self, src_dir: str, s3_uri: str):
        pass

    def list_directory(self, bucket_name: str, src_dir: str) -> Iterable[str]:
        pass
