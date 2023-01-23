from dataclasses import dataclass
from backend.common.utils.dl_sources.url import URL, MissingHeaderException, from_url
import requests


@dataclass
class FileInfo:
    size: int
    name: str


class FileInfoException(Exception):
    pass


class UriProviderInterface:
    def validate(self, uri: str) -> bool:
        pass

    def parse(self, uri: str) -> URL:
        pass

    def get_file_info(self, uri: str) -> FileInfo:
        pass


class UriProvider(UriProviderInterface):
    def validate(self, uri: str) -> bool:
        """
        Returns true if and only if `uri` is a valid URI with respect to the registered sources.
        Currently we support:
        - Dropbox URLs
        - S3 URLs and URIs
        """
        link = from_url(uri)
        return link is not None

    def parse(self, uri: str) -> URL:
        """
        Returns a parsed URL object
        """
        return from_url(uri)

    def get_file_info(self, uri: str) -> FileInfo:
        """
        Returns the size and the name of the file specified in `uri`, as indicated by its provider
        """
        try:
            link = from_url(uri)
            file_info = link.file_info()

            return FileInfo(
                file_info["size"],
                file_info["name"],
            )

        except requests.HTTPError:
            raise FileInfoException("The URL provided causes an error with Dropbox.")
        except MissingHeaderException as ex:
            raise FileInfoException from ex
