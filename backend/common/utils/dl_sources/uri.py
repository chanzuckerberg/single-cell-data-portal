import shutil
import typing
from abc import ABC, abstractmethod
from urllib.parse import ParseResult, urlparse

import boto3
import requests

from backend.common.utils import downloader
from backend.layers.thirdparty.s3_provider import S3Provider


class DownloadFailed(Exception):
    pass


class MissingHeaderException(Exception):
    def __init__(self, detail: str = "", *args, **kwargs) -> None:
        self.detail = "Missing header from response. " + detail


class URI(ABC):
    """Define the abstract base class to support different download sources."""

    def __init__(self, uri, parsed_uri: ParseResult):
        self.uri: str = uri
        self.parsed_uri: ParseResult = parsed_uri

    @classmethod
    @abstractmethod
    def validate(cls, uri: str) -> typing.Optional["URI"]:
        """Validates the URI matches the expected format, and returns a new class object if valid.."""

    @abstractmethod
    def file_info(self) -> dict:
        """
        Extract information about a file from a URI.
        """

    @property
    def scheme(self):
        return self.parsed_uri.scheme

    @property
    def netloc(self):
        return self.parsed_uri.netloc

    @property
    def path(self):
        return self.parsed_uri.path

    def _get_key(self, headers: dict, key: str) -> str:
        try:
            return headers[key]
        except KeyError:
            raise MissingHeaderException(
                f"{self.__class__.__name__}:URI({self.uri}) failed request. '{key}' not present in the header."
            ) from None

    def _get_key_with_fallback(self, headers: dict, key: str, fallback_key: str) -> str:
        try:
            return headers.get(key) or headers[fallback_key]
        except KeyError:
            raise MissingHeaderException(
                f"""{self.__class__.__name__}:URI({self.uri}) failed request.
                Neither '{key}' nor '{fallback_key}' are present in the header.
                """
            ) from None

    def _disk_space_check(self):
        file_size = self.file_info()["size"]
        if file_size and file_size >= shutil.disk_usage("/")[2]:
            raise DownloadFailed("Insufficient disk space.")

    @abstractmethod
    def download(self, local_file_name: str) -> None:
        pass


class DropBoxURL(URI):
    """Supports download URIs from a DropBox share link."""

    @classmethod
    def validate(cls, uri: str) -> typing.Optional["URI"]:
        """Converts a valid DropBox URI into a direct download link. If the uri is not a valid DropBox URI, none is
        returned. Otherwise, the converted URI is returned.
        """

        parsed_uri = urlparse(uri)
        if parsed_uri.scheme != "https" or parsed_uri.netloc != "www.dropbox.com":
            return None
        # dl=0 will show the file in the preview page. A link with ? dl=1 will force the file to download.
        if "dl=0" in parsed_uri.query:
            new_query = parsed_uri.query.replace("dl=0", "dl=1")
        elif not parsed_uri.query:
            new_query = "dl=1"
        elif "dl=1" in parsed_uri.query:
            new_query = parsed_uri.query
        else:
            new_query = parsed_uri.query + "&dl=1"

        parsed_uri = parsed_uri._replace(query=new_query)
        return cls(parsed_uri.geturl(), parsed_uri)

    def file_info(self) -> dict:
        """
        Extract information about a file from a DropBox URI.
        :param uri: a DropBox URI leading to a file.
        :return: The file name and size of the file.
        """
        try:
            resp = requests.head(self.uri, allow_redirects=True)
            resp.raise_for_status()
        except Exception:
            return {"size": None, "name": None}

        try:
            size = int(self._get_key_with_fallback(resp.headers, "content-length", "x-dropbox-content-length"))
        except Exception:
            size = None

        return {
            "size": size,
            "name": self._get_key(resp.headers, "content-disposition").split(";")[1].split("=", 1)[1][1:-1],
        }

    def download(self, local_file_name: str):
        self._disk_space_check()
        downloader.download(self.uri, local_file_name)


class S3URL(URI):
    """Supports presigned URLs from an AWS S3 bucket."""

    _netloc = "amazonaws.com"
    _scheme = "https"

    @classmethod
    def validate(cls, uri: str):
        parsed_uri = urlparse(uri)
        return (
            cls(uri, parsed_uri)
            if parsed_uri.scheme == cls._scheme and parsed_uri.netloc.endswith(cls._netloc)
            else None
        )

    def file_info(self) -> dict:
        resp = requests.get(self.uri, headers={"Range": "bytes=0-0"})
        resp.raise_for_status()

        return {
            "size": int(self._get_key(resp.headers, "content-range").split("/")[1]),
            "name": self.parsed_uri.path,
        }

    def download(self, local_file_name: str):
        self._disk_space_check()
        downloader.download(self.uri, local_file_name)


class S3URI(URI):
    """
    Handles S3 URIs: s3://<bucket>/<key>
    """

    def __init__(self, uri, parsed_uri: ParseResult):
        super().__init__(uri, parsed_uri)
        self.s3_provider = S3Provider()

    @classmethod
    def validate(cls, uri: str) -> typing.Optional["URI"]:
        parsed = urlparse(uri)
        bucket_name = parsed.netloc
        key = parsed.path
        if parsed.scheme == "s3" and bucket_name and key:
            return cls(parsed.geturl(), parsed)
        else:
            return None

    def file_info(self) -> dict:
        s3 = boto3.resource("s3")
        s3_object = s3.Object(self.bucket_name, self.key)
        return {"name": self.key, "size": s3_object.content_length}

    @property
    def bucket_name(self):
        return self.parsed_uri.netloc

    @property
    def key(self):
        return self.parsed_uri.path

    def download(self, local_file_name: str):
        self._disk_space_check()
        self.s3_provider.download_file(self.bucket_name, self.key, local_file_name)


class RegisteredSources:
    """Manages all of the download sources."""

    _registered: typing.Set[typing.Type[URI]] = set()

    @classmethod
    def add(cls, parser: typing.Type[URI]):
        if issubclass(parser, URI):
            cls._registered.add(parser)
        else:
            raise TypeError(f"subclass type {URI.__name__} expected")

    @classmethod
    def remove(cls, parser: typing.Type[URI]):
        cls._registered.remove(parser)

    @classmethod
    def get(cls) -> typing.Iterable:
        return cls._registered


def from_uri(uri: str) -> typing.Optional[URI]:
    """Given a URI return a object that can be used by the processing container to download data."""
    for source in RegisteredSources.get():
        uri_obj = source.validate(uri)
        if uri_obj:
            return uri_obj
    return None


# RegisteredSources are processed in the order registered and returns the first match.
RegisteredSources.add(DropBoxURL)
RegisteredSources.add(S3URL)
RegisteredSources.add(S3URI)
