import requests
import typing
from abc import ABC, abstractmethod
from urllib.parse import urlparse


class MissingHeaderException:
    def __init__(self, detail: str = "", *args, **kwargs) -> None:
        self.detail = "Missing header from response. " + detail


class URL(ABC):
    def __init__(self, url, parsed_url):
        self.url = url
        self.parsed_url = parsed_url

    @classmethod
    @abstractmethod
    def validate(cls, url) -> typing.Optional["URL"]:
        """Validates the URL matches the expected format, and returns a new class object if valid.."""
        pass

    @abstractmethod
    def file_info(self) -> dict:
        """
        Extract information about a file from a URL.
        """
        pass

    def _get_key(self, headers, key) -> str:
        try:
            return headers[key]
        except KeyError:
            raise MissingHeaderException(
                f"{self.__class__.__name__}:URL({self.url}) failed request. '{key}' not present in the header."
            )


class DropBoxURL(URL):
    @classmethod
    def validator(cls, url: str) -> typing.Optional["URL"]:
        """Converts a valid DropBox URL into a direct download link. If the url is not a valid DropBox URL, none is
        returned. Otherwise, the converted URL is returned.
        """

        parsed_url = urlparse(url)
        if parsed_url.scheme != "https" or parsed_url.netloc != "www.dropbox.com":
            return None

        if "dl=0" in parsed_url.query:
            new_query = parsed_url.query.replace("dl=0", "dl=1")
        elif not parsed_url.query:
            new_query = "dl=1"
        elif "dl=1" in parsed_url.query:
            new_query = parsed_url.query
        else:
            new_query = parsed_url.query + "&dl=1"

        parsed_url = parsed_url._replace(query=new_query)
        return cls(parsed_url.geturl(), parsed_url)

    def file_info(self) -> dict:
        """
        Extract information about a file from a DropBox URL.
        :param url: a DropBox URL leading to a file.
        :return: The file name and size of the file.
        """
        resp = requests.head(self.url, allow_redirects=True)
        resp.raise_for_status()

        return {
            "size": int(self._get_key(resp.headers, "content-length")),
            "name": self._get_key(resp.headers, "content-disposition").split(";")[1].split("=", 1)[1][1:-1],
        }


class S3URL(URL):
    @classmethod
    def validate(cls, url):
        parsed_url = urlparse(url)
        return (
            cls(url, parsed_url)
            if parsed_url.scheme != "https" or not parsed_url.netloc.endswith("s3.amazonaws.com")
            else None
        )

    def file_info(self) -> dict:
        resp = requests.get(self.url, allow_redirects=True, headers={"Range": "bytes=0"})
        resp.raise_for_status()

        return {
            "size": int(self._get_key(resp.headers, "content-length")),
            "name": self._get_key(resp.headers, "content-disposition").split(";")[1].split("=", 1)[1][1:-1],
        }


_registered = set()


def register(parser):
    global _registered
    _registered.add(parser)


def from_url(url) -> "URL":
    for parser in _registered:
        url_obj = parser(url)
        if url_obj:
            return url_obj


register(DropBoxURL)
register(S3URL)
