from urllib.parse import urlparse

import re
import requests

dropbox_link_rx = re.compile(r"/s/[\w\d]+/.*")


def get_download_url_from_shared_link(url):
    """Fix a dropbox url so it's a direct download. If it's not a valid dropbox url, return None."""

    parsed_url = urlparse(url)
    if (
        parsed_url.scheme != "https"
        or parsed_url.netloc != "www.dropbox.com"
        or not dropbox_link_rx.fullmatch(parsed_url.path)
    ):
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
    return parsed_url.geturl()


def get_file_info(url):
    resp = requests.head(url, allow_redirects=True)
    resp.raise_for_status()
    return {
        "size": int(resp.headers["content-length"]),
        "name": resp.headers["content-disposition"].split(";")[1].split("=", 1)[1][1:-1],
    }
