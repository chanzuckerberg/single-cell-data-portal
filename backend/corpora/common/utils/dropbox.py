from urllib.parse import urlparse, parse_qs, urlencode, urlunparse

import regex
import requests


def get_download_url_from_shared_link(url):
    parsed_url = urlparse(url)
    query = parse_qs(parsed_url.query)
    params = {"dl": "1"}
    query.update(params)

    url_parts = list(parsed_url)
    url_parts[4] = urlencode(query)
    result = urlunparse(url_parts)
    return result


dropbox_link_rx = regex.compile(r"https://www.dropbox.com/s/[\w\d]+/.*\?dl=[01]")


def verify(url):
    return True if dropbox_link_rx.fullmatch(url) else False


def get_file_info(url):
    resp = requests.head(url, allow_redirects=True)
    return {
        "size": int(resp.headers["content-length"]),
        "name": resp.headers["content-disposition"].split(";")[1].split("=", 1)[1][1:-1],
    }
