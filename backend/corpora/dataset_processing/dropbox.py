import os
import subprocess
import urllib

import requests


def fix_dropbox_url(url):
    pr = urllib.parse.urlparse(url)

    if pr.scheme != "https":
        return None

    if pr.netloc != "www.dropbox.com":
        return None

    if "dl=0" in pr.query:
        new_query = pr.query.replace("dl=0", "dl=1")
    elif not pr.query:
        new_query = "dl=1"
    else:
        new_query = pr.query + "?dl=1"

    pr = pr._replace(query=new_query)

    return pr.geturl()


def fetch_dropbox_url(dropbox_url):

    fixed_dropbox_url = fix_dropbox_url(dropbox_url)

    if not fixed_dropbox_url:
        raise ValueError(f"Malformed Dropbox URL: {dropbox_url}")

    subprocess.run(["wget", fixed_dropbox_url, "-O", "local.h5ad"], check=True)

    resp = requests.head(fixed_dropbox_url, allow_redirects=True)
    content_length = resp.headers["Content-Length"]

    if os.path.getsize("local.h5ad") != content_length:
        raise RuntimeError("Downloaded file isn't the correct size!")

    return "local.h5ad"
