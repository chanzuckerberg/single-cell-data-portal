import subprocess
import urllib


def fix_dropbox_url(url):
    """Fix a dropbox url so it's a direct download. If it's not a valid dropbox url, return None."""

    pr = urllib.parse.urlparse(url)

    if pr.scheme != "https":
        return None

    if pr.netloc != "www.dropbox.com":
        return None

    if "dl=0" in pr.query:
        new_query = pr.query.replace("dl=0", "dl=1")
    elif not pr.query:
        new_query = "dl=1"
    elif "dl=1" in pr.query:
        new_query = pr.query
    else:
        new_query = pr.query + "&dl=1"

    pr = pr._replace(query=new_query)

    return pr.geturl()


def fetch_dropbox_url(dropbox_url, local_path):
    """Given a dropbox url, download it to local_path.

    Handles fixing the url so it downloads directly.
    """

    fixed_dropbox_url = fix_dropbox_url(dropbox_url)

    if not fixed_dropbox_url:
        raise ValueError(f"Malformed Dropbox URL: {dropbox_url}")

    subprocess.run(["wget", fixed_dropbox_url, "-O", local_path], check=True)

    return local_path
