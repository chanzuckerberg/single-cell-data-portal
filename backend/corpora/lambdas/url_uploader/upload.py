import boto3
import requests
from urllib import parse

s3 = boto3.client("s3")


def upload_from_url(url: str, file_name: str):
    with requests.get(url, stream=True, allow_redirects=True) as r:
        r.raise_for_status()
        file_content = r.raw
        s3.upload_fileobj(file_content, bucket_name, file_name)


def download_file(url: str, loca_file_name: str):
    # NOTE the stream=True parameter below
    with requests.get(url, stream=True, allow_redirects=True) as r:
        r.raise_for_status()
        with open(loca_file_name, "wb") as f:
            for chunk in r.iter_content(chunk_size=8192):
                # If you have chunk encoded response uncomment if
                # and set chunk_size parameter to None.
                # if chunk:
                f.write(chunk)
    return loca_file_name


url_readme = ("https://www.dropbox.com/s/3rer8zevxprds8o/README.md?dl=1", "readme.md")
url_18MB = (
    "https://www.dropbox.com/s/wd8ydc6ffiikqgv/The%20Complete%20Visual%20Guide%20to%20Building%20a%20House%20by%20Chuck%20Lockhart%20and%20John%20Carroll.epub?dl=1",
    "house.epub",
)
url_450MB = ("https://www.dropbox.com/s/h0wz5mfvee2vtmv/IMG_0295.MOV?dl=1", "video.mov")
url_gdrive = (
    "https://drive.google.com/file/d/0B-I5F-7UYWNndVdUdGpud1h3V1lxRWNoQXd1eDFWOXpSOWF3/view?usp=sharing",
    "gdrive.mp4",
)
"https://drive.google.com/file/d/0B-I5F-7UYWNnaG9tc2p6OW83SzlVTzVEZGYyX3J4Qy00RWtN/view?usp=sharing"
"https://drive.google.com/file/d/0B-I5F-7UYWNnaG9tc2p6OW83SzlVTzVEZGYyX3J4Qy00RWtN/download?usp=sharing"
file_id = "0B-I5F-7UYWNnaG9tc2p6OW83SzlVTzVEZGYyX3J4Qy00RWtN"
f"https://drive.google.com/uc?export=download&id=0B-I5F-7UYWNnaG9tc2p6OW83SzlVTzVEZGYyX3J4Qy00RWtN"
f"https://drive.google.com/uc?export=download&id={file_id}"
bucket_name = "tsmith-url-upload"
# file_name = 'delete_me.mov'


# download_file(url_readme, "test_readme.md")
# upload_from_url(*url_readme)
# upload_from_url(*url_18MB)
upload_from_url(*url_gdrive)


def google_drive_get_download_url_from_shared_link(url):
    parsed_url = parse.urlparse(url)
    file_id = parsed_url.path.rsplit("/", maxsplit=1)[-1][-1]
    return f"https://drive.google.com/uc?export=download&id={file_id}"


def dropbox_get_download_url_from_shared_link(url):
    parsed_url = parse.urlparse(url)
    query = parse.parse_qs(parsed_url)
    params = {"dl": "1"}
    query.update(params)

    url_parts = list(parsed_url)
    url_parts[4] = parse.urlencode(query)
    result = parse.urlunparse(url_parts)
    return result
