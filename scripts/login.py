import json
import re
from urllib import parse

import requests  # type: ignore


def get_cxguser_cookie():
    # suppress https warnings.
    requests.packages.urllib3.disable_warnings(requests.packages.urllib3.exceptions.InsecureRequestWarning)

    OIDC_LOGIN_URL = "https://backend.corporanet.local:5000/dp/v1/login"
    USERS_CONFIGURATION_PATH = "oauth/users.json"
    CLIENTS_CONFIGURATION_PATH = "oauth/clients-config.json"

    with open(USERS_CONFIGURATION_PATH) as f:
        data = json.load(f)

    user_data = data[0]
    user_name = user_data["Username"]
    password = user_data["Password"]

    with open(CLIENTS_CONFIGURATION_PATH) as f:
        data = json.load(f)
    client_data = data[0]
    client_id = client_data["ClientId"]

    # Calling https://backend.corporanet.local:5000/dp/v1/login
    session = requests.session()
    response = session.get(OIDC_LOGIN_URL, verify=False, allow_redirects=False)
    location = response.headers["Location"]
    login_cookies = response.cookies
    params = dict(parse.parse_qsl(parse.urlsplit(location).query))

    # Calling /connect/authorize/
    location = location.replace("https://localhost:8443", "https://oidc")
    response = session.get(location, verify=False, allow_redirects=False)
    location = response.headers["Location"]

    # Calling Account/login
    response = session.get(location, verify=False, allow_redirects=False)

    # Parsing login page
    page_content = response.content
    reg = re.compile(r'input name="__RequestVerificationToken" type="hidden" value="([A-Za-z0-9_\-]*)" />')
    verification_token = reg.search(page_content.decode("utf-8")).groups()[0]

    login_payload = {
        "Username": user_name,
        "Password": password,
        "__RequestVerificationToken": verification_token,
        "button": "login",
        "RememberLogin": False,
        "client_id": client_id,
        "scope": "openid profile email offline_access",
        "state": params["state"],
        "nonce": params["nonce"],
        "code_challenge": params["code_challenge"],
        "code_challenge_method": params["code_challenge_method"],
        "response_type": params["response_type"],
        "redirect_uri": params["redirect_uri"],
    }

    # calling /Account/login with all the login details.
    response = session.post(location, data=login_payload, verify=False, allow_redirects=False)

    # Calling /connect/authorize/
    location = response.headers["Location"]
    response = session.get("https://oidc" + location, verify=False, allow_redirects=False)
    location = response.headers["Location"]

    session.cookies = login_cookies
    response = session.get(location, verify=False, allow_redirects=False)

    # final cookie with which we can make further API calls.
    return response.cookies


if __name__ == "__main__":
    cxguser_cookie = get_cxguser_cookie()
    print(f"cxguser_cookie:{cxguser_cookie}")
