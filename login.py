import base64
import requests
import json
import re
from urllib import parse

requests.packages.urllib3.disable_warnings(requests.packages.urllib3.exceptions.InsecureRequestWarning)

OIDC_LOGIN_URL="http://localhost:5000/dp/v1/login"
USERS_CONFIGURATION_PATH="oauth/users.json"
CLIENTS_CONFIGURATION_PATH="oauth/clients-config.json"

f = open(USERS_CONFIGURATION_PATH)
data = json.load(f) 

user_data = data[0]
user_name = user_data['Username']
password = user_data['Password']

f = open(CLIENTS_CONFIGURATION_PATH)
data = json.load(f)
client_data = data[0]
client_id = client_data['ClientId']
client_secret = client_data['ClientSecrets'][0]


# Calling http://localhost:5000/dp/v1/login
session = requests.session()
response = session.get(OIDC_LOGIN_URL, verify=False, allow_redirects = False)
location = response.headers['Location']
login_cookies = response.cookies
params = dict(parse.parse_qsl(parse.urlsplit(location).query))


# Calling /connect/authorize/
session = requests.session()
location = location.replace("https://localhost:8443", "https://oidc")
response = session.get(location, verify=False, allow_redirects = False, cookies=login_cookies)
location = response.headers['Location']

# Calling Account/login
session = requests.session()
response = session.get(location, verify=False, allow_redirects = False, cookies = login_cookies)
authorize_cookies = response.cookies

# Parsining login page
page_content = response.content

reg = re.compile('input name="__RequestVerificationToken" type="hidden" value="([A-Za-z0-9_\-]*)" />')
verification_token = reg.search(page_content.decode('utf-8')).groups()[0]

# calling /Account/login with all the login details.
session = requests.session()
session.headers.update({'content-type': 'application/x-www-form-urlencoded'})

login_payload = {
  'Username': user_name,
  'Password': password,
  '__RequestVerificationToken': verification_token,
  'button': 'login',
  'RememberLogin': False,
  'client_id': client_id,
  'scope': 'openid profile email offline_access',
  'state': params['state'],
  'nonce': params['nonce'],
  'code_challenge': params['code_challenge'],
  'code_challenge_method': params['code_challenge_method'],
  'response_type': params['response_type'],
  'redirect_uri': params['redirect_uri']
}

combine_cookies = {}
for c in login_cookies:
  combine_cookies[c.name] = c.value
for c in authorize_cookies:
  combine_cookies[c.name] = c.value

response = session.post(location, data=login_payload, verify=False, allow_redirects = False, cookies = combine_cookies)

# Calling /connect/authorize/
location = response.headers['Location']
for c in response.cookies:
  combine_cookies[c.name] = c.value

session = requests.session()
response = session.get("https://oidc" + location, verify=False, allow_redirects = False, cookies=combine_cookies)
location = response.headers['Location']


location = response.headers['Location']
session = requests.session()
response = session.get(location, verify=False, allow_redirects = False, cookies=login_cookies)

# Get the final cookie with which we can make further API calls.
cxguser_cookie = response.cookies

print(f"cxguser_cookie:{cxguser_cookie}")

# Sample call to get userinfo
response = requests.get("http://localhost:5000/dp/v1/userinfo", cookies=cxguser_cookie)
print(response, response.text)

# Sample call to get collections
response = requests.get("http://localhost:5000/dp/v1/collections", cookies=cxguser_cookie)
print(response, response.text)
