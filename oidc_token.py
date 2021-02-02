import base64
import requests
import json

OIDC_TOKEN_URL="http://oidc/connect/token"
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

token_payload = {
    'client_id': client_id, 
    'client_secret': client_secret, 
    'grant_type': 'password', 
    'Username': user_name, 
    'password': password,
}

response = requests.post(OIDC_TOKEN_URL, data=token_payload)
json_data = json.loads(response.text)

json_data["id_token"] = json_data["access_token"]
access_token = json_data["access_token"]
id_token = json_data["id_token"]
refresh_token = json_data["refresh_token"]


print(f"\naccess_token: {access_token}\n")
print(f"id_token: {id_token}\n")
print(f"refresh_token: {refresh_token}\n")

cookies = {
    'cxguser' : base64.b64encode(json.dumps(json_data).encode("utf-8")).decode()
}

response = requests.get("http://localhost:5000/dp/v1/collections", cookies=cookies)
json_data = json.loads(response.text)
print("Collections:", json_data['collections'])
