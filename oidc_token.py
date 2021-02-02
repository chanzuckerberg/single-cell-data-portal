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

print(client_id, client_secret, user_name, password)

token_payload = {
    'client_id': client_id, 
    'client_secret': client_secret, 
    'grant_type': 'password', 
    'Username': user_name, 
    'password': password
}

response = requests.post(OIDC_TOKEN_URL, data=token_payload)
json_data = json.loads(response.text)
print("access_token:", json_data["access_token"])
print("refresh_token:", json_data["refresh_token"])