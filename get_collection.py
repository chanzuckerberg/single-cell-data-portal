import requests

token = "eyJhbGciOiJSUzI1NiIsInR5cCI6IkpXVCIsImtpZCI6Ik5iRUNKeDlTVEhzNnEyX1JoNEVFLSJ9.eyJpc3MiOiJodHRwczovL2N6aS1jZWxseGdlbmUtZGV2LnVzLmF1dGgwLmNvbS8iLCJzdWIiOiJnb29nbGUtb2F1dGgyfDExNTg5ODU5MDIyODY2Mjg3ODYzMCIsImF1ZCI6Imh0dHBzOi8vYXBpLmNlbGx4Z2VuZS5kZXYuc2luZ2xlLWNlbGwuY3ppLnRlY2hub2xvZ3kvZHAvdjEvY3VyYXRvciIsImlhdCI6MTY1MTEzMTA2MSwiZXhwIjoxNjUxMjE3NDYxLCJhenAiOiJHZUNTaVdyOFNBbnQwbFlWbVF1UktjaXJrcVMyRFRjeSIsInNjb3BlIjoicHJvZmlsZSBlbWFpbCIsImd0eSI6InBhc3N3b3JkIn0.I5z0ME6ZutPULu9RFTlQiKXblfZen9DKslL1O9pApLvdt3fzXIhDMyvaA4wC9KWvjl2hb4swoo4oQNRT_iI7pVN9IEH_fK-v7P7SNUAalptiWjeEfRT54uBgktO_hHDxpuVa2SVEAUPKEbIRXlye7U_qhsmlAAllA56kIxIAfeQ73ug2r1X1l_HB76O8iCdBb5X9_f7ZEWpcshim_Mc21WO_Ad70ohtjpHu_9J-zVMJ4v72Mlw0qAxxvFiIIUxocprVdyZeyOXZq2Q6xqrOxCz4OWnv5QuSRqKq0R0w56pWzN9AvumFpU-VBwhBxjdMyKoAsLbA_PWBHYA88uxXk7g"

# url = "https://api.cellxgene.dev.single-cell.czi.technology/curation/v1/collections/77f91594-9d5e-48cd-b3f8-48c2fbbd36ff/datasets/s3-upload-credentials"
# url = "https://dan-rdev-s3-backend.rdev.single-cell.czi.technology/curation/v1/collections/77f91594-9d5e-48cd-b3f8-48c2fbbd36ff/datasets/s3-upload-credentials"
url = (
    "https://api.cellxgene.dev.single-cell.czi.technology/curation/v1/collections/0076fb88-073b-4436-ae1e-abae6868f885"
)

# headers = {"Authorization": f"{token}"}
headers = {}
resp = requests.get(url, headers=headers)
print(resp.json())
