import time
from jose import jws, jwt
from jose.constants import ALGORITHMS


seconds_in_day = 86400
days = 60


def generate(user_id, secret):
    iat = time.time()
    exp = iat + (seconds_in_day * days)
    sub = user_id
    return jws.sign({"exp": exp, "ait": iat, "sub": sub}, secret, algorithm=ALGORITHMS.HS256)


def verify(token, secret):
    jwt.decode(token, secret, algorithms=["HS256"])
