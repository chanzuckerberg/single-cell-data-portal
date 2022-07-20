import time
from jose import jws, jwt
from jose.constants import ALGORITHMS


seconds_in_day = 86400


def generate(user_id: str, secret: str, days_to_live: int = 1) -> str:
    iat = time.time()
    exp = iat + (seconds_in_day * days_to_live)
    sub = user_id
    return jws.sign({"exp": exp, "iat": iat, "sub": sub}, secret, algorithm=ALGORITHMS.HS256)


def verify(token: str, secret: str):
    return jwt.decode(token, secret, algorithms=["HS256"])
