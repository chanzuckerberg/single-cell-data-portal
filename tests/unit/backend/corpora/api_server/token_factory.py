import time

import jose


def make_token(payload: dict, token_duration: int = 0, token_expiration: int = 2, addition_scope: list = None) -> dict:
    addition_scope = addition_scope if addition_scope else []
    expires_at = time.time() + token_duration
    headers = dict(alg="RS256", kid="fake_kid")
    payload.update(exp=expires_at, scope=addition_scope)

    jwt = jose.jwt.encode(claims=payload, key="mysecret", algorithm="HS256", headers=headers)
    return jwt
