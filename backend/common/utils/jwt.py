from jose import ExpiredSignatureError, jwt
from jose.exceptions import JWTClaimsError, JWTError

from backend.common.utils.http_exceptions import ExpiredCredentialsError, UnauthorizedError


def jwt_decode(*args, **kwargs) -> dict:
    try:
        return jwt.decode(*args, **kwargs)
    except ExpiredSignatureError:
        raise ExpiredCredentialsError(
            detail="Your authentication credentials have expired. Please provide updated "
            "credentials before trying again."
        )
    except JWTClaimsError:
        raise UnauthorizedError(detail="Incorrect claims, please check the audience and issuer.")
    except Exception:
        raise UnauthorizedError(detail="Unable to parse authentication token.")


def get_unverified_header(token: str) -> dict:
    try:
        return jwt.get_unverified_header(token)
    except JWTError:
        raise UnauthorizedError(detail="Unable to parse authentication token.")
