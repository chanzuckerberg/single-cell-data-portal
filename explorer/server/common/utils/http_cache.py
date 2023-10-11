from functools import wraps

from flask import Blueprint, current_app, make_response

webbp = Blueprint("webapp", "server.common.web", template_folder="templates")
ONE_WEEK = 7 * 24 * 60 * 60
ONE_YEAR = 365 * 24 * 60 * 60


def _cache_control(always, **cache_kwargs):
    """
    Used to easily manage cache control headers on responses.
    See Werkzeug for attributes that can be set, eg, no_cache, private, max_age, etc.
    https://werkzeug.palletsprojects.com/en/1.0.x/datastructures/#werkzeug.datastructures.ResponseCacheControl
    """

    def inner_cache_control(f):
        @wraps(f)
        def wrapper(*args, **kwargs):
            response = make_response(f(*args, **kwargs))
            if not always and not current_app.app_config.server__app__generate_cache_control_headers:
                return response
            if response.status_code >= 400:
                return response
            for k, v in cache_kwargs.items():
                setattr(response.cache_control, k, v)
            return response

        return wrapper

    return inner_cache_control


def cache_control(**cache_kwargs):
    """config driven"""
    return _cache_control(False, **cache_kwargs)


def cache_control_always(**cache_kwargs):
    """always generate headers, regardless of the config"""
    return _cache_control(True, **cache_kwargs)
