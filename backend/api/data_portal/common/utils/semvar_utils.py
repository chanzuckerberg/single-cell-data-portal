import re

# Official SemVer regex: https://semver.org/
SEMVER_FORMAT = re.compile(
    r"^(?P<major>0|[1-9]\d*)\.(?P<minor>0|[1-9]\d*)\.(?P<patch>0|[1-9]\d*)(?:-(?P<prerelease>(?:0|[1-9]\d*|\d*["
    r"a-zA-Z-][0-9a-zA-Z-]*)(?:\.(?:0|[1-9]\d*|\d*[a-zA-Z-][0-9a-zA-Z-]*))*))?(?:\+(?P<buildmetadata>[0-9a-zA-Z-]+("
    r"?:\.[0-9a-zA-Z-]+)*))?$"
)


def validate_version_str(version_str, release_only=True):
    """
    Test if a string conforms to SemVer format (https://semver.org/)
    :param version_str: a string to be validated
    :param release_only: only declare releases (not prereleases) valid
    :return: True if the version string is of a valid SemVer format else False
    """

    match = SEMVER_FORMAT.match(version_str)
    has_match = match is not None
    if has_match and release_only:
        return not match.group("prerelease")
    return has_match
