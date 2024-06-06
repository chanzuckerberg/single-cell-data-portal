def format_citation_dp(message: dict) -> str:
    """
    Formats the citation message.

    Parameters
    ----------
    message : dict
        The message containing publisher_metadata from the /collections API.

    Returns
    -------
    str
        The formatted citation string.
    """

    first_author = message["authors"][0]
    author_str = f"{first_author['family']} et al." if "family" in first_author else f"{first_author['name']} et al."

    journal = " " + message["journal"] if message["journal"] else ""
    year = f"{message['published_year']}"

    return f"{author_str} ({year}) {journal}"


def format_citation_crossref(message: dict) -> str:
    """
    Formats the citation message.

    Parameters
    ----------
    message : dict
        The message containing citation details output from CrossRef.

    Returns
    -------
    str
        The formatted citation string.
    """

    first_author = message["author"][0]
    author_str = f"{first_author['family']} et al." if "family" in first_author else f"{first_author['name']} et al."

    journal = " " + message["container-title"][0] if len(message["container-title"]) else ""
    year = message["created"]["date-parts"][0][0]

    return f"{author_str} ({year}) {journal}"
