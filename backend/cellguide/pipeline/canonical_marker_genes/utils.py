from backend.wmg.data.utils import setup_retry_session


def clean_doi(doi: str) -> str:
    """
    Cleans the DOI string.

    Parameters
    ----------
    doi : str
        The DOI string to be cleaned.

    Returns
    -------
    str
        The cleaned DOI string.
    """
    doi = doi.strip()
    if doi != "" and doi[-1] == ".":
        doi = doi[:-1]
    if " " in doi:
        doi = doi.split(" ")[1]  # this handles cases where the DOI string is "DOI: {doi}"
    doi = doi.strip()
    return doi


def get_title_and_citation_from_doi(doi: str) -> str:
    """
    Retrieves the title and citation from a DOI.

    Parameters
    ----------
    doi : str
        The DOI string.

    Returns
    -------
    str
        The title and citation associated with the DOI.
    """

    url = f"https://api.crossref.org/works/{doi}"

    # Send a GET request to the API
    session = setup_retry_session()
    response = session.get(url)

    # If the GET request is successful, the status code will be 200
    if response.status_code == 200:
        # Get the response data
        data = response.json()

        # Get the title and citation count from the data
        try:
            title = data["message"]["title"][0]
            citation = format_citation_crossref(data["message"])
        except Exception:
            try:
                title = data["message"]["items"][0]["title"][0]
                citation = format_citation_crossref(data["message"]["items"][0])
            except Exception:
                return doi
        return f"{title}\n\n - {citation}"
    else:
        return doi


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
    if "family" in first_author:
        author_str = f"{first_author['family']}, {first_author['given']} et al."
    else:
        author_str = f"{first_author['name']} et al."

    journal = " " + message["journal"] if message["journal"] else ""
    year = f"{message['published_year']}"

    return f"{author_str} ({year}){journal}"


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
    if "family" in first_author:
        author_str = f"{first_author['family']}, {first_author['given']} et al."
    else:
        author_str = f"{first_author['name']} et al."

    journal = " " + message["container-title"][0] if len(message["container-title"]) else ""
    year = message["created"]["date-parts"][0][0]

    return f"{author_str} ({year}){journal}"
