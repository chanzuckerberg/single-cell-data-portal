"""
This is a script to upload cellguide descriptions using the /cellguide/v1/upload endpoint.

It requires a csv file to be present in the scripts directory with the following columns:
Cell Ontology ID, Final version (QC'd), Supporting reference 1, Supporting reference 2, Supporting reference 3, 
Supporting reference 4, Supporting reference 5

The script will convert the csv to json and post each row to the endpoint.
Step 1 - Login: python upload_cellguide_descriptions.py --deployment staging login --api_key <api_key>
Step 2 - Upload: python upload_cellguide_descriptions.py --deployment staging upload --path cl_descriptions.csv
"""

import csv
import json
import os

import click
import requests


@click.group()
@click.option("--deployment", default="test", show_default=True, help="The name of the deployment to target")
@click.pass_context
def cli(ctx, deployment):
    """
    Set deployment stage and api base url
    """
    os.environ["DEPLOYMENT_STAGE"] = deployment
    api_base = {
        "rdev": "https://pr-6426-backend.rdev.single-cell.czi.technology",
        "test" or "dev": "https://api.cellxgene.dev.single-cell.czi.technology",
        "staging": "https://api.cellxgene.staging.single-cell.czi.technology",
        "production": "https://api.cellxgene.cziscience.com",
    }
    ctx.obj = {"api_base": api_base[deployment]}


@cli.command()
@click.pass_context
@click.option("--api_key", default="", help="curator api key")
def login(ctx, api_key):
    """
    Login using the curator api key
    """
    login_url = f"{ctx.obj['api_base']}/curation/v1/auth/token"
    headers = {"Content-Type": "application/json", "x-api-key": api_key}
    response = requests.post(login_url, headers=headers)
    os.environ["API_TOKEN"] = json.loads(response.text)["access_token"]
    click.echo(f"API_TOKEN has been set to {os.environ['API_TOKEN']}!")


@cli.command()
@click.option("--path", default=".", help="path to csv")
@click.pass_context
def upload(ctx, path):
    """
    Upload cellguide descriptions
    """
    api_token = os.environ.get("API_TOKEN")
    if not api_token:
        click.echo("API_TOKEN is not set. Please run login first.")
        return
    click.echo(api_token)

    url = f"{ctx.obj['api_base']}/cellguide/v1/upload"
    headers = {"Content-Type": "application/json", "Authorization": f"Bearer {api_token}"}

    with open(path, "r") as file:
        reader = csv.DictReader(file)
        for row in reader:
            json_row = convert_to_json(row)
            response = requests.post(url, headers=headers, data=json_row)
            click.echo(response.text)


def convert_to_json(csv_row):
    """
    Convert csv row to json
    """
    json_row = {
        "cell_ontology_term_id": csv_row["Cell Ontology ID"].replace(":", "_"),
        "description": csv_row["Final version (QC'd)"],
        "references": format_references(csv_row),
    }
    return json.dumps(json_row)


def format_references(csv_row):
    """
    Format references list
    """
    references = []
    for i in range(1, 6):
        if csv_row[f"Supporting reference {i}"]:
            url = csv_row[f"Supporting reference {i}"]
            if not (url.startswith("http://") or url.startswith("https://")):
                url = "https://www.doi.org/" + url
            references.append(url)
    return references


if __name__ == "__main__":
    cli()
