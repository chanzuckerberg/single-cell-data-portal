import csv
import json
import os
from time import sleep
import click
import requests

@click.group()
@click.option('--deployment', default="test", show_default=True, help='The name of the deployment to target')
@click.pass_context
def cli(ctx, deployment):
    """set deployment state"""
    os.environ["DEPLOYMENT_STAGE"] = deployment
    api_base = {
        "dev": "https://api.cellxgene.dev.single-cell.czi.technology",
        "staging": "https://api.cellxgene.staging.single-cell.czi.technology",
        "production": "https://api.cellxgene.cziscience.com",
    }
    ctx.obj = {'api_base': api_base[deployment]}

@cli.command()
@click.pass_context
@click.option('--api_key', default='', help='api key')
def login(ctx, api_key):
    """login using the curator api key"""

    login_url = f"{ctx.obj['api_base']}/curation/v1/auth/token"
    headers = {"Content-Type": "application/json", "x-api-key": api_key}
    response = requests.post(login_url, headers=headers)
    
    os.environ['API_TOKEN'] = json.loads(response.text)['access_token']
    click.echo(f"API_TOKEN has been set to {os.environ['API_TOKEN']}!")

@cli.command()
@click.option('--path', default='.', help='path to csv')
@click.pass_context
def upload(ctx, path):
    """Upload cellguide descriptions"""
    api_token = os.environ.get('API_TOKEN')
    if not api_token:
        click.echo('API_TOKEN is not set. Please run `cellxgene login` first.')
        return
    click.echo(api_token)
    url = f"{ctx.obj['api_base']}/cellguide/v1/upload"
    headers = {"Content-Type": "application/json", "Authorization": f"Bearer {api_token}"}
    with open(path, 'r') as file:
        reader = csv.DictReader(file)
        for row in reader:
            json_row = convert_to_json(row)
            print(json_row)
            sleep(10)
            print(url)
            response = requests.post(url, headers=headers, data=json_row)
            click.echo(response.text)
        # response = requests.post(url, headers=headers, data=file.read())
        # click.echo(response.text)
    
def convert_to_json(csv_row):
    """Convert csv row to json"""
    json_row = {
            'cell_ontology_term_id': csv_row['Cell Ontology ID'].replace(":", "_"),
            'description': csv_row['Final version (QC\'d)'],
            'references': format_references(csv_row),
    }
    return json.dumps(json_row)

def format_references(csv_row):
    """Format references list"""
    references = []
    for i in range(1, 6):
        if csv_row[f'Supporting reference {i}']:
            url = csv_row[f'Supporting reference {i}']
            if not (url.startswith('http://') or url.startswith('https://')):
                url = 'https://www.doi.org/' + url
            references.append(url)
    return references

if __name__ == '__main__':
    cli()
