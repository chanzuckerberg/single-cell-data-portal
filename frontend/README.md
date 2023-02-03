# Data Portal Front End

## Development

See [DEV_ENV.md](DEV_ENV.md) for local dev quickstart

NOTE: If you want to `npm i` in your local machine, please make sure to install [`nvm`](https://github.com/nvm-sh/nvm), since we use `.nvmrc` to specify which node version to install the dependencies

## Environment Variables

The environment variables for the web application. The variables are stored in /frontend/configs/\*. E.g., `frontend/configs/local.js`

For local development, please copy `local.js` to a new file named `configs.js`
in the same directory (`frontend/configs/configs.js`)

WARNING: Do not store sensitive data in the environment variables.

| Name            | Description                                          |
| --------------- | ---------------------------------------------------- |
| AUTH0_DOMAIN    | The hosted Auth0 domain used for Authentication      |
| AUTH0_CLIENT_ID | The client id of the Auth0 application for this site |
| AUDIENCE        | The domain of the corpora api                        |
| API_URL         | The URL to the corpora api                           |

## Deployment

1. Ensure your `awscli` is configured with the
   [required credentials and profiles](../docs/awscli.md).
   Set the appropriate `AWS_PROFILE`.

   ```shell
   export AWS_PROFILE=single-cell-dev
   ```

1. **Specify deployment.**

   Set the `DEPLOYMENT_STAGE` environment variable to a valid deployed environment: `test`, `staging`

   ```shell
   export DEPLOYMENT_STAGE=test
   ```

1. **Deploy.**

   Files are deployed to a publicly accessible bucket. Do not include sensitive data in the deployed files.

   ```shell
   make deploy
   ```

## e2e Tests

See [here](tests/README.md) for details
