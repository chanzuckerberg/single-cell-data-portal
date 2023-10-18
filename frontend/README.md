# Data Portal Front End

## Mise en place

1. Install recommended VSCode extensions ([source](https://dev.to/askrishnapravin/recommend-vs-code-extensions-to-your-future-teammates-4gkb))
   1. Open the repo in VSCode
   1. You should get a notification to install the recommended extensions
      ![Capture-2023-09-21-102235](https://github.com/chanzuckerberg/single-cell-data-portal/assets/6309723/07eaff7b-420e-457e-942d-6f02c1609660)
   1. If not, go to the VSCode extension tab in the editor, and type `@recommended` in the searchbox to see the list.
      ![search for recommended extensions](https://github.com/chanzuckerberg/single-cell-data-portal/assets/6309723/a765a536-b54d-4ca1-9072-68f43bf9c09b)
   1. See [here](https://docs.google.com/document/d/1qveZszisGdH6FvP5XI6y5re93A0ZI9aX_k45veiyWMY/edit#bookmark=id.cdwlzmjo98io) for why the extensions are recommended

## Development

The following steps will start a FE server that connects to `dev` API. (See [useful tips](#useful-tips) section to connect to a different API)

1. Install [`nvm`](https://github.com/nvm-sh/nvm)
   - Example: `brew install nvm`
1. Check `.nvmrc` to see which version of node to download.
   - Example: `nvm install 16.14.2 && nvm use 16.14.2`
1. Install npm packages and playwright browsers
   - Example: `npm i && npx playwright install`
1. Copy configs file `frontend/src/configs/local.js` to `frontend/src/configs/configs.js`
   - Example: `cp ./src/configs/local.js ./src/configs/configs.js`
1. Start FE server
   - Example: `npm run dev`
1. Navigate to `http://localhost:3000/`

### Option: With Docker

See [DEV_ENV.md](../DEV_ENV.md) for local dev quick start using Docker containers

### Useful tips

1. Test FE app production build locally: `npm run build && npm run serve`
1. Connect FE app to a different deployed env API
   1. Go to `frontend/src/configs/configs.js`, uncomment the env API you want to use

## Environment Variables

The environment variables for the web application. The variables are stored in /frontend/configs/\*. E.g., `frontend/configs/local.js`

For local development, please copy `local.js` to a new file named `configs.js`
in the same directory (`frontend/src/configs/configs.js`)

WARNING: Do not store sensitive data in the environment variables.

| Name            | Description                                          |
| --------------- | ---------------------------------------------------- |
| AUTH0_DOMAIN    | The hosted Auth0 domain used for Authentication      |
| AUTH0_CLIENT_ID | The client id of the Auth0 application for this site |
| AUDIENCE        | The domain of the corpora api                        |
| API_URL         | The URL to the corpora api                           |

## Deployment

Steps are described in [On Call Check List](https://docs.google.com/document/d/1G2NTjXTJJeHyhqvnyzYmcO0Um24Ph0dCLUyMIWZvLfg/edit#)

## e2e Tests

See [tests documentation](tests/README.md) for details
