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

## Smoke Test

See [here](tests/README.md) for details

## 🧐 What's inside

A quick look at the top-level files and directories you'll see in a Gatsby project.

.
├── node_modules
├── src
├── .gitignore
├── .prettierrc
├── gatsby-browser.js
├── gatsby-config.js
├── gatsby-node.js
├── gatsby-ssr.js
├── LICENSE
├── package-lock.json
├── package.json
└── README.md

1. **`/node_modules`**: This directory contains all of the modules of code that your project depends on (npm packages) are automatically installed.

2. **`/src`**: This directory will contain all of the code related to what you will see on the front-end of your site (what you see in the browser) such as your site header or a page template. `src` is a convention for “source code”.

3. **`.gitignore`**: This file tells git which files it should not track / not maintain a version history for.

4. **`.prettierrc`**: This is a configuration file for [Prettier](https://prettier.io/). Prettier is a tool to help keep the formatting of your code consistent.

5. **`gatsby-browser.js`**: This file is where Gatsby expects to find any usage of the [Gatsby browser APIs](https://www.gatsbyjs.org/docs/browser-apis/) (if any). These allow customization/extension of default Gatsby settings affecting the browser.

6. **`gatsby-config.js`**: This is the main configuration file for a Gatsby site. This is where you can specify information about your site (metadata) like the site title and description, which Gatsby plugins you’d like to include, etc. (Check out the [config docs](https://www.gatsbyjs.org/docs/gatsby-config/) for more detail).

7. **`gatsby-node.js`**: This file is where Gatsby expects to find any usage of the [Gatsby Node APIs](https://www.gatsbyjs.org/docs/node-apis/) (if any). These allow customization/extension of default Gatsby settings affecting pieces of the site build process.

8. **`gatsby-ssr.js`**: This file is where Gatsby expects to find any usage of the [Gatsby server-side rendering APIs](https://www.gatsbyjs.org/docs/ssr-apis/) (if any). These allow customization of default Gatsby settings affecting server-side rendering.

9. **`LICENSE`**: Gatsby is licensed under the MIT license.

10. **`package-lock.json`** (See `package.json` below, first). This is an automatically generated file based on the exact versions of your npm dependencies that were installed for your project. **(You won’t change this file directly).**

11. **`package.json`**: A manifest file for Node.js projects, which includes things like metadata (the project’s name, author, etc). This manifest is how npm knows which packages to install for your project.

12. **`README.md`**: A text file containing useful reference information about your project.

## 🎓 Learning Gatsby

Looking for more guidance? Full documentation for Gatsby lives [on the website](https://www.gatsbyjs.org/). Here are some places to start:

- **For most developers, we recommend starting with our [in-depth tutorial for creating a site with Gatsby](https://www.gatsbyjs.org/tutorial/).** It starts with zero assumptions about your level of ability and walks through every step of the process.

- **To dive straight into code samples, head [to our documentation](https://www.gatsbyjs.org/docs/).** In particular, check out the _Guides_, _API Reference_, and _Advanced Tutorials_ sections in the sidebar.
