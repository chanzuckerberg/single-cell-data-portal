# Workflow Taxonomy

## Where and what tests are run

---

Deployments do not run in GitHub Actions. Follow [the manual deployment runbook](../../docs/manual-deployment.md).

### Workflow: `lint-pr.yml`

#### Ran on:

- `pull_request_target`

  - opening, editing, and syncing a pull request

  #### Summary: Runs a GitHub action that lints the PR commit message according to conventional commit standards

### Workflow: `push-tests.yml`

#### Ran on:

- `push` to `main`, `staging`, `prod`
- `pull_request` targeting `"*"`

#### Jobs:

- ##### `lint`:

  - Inside backend container: `make lint`
    - https://github.com/chanzuckerberg/single-cell-data-portal/blob/6a423183c255737d2a44e40447a91d0ece041a41/Makefile#L19-L20
  - Inside frontend container: (frontend/Makefile)`make lint`
    - https://github.com/chanzuckerberg/single-cell-data-portal/blob/6a423183c255737d2a44e40447a91d0ece041a41/frontend/Makefile#L18-L19
      - https://github.com/chanzuckerberg/single-cell-data-portal/blob/6a423183c255737d2a44e40447a91d0ece041a41/frontend/package.json#L111
      - https://github.com/chanzuckerberg/single-cell-data-portal/blob/6a423183c255737d2a44e40447a91d0ece041a41/frontend/package.json#L110
        - `node_modules/.bin/next lint` -> runs eslint with config [next.js docs](https://nextjs.org/docs/basic-features/eslint)
        - `node_modules/.bin/stylelint --fix '**/*.{js,ts,tsx,css}'` -> runs stylelint with config
        - https://github.com/chanzuckerberg/single-cell-data-portal/blob/6a423183c255737d2a44e40447a91d0ece041a41/frontend/package.json#L109

- ##### `e2e-test`:

  - Installs dependencies and runs local frontend server pointing at dev env BE API:
    - https://github.com/chanzuckerberg/single-cell-data-portal/blob/6a423183c255737d2a44e40447a91d0ece041a41/.github/workflows/push-tests.yml#L80-L83
  - Runs e2e tests:
    - `npm run e2e` ->`playwright test`
      - runs all tests in `frontend/tests` using playwright config, no TEST_ENV provided defaults to `local`
  - Pushes images to RCS

- ##### `build-extra-images`: TODO

- ##### `push-prod-images`: TODO

- ##### `backend-unit-test`:

  - Checks if containers need to be rebuilt based on diffs on docker + requirements files
  - Runs tests in docker-compose: https://github.com/chanzuckerberg/single-cell-data-portal/blob/6a423183c255737d2a44e40447a91d0ece041a41/.github/workflows/push-tests.yml#L214-L216
    - https://github.com/chanzuckerberg/single-cell-data-portal/blob/6a423183c255737d2a44e40447a91d0ece041a41/Makefile#L106-L107
    - https://github.com/chanzuckerberg/single-cell-data-portal/blob/6a423183c255737d2a44e40447a91d0ece041a41/Makefile#L182-L185
      - https://github.com/chanzuckerberg/single-cell-data-portal/blob/6a423183c255737d2a44e40447a91d0ece041a41/Makefile#L27-L30

- ##### `processing-unit-test`:

  - Checks if containers need to be rebuilt based on diffs on docker + requirements files
  - Runs tests in docker-compose: https://github.com/chanzuckerberg/single-cell-data-portal/blob/6a423183c255737d2a44e40447a91d0ece041a41/.github/workflows/push-tests.yml#L260-L262
    - https://github.com/chanzuckerberg/single-cell-data-portal/blob/6a423183c255737d2a44e40447a91d0ece041a41/Makefile#L106-L107
    - https://github.com/chanzuckerberg/single-cell-data-portal/blob/6a423183c255737d2a44e40447a91d0ece041a41/Makefile#L189-L192
      - https://github.com/chanzuckerberg/single-cell-data-portal/blob/6a423183c255737d2a44e40447a91d0ece041a41/Makefile#L33-L36

- ##### `wmg-processing-unit-test`:

  - Checks if containers need to be rebuilt based on diffs on docker + requirements files
  - Runs tests in docker-compose: https://github.com/chanzuckerberg/single-cell-data-portal/blob/6a423183c255737d2a44e40447a91d0ece041a41/.github/workflows/push-tests.yml#L306-L307
    - `https://github.com/chanzuckerberg/single-cell-data-portal/blob/6a423183c255737d2a44e40447a91d0ece041a41/Makefile#L196-L199
      - https://github.com/chanzuckerberg/single-cell-data-portal/blob/6a423183c255737d2a44e40447a91d0ece041a41/Makefile#L39-L42

- ##### `push-image`: TODO
- ##### `create_deployment`: TODO

### Workflow: `scale-test.yml`: TODO
