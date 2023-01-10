# Workflow Taxonomy

## Where and what tests are run

---

### Workflow: `deploy-happy-stack.yml`

#### Ran on:

- `deployment`: (deployments are triggered by `push-tests`)

#### Jobs:

- ##### `upgrade`: TODO
- ##### `e2e-test`:
  - Runs playwright against deployed frontend app uploading any artifacts
    - https://github.com/chanzuckerberg/single-cell-data-portal/blob/6a423183c255737d2a44e40447a91d0ece041a41/.github/workflows/deploy-happy-stack.yml#L133-L144
- ##### `functional-test`:
  - If environment is not `prod` start up backend container and run `make local-functional-test`
  - https://github.com/chanzuckerberg/single-cell-data-portal/blob/6a423183c255737d2a44e40447a91d0ece041a41/Makefile#L215-L223
    - https://github.com/chanzuckerberg/single-cell-data-portal/blob/6a423183c255737d2a44e40447a91d0ece041a41/Makefile#L49-L51
- ##### `performance-test`:
  - If environment is `prod` run `make prod-performance-test`
  - https://github.com/chanzuckerberg/single-cell-data-portal/blob/6a423183c255737d2a44e40447a91d0ece041a41/Makefile#L54-L55

### Workflow: `push-processing-base.yml`

TODO

### Workflow: `push-rdev.yml`

TODO

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
