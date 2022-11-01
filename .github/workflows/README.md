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
    - ```bash
            if [ "${DEPLOYMENT_STAGE}" == "stage" ]; then
            export DEPLOYMENT_STAGE=staging
          fi
          if [ "${DEPLOYMENT_STAGE}" != "prod" ]; then
            pip3 install -r scripts/smoke_tests/requirements.txt
            python3 -m scripts.smoke_tests.setup
          fi
          cd frontend
          npm ci
          npx playwright install --with-deps
          cp src/configs/${DEPLOYMENT_STAGE}.js src/configs/configs.js
          DEBUG=pw:api npm run e2e-${DEPLOYMENT_STAGE}
      ```
- ##### `functional-test`:
  - If environment is not `prod` start up backend container and run `make local-functional-test`
  - `make local-functional-test` -> `python3 -m unittest discover --start-directory tests/functional --top-level-directory . --verbose`
- ##### `performance-test`:
  - If environment is `prod` run `make prod-performance-test`
  - `make prod-performance-test` -> `python3 -m unittest discover --start-directory tests/performance --top-level-directory . --verbose`

### Workflow: `lint-pr.yml`

#### Ran on:

- `pull_request_target`

  - opening, editing, and syncing a pull request

  #### Summary: Runs a GitHub action that lints the PR commit message

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

  - Inside backend container: `make lint` -> `flake8 backend tests`
  - Inside frontend container: (frontend/Makefile)`make lint` -> `npm run prettier-check & npm run lint`
    - `npm run prettier-check` -> `prettier --check .`
    - `npm run lint` -> `concurrently \"node_modules/.bin/next lint\" \"node_modules/.bin/stylelint --fix '**/*.{js,ts,tsx,css}'\" \"npm run type-check\"`
      - `node_modules/.bin/next lint` -> runs eslint with config [next.js docs](https://nextjs.org/docs/basic-features/eslint)
      - `node_modules/.bin/stylelint --fix '**/*.{js,ts,tsx,css}'` -> runs stylelint with config
      - `npm run type-check` -> `tsc --noEmit`

- ##### `e2e-test`:

  - Installs dependencies and runs local frontend server pointing at dev env BE API:
    - `npm ci &&npx playwright install --with-deps && cp src/configs/dev.js src/configs/configs.js && npm run dev&`
  - Runs e2e tests:
    - `npm run e2e` ->`playwright test`
      - runs all tests in `frontend/tests` using playwright config, no TEST_ENV provided defaults to `local`
  - Pushes images to RCS

- ##### `build-extra-images`: TODO

- ##### `push-prod-images`: TODO

- ##### `backend-unit-test`:

  - Checks if containers need to be rebuilt based on diffs on docker + requirements files
  - Runs tests in docker-compose: `make local-init-test-data & make all-local-unit-test-backend`
    - `make local-init-test-data` -> `docker-compose $(COMPOSE_OPTS) run --rm -T backend /bin/bash -c "pip3 install awscli && cd /single-cell-data-portal && scripts/setup_dev_data.sh"`
    - `make all-local-unit-test-backend` -> (inside backend container) `make container-unittest;`
      - `make container-unittest` -> `DEPLOYMENT_STAGE=test PYTHONWARNINGS=ignore:ResourceWarning python3 \ -m unittest discover --start-directory tests/unit/backend --top-level-directory . --verbose;`

- ##### `processing-unit-test`:

  - Checks if containers need to be rebuilt based on diffs on docker + requirements files
  - Runs tests in docker-compose
    - `make local-init-test-data` -> `docker-compose $(COMPOSE_OPTS) run --rm -T backend /bin/bash -c "pip3 install awscli && cd /single-cell-data-portal && scripts/setup_dev_data.sh"`
    - `make local-unit-test-processing` -> (inside processing container) `make processing-unittest;`
      - `make processing-unittest` -> `DEPLOYMENT_STAGE=test PYTHONWARNINGS=ignore:ResourceWarning python3 \ -m unittest discover --start-directory tests/unit/processing_container --top-level-directory . --verbose;`

- ##### `wmg-processing-unit-test`:

  - Checks if containers need to be rebuilt based on diffs on docker + requirements files
  - Runs tests in docker-compose:
    - `make local-unit-test-wmg-processing` -> (inside wmg_processing container)`make wmg-processing-unittest`
      - `make wmg-processing-unittest` -> `DEPLOYMENT_STAGE=test PYTHONWARNINGS=ignore:ResourceWarning python3 \ -m unittest discover --start-directory tests/unit/wmg_processing --top-level-directory . --verbose;`

- ##### `push-image`: TODO
- ##### `create_deployment`: TODO

### Workflow: `scale-test.yml`: TODO
