# Local Development Environment

## Development quickstart
1. [install docker](https://docs.docker.com/get-docker/). If brew is installed run `brew install docker`.
1. From the root of this repository, run `make dev-init` to build and run the dev environment. The first build takes awhile, but subsequent runs will use cached artifacts.
1. Visit [http://localhost:5000](http://localhost:5000) to view the backend, and [http://localhost:8000](http://localhost:8000) for the frontend.
1. **Open the source code and start editing!**
  - Modify code in the `frontend/src` directory, save your changes and the browser will update in real time.
  - Modify code in the `backend` directory, and the data portal api will reload automatically.

### Containers managed by the dev environment
The data portal dev environment is a set of containers defined in [docker-compose.yml](docker-compose.yml). The [backend docker image](backend/Dockerfile) and [frontend docker image](frontend/Dockerfile) are built locally. Update any of these files as necessary and run `make dev-sync` to sync your dev environment with these configs.
![Dev Environment Containers](docs/dev_env.jpg)

### Updating frontend/backend dependencies
Both the data portal frontend and backend services will automatically reload when their source code is modified, but they won't automatically rebuild when their dependencies (such as npm or pip package lists) change.

To update the dev env to reflect changes to [frontend/package.json](frontend/package.json) or [backend/chalice/api_server/requirements.txt](backend/chalice/api_server/requirements.txt), run `make dev-sync`

### Update Dev Data
The dev environment is initialized with AWS Secrets/S3 data in the [scripts/setup_dev_data.sh](scripts/setup_dev_data.sh) script, as well as DB data from [tests/unit/backend/corpora/fixtures/database/__init__.py](tests/unit/backend/corpora/fixtures/database/__init__.py). To add more data, modify these scripts and run `make dev-init` to reload the dev environment's data stores.

### Make targets for managing dev:

| Command                 | Description                                                                          | Notes                                                                                                |
| ----------------------- | ------------------------------------------------------------------------------------ | ---------------------------------------------------------------------------------------------------- |
| `make dev-init`         | Launch a new local dev env and populate it with test data.                           |                                                          |
| `make dev-start`        | Start a local dev environment that's been stopped.                                   |                                                          |
| `make dev-stop`         | Stop the local dev environment.                                                      |                                                          |
| `make dev-logs`         | Tail the logs of the dev env containers.                                             | Run `make dev-logs CONTAINER=backend` to tail the logs of a specific container. Dev containers are: backend, frontend, localstack, database, oidc |
| `make dev-shell CONTAINER=frontend`  | Open a command shell in one of the dev containers                       | Dev containers are: backend, frontend, localstack, database, oidc |
| `make dev-status`       | Show the status of the containers in the dev environment.                            |                                                          |
| `make dev-clean`        | Remove everything related to the local dev environment (including db data!)          |                                                          |
| `make dev-sync`         | Re-sync the dev-environment state after modifying library deps or docker configs     |                                                          |

### Make targets for running tests in dev
| Command                 | Description                                                                          | Notes                                                                                                |
| ----------------------- | ------------------------------------------------------------------------------------ | ---------------------------------------------------------------------------------------------------- |
| `make dev-unit-test`    | Run backend tests in the dev environment                                             |                                                          |
| `make dev-functional-test` | Run functional tests in the dev environment                                       |                                                          |
| `make dev-smoke-test`   | Run frontend/e2e tests in the dev environment                                        |                                                          |

### External dependencies
The dev environment has no network dependencies, but it launches some extra containers to mock external dependencies:
 - [LocalStack](https://github.com/localstack/localstack) to mock AWS
 - [OIDC server mock](https://github.com/Soluto/oidc-server-mock) in place of Auth0.
 - [postgres](https://hub.docker.com/_/postgres) in place of RDS.

#### TLS Certificate for mock authentication service
Due to browser security considerations, we must run the mock authentication
service using a self-signed certificate. The dev-init and dev-clean make targets
handle managing a keypair/certificate for each dev env and installing it in the
OSX system keychain.

Details: OIDC requires setting a token, and requires the cookie storing that
token to be stored with samesite=None to work properly. Recent versions of
browsers such as Chrome intentionally only allow samesite=None if the connection
is over a secure network connection i.e. TLS. Thus we need to run even a local
development auth service behind a certificate. We bundle a pre-generated
self-signed cert in for convenience.
