# Tests documentation

## How to run E2E tests locally

### What

Running E2E tests locally means the tests are run in your local machine against the web app that's either running locally as well (localhost, docker) OR a deployed environment (dev, staging, prod)

### When

There are a few scenarios you might want to run tests locally:

#### Against local app (`make local-init` or `npm run dev`)

1. You are preparing a new PR and you want to make sure all existing tests pass, and/or you want to add new tests to cover your new feature or bug fix

2. A PR E2E test failed, and you want to see if the test failure is environment dependent.

   E.g., there are two ways to run your FE app locally, either in a Docker container through make commands (`make local-init`) or straight up using `npm run dev` in the `/frontend` directory. PR E2E test runs against local app in a Docker container, so if you run local tests against `npm run dev` started local app and everything is passing, but failing against `make` started Docker container local app, then likely the issue is related to the FE app inside a Docker container.

3. A GHA E2E test against a deployed environment (dev, staging, prod) failed, and you want to see if the test failure is environment dependent.

   E.g., when a test fails against a deployed environment, it's helpful to run the same test locally against your local app to see if it's passing. Because if it is passing, the test failure could be environment dependent, such as env dependent data, data size, AWS machine size, third party rate limiting and/or outage

#### Against a deployed app (dev, staging, prod, rdev)

1. A GHA E2E test against a deployed environment failed, and you want to see if it's just flaky.

   E.g., since tests run a lot slower in GHA vs. locally due to several reasons, such as not as powerful machine, long GHA test queue, and that GHA can't run tests in parallel in the same machine (not many CPUs to utilize), it's sometimes just a lot faster to run that one failed test locally against the same deployed environment

2. A GHA E2E test against a deployed environment failed, and you want to debug the test code.

   E.g., sometimes a test failure is only reproducible in a deployed environment, in such case, you have no choice but to debug the test against the deployed environment

### Where

All the E2E test commands can be found in `frontend/Makefile` and `frontend/package.json`. The `frontend/Makefile` commands are wrappers of `frontend/package.json` commands, so you can use either.

### How

1. Set up your local FE environment before running any `npm` commands:

   1. In `frontend/` directory, run `npm i`

1. local -> local `npm run dev`

   1. Make sure your local app is connected to a working BE API service.

      1. If your `frontend/src/configs/configs.js` is pointing to `API_URL: "http://backend.corporanet.local:5000"`, you need BE running in Docker containers

      1. If your `frontend/src/configs/configs.js` is pointing to `API_URL: "https://api.cellxgene.dev.single-cell.czi.technology"`, or any other deployed BE environments, make sure they are working

   1. Make sure you have your local app running already on http://localhost:3000. If not, in `frontend/` directory, run `npm run dev`

   1. In `frontend/` directory, run `npm run e2e`

1. local -> local container `make local-init`

   1. Make sure you have your local app running already on http://frontend.corporanet.local:3000. If not, in the root directory of the repo, run `make local-init` to start all containers
   1. In `frontend/` directory, run `npm run e2e`

1. local -> dev

   1. Manually check https://cellxgene.dev.single-cell.czi.technology is working
   1. Ask the team what the `TEST_ACCOUNT_PASS` is
   1. In `frontend/` directory, run `TEST_ACCOUNT_PASS=PUT_PASSWORD_HERE npm run e2e-dev`

1. local -> staging

   1. Manually check https://cellxgene.staging.single-cell.czi.technology is working
   1. Ask the team what the `TEST_ACCOUNT_PASS` is
   1. In `frontend/` directory, run `TEST_ACCOUNT_PASS=PUT_PASSWORD_HERE npm run e2e-staging`

1. local -> prod
   1. Manually check https://cellxgene.cziscience.com is working
   1. In `frontend/` directory, run `npm run e2e-prod`
   - NOTE: we don't run logged in tests in prod, since we don't want to accidentally
     add test data to prod database

## How to debug Github action failed E2E tests

1. The following steps will use [this GHA page](https://github.com/chanzuckerberg/single-cell-data-portal/actions/runs/3276702818) as an example
1. Go to the failed Github Action summary page
   ![summary section](https://user-images.githubusercontent.com/6309723/196766311-de738283-b0c0-4de3-bdb0-1b9df1396818.png)
1. Scroll all the way to the `Annotations` section to see what tests failed
   ![annotations section](https://user-images.githubusercontent.com/6309723/196766314-be03a71f-0427-4bf5-89eb-91f6d090a338.png)
1. Now scroll to the bottom to find `Artifacts` section to download screenshots, videos, and network records and responses
   ![artifacts section](https://user-images.githubusercontent.com/6309723/196766315-cb613ab3-7902-4560-bc15-2314d361721c.png)
1. Unzip `test-results`, and you'll find directories of failed tests and their corresponding test artifacts
   1. `test-failed-*.png` shows the screenshot of the app state before the test exists
   1. `video.webm` shows the whole test session, so you can see what the test actually did
   1. `trace.zip` use [Playwright Trace Viewer](https://playwright.dev/docs/trace-viewer) to inspect what the test did play by play along with the network responses at any given time.
      1. In the test directory, run `npx playwright show-trace trace.zip`
         ![test results](https://user-images.githubusercontent.com/6309723/196766317-95afbb3c-2890-42af-a42e-0b1f7702d73c.png)
         ![Trace Viewer](https://user-images.githubusercontent.com/6309723/196768996-899f1086-13e8-4c3d-804a-eb050bfa6c71.gif)
1. If the information above is not enough to help you pinpoint the root cause of the test failure, we can use [Playwright Inspector](https://playwright.dev/docs/debug) to rerun the [test locally](#how-to-run-e2e-tests-locally) and debug the test steps live

   1. Find the test you want to debug.
      The `test-results` directory name is part of the test name, so you can use that to find the actual test case.

      For example, `test-results` directory name `features-collection-collection-Collection-Depl-5df78-te-a-collection-with-a-DOI-in-an-invalid-format-chromium` points us to `frontend/tests/features/collection/collection.test.ts` and test case `doesn't create a collection with a DOI in an invalid format`

   1. Once you find the test, add `.only` to the `test()` call, so Playwright will only run that test in the file
      ![Add .only](https://user-images.githubusercontent.com/6309723/196773185-8553273c-0cfd-4861-ac57-e2e21abfb76a.png)

   1. Use `await page.pause()` anywhere in the test code to pause Playwright execution

      ![page.pause()](https://user-images.githubusercontent.com/6309723/196776657-ac957fd6-935f-44a1-9188-490e85ad3e74.png)

   1. Since we're debugging Dev environment, we need to use the `local -> dev` command as shown in the [How](#how) section above.

   1. In `single-cell-data-portal/frontend/`, run `TEST_ACCOUNT_PASS=PUT_PASSWORD_HERE npm run e2e-dev -- /YOUR_DIRECTORY_PATH/single-cell-data-portal/frontend/tests/features/collection/collection.test.ts -- --debug`.

      NOTE: `PUT_PASSWORD_HERE` needs to be replaced with the actual test account Auth0 password, and `YOUR_DIRECTORY_PATH` needs to be where your `single-cell-data-portal` directory lives

   1. Playwright will now spin up a browser window and a debugging console for you to inspect the test live!

      1. You will find the play button in the console, so you can play/pause the test as you see fit
      1. For more information, visit [Playwright Inspector](https://playwright.dev/docs/debug)

      ![Playwright Inspector](https://user-images.githubusercontent.com/6309723/196775284-81bd4853-b4c6-45c4-827b-22e97b97ee06.png)
