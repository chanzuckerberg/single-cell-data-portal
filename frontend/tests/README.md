# Tests documentation

## WARNING

1. Do **NOT** import `.tsx` files in your tests, otherwise we'll get weird parsing errors, due to Playwright parser not parsing SVG files at the moment.

1. Use `test` from `import { test } from "tests/common/test"` instead of `import { test } from "@playwright/test"`

## How to run e2e tests from your local machine

See the [#how](#how) section

## Flags

1. `LOGIN` and `SKIP_LOGIN`: Add `LOGIN=false` or `SKIP_LOGIN=true` to your test command if you don't need Playwright to log into Data Portal
   - Example: `LOGIN=false npm run e2e`
1. `HEADLESS` and `HEADFUL`: Add `HEADLESS=false` or `HEADFUL=true` to your test command to launch browser
   - Example: `HEADFUL=true npm run e2e`
1. `RETRY`: Add `RETRY=false` if you don't want to retry your test. This is good for failing fast when you're writing tests
   - Example: `RETRY=false npm run e2e`
1. `USE_COOKIE`: Manually use your own cookie for authenticated tests. Should only be used locally. **The cookie value is set by modifying the `MANUAL_COOKIE` variable in `playwright.config.ts`**
   - Example: `USE_COOKIE=true npm run e2e`
1. `RDEV_TOKEN`: Fetch and use access token to access an rdev BE in local FE. NOTE: This requires manually changing your local FE to hit an rdev BE URL
   - Example: `RDEV_TOKEN=true npm run e2e`

## Cheat Sheet

1. Run individual test file: Instead of executing the all test suite, you can pass a test file to the command to speed up your debug cycle

   - Example: `npm run e2e -- -- TEST_FILE_PATH`

1. `.only()`: Add `.only()` to a `test()` call, so Playwright will only run that test
   ![Add .only](https://user-images.githubusercontent.com/6309723/196773185-8553273c-0cfd-4861-ac57-e2e21abfb76a.png)

1. `--ui`: UI Mode lets you explore, run, and debug tests with a time travel experience complete with watch mode. All test files are loaded into the testing sidebar where you can expand each file and describe block to individually run, view, watch and debug each test.

   - [Source](https://playwright.dev/docs/test-ui-mode#running-tests-in-ui-mode)

   - Example: `npm run e2e -- --ui -- TEST_FILE_PATH`

   ![Image](https://user-images.githubusercontent.com/13063165/234295914-f7ee3d8b-33a7-41b3-bc91-d363baaa7305.png)

1. `--debug`: Debug mode launches Playwright Inspector, which lets you play, pause, or step through each action of your test using the toolbar at the top of the Inspector. You can see the current action highlighted in the test code, and matching elements highlighted in the browser window.

   - [Source](https://playwright.dev/docs/debug#playwright-inspector)

   - Example: `npm run e2e -- --debug -- TEST_FILE_PATH`

   ![Image](https://user-images.githubusercontent.com/13063165/212924587-4b84e5f6-b147-40e9-8c75-d7b9ab6b7ca1.png)

1. `npm run e2e-trace`: Use [Playwright Trace Viewer](https://playwright.dev/docs/trace-viewer) to inspect what the test did play by play along with the network responses at any given time

   - [Source](https://playwright.dev/docs/trace-viewer)

   - Example: `npm run e2e-trace PATH_TO_TRACE.zip`

     ![Trace Viewer](https://user-images.githubusercontent.com/6309723/196768996-899f1086-13e8-4c3d-804a-eb050bfa6c71.gif)

### What

Running e2e tests locally means the tests are run in your local machine against the web app that's either running locally as well (localhost, docker) OR a deployed environment (dev, staging, prod, and rdev)

### When

There are a few scenarios you might want to run tests locally:

#### Against local app (`npm run dev` (recommended) or `make local-init`)

1. You are preparing a new PR and you want to make sure all existing tests pass, and/or you want to add new tests to cover your new feature or bug fix

2. A PR e2e test failed, and you want to see if the test failure is BE environment dependent.

   Since PR e2e tests run against their own PR rdev, the quickest way to test if the test failure is BE environment dependent is to test out your local FE with different BE environments and see if the test fail in one BE environment, but not in another

   For example, the following steps will test local FE + rdev BE and local FE + staging BE:

   local FE + rdev BE:

   1. Go to `frontend/src/configs/configs.js` and change API_URL to your rdev BE URL.
      e.g., `API_URL: "https://pr-6062-backend.rdev.single-cell.czi.technology",`

   1. Restart your local FE server, so the local FE app will connect to the rdev BE

   1. Visit https://localhost:3000/ from your favorite browser and make sure the app works

   1. Run `RDEV_TOKEN=true npm run e2e`
      1. IMPORTANT: `RDEV_TOKEN` flag needs to be set, so the test suite setup will fetch and use rdev token. **Otherwise the browsers that Playwright runs will get CORS error when hitting rdev BE endpoints from your local FE app**. This can be confusing, especially since testing the app from the step above will work with your favorite browser, since your browser likely has cached an access token in your browser cookies, which is NOT available in the browsers Playwright spins up
      1. You'll likely want to add `.only` to specific test(s) instead of running the whole test suite. See [Cheat Sheet](#cheat-sheet) section for details

   local FE + staging BE:

   1. Go to `frontend/src/configs/configs.js` and change API_URL to staging BE URL.
      e.g., `API_URL: "https://api.cellxgene.staging.single-cell.czi.technology",`

   1. Restart your local FE server, so the local FE app will connect to the staging BE

   1. Visit https://localhost:3000/ from your favorite browser and make sure the app works

   1. Run `npm run e2e`
      1. NOTE: We don't need to pass `RDEV_TOKEN=true` here, since we're hitting staging BE

3. A GHA e2e test against a deployed environment (dev, staging, prod, rdev) failed, and you want to see if the test failure is environment dependent.

   E.g., when a test fails against a deployed environment, it's helpful to run the same test locally against your local app to see if it's passing. Because if it is passing, the test failure could be environment dependent, such as env dependent data, data size, AWS machine size, third party rate limiting and/or outage

#### Against a deployed app (dev, staging, prod, rdev)

1. A GHA e2e test against a deployed environment failed, and you want to see if it's just flaky.

   E.g., since tests run a lot slower in GHA vs. locally due to several reasons, such as not as powerful machine, long GHA test queue, and that GHA can't run tests in parallel in the same machine (not many CPUs to utilize), it's sometimes just a lot faster to run that one failed test locally against the same deployed environment

2. A GHA e2e test against a deployed environment failed, and you want to debug the test code.

   E.g., sometimes a test failure is only reproducible in a deployed environment, in such case, you have no choice but to debug the test against the deployed environment

### Where

All the e2e test commands can be found in `frontend/Makefile` and `frontend/package.json`. The `frontend/Makefile` commands are wrappers of `frontend/package.json` commands, so you can use either.

### How

#### Set up

Before running any tests, start with [mise en place](../README.md#mise-en-place)

#### Options

1. local -> local (app started with `npm run dev`)

   1. Make sure you have your local app running already on https://localhost:3000. If not, in `frontend/` directory, run `npm run dev`

   1. In `frontend/` directory, run `npm run e2e`

   - NOTE: `SKIP_LOGIN=true npm run e2e` if login is not required for the tests.

1. local -> local container (app started with `make local-init`)

   1. Make sure you have your local app running already on https://frontend.corporanet.local:3000. If not, in the root directory of the repo, run `make local-init` to start all containers

   1. In `frontend/` directory, run `npm run e2e`

   - NOTE: `SKIP_LOGIN=true npm run e2e` if login is not required for the tests.

1. local -> dev

   1. Manually check https://cellxgene.dev.single-cell.czi.technology is working

   1. In `frontend/` directory, run `TEST_ACCOUNT_PASS=PUT_PASSWORD_HERE npm run e2e-dev`

   - NOTE: Replace `PUT_PASSWORD_HERE` with `corpora/backend/dev/auth0-secret.test_account_password` in AWS Secret Manager
   - NOTE: To run in one specific browser, run `npm run e2e-dev -- --project chromium` (chromium|firefox|edge)

1. local -> staging

   1. Manually check https://cellxgene.staging.single-cell.czi.technology is working
   1. In `frontend/` directory, run `TEST_ACCOUNT_PASS=PUT_PASSWORD_HERE npm run e2e-staging`

   - NOTE: Replace `PUT_PASSWORD_HERE` with `corpora/backend/staging/auth0-secret.test_account_password` in AWS Secret Manager

1. local -> prod

   1. Manually check https://cellxgene.cziscience.com is working
   1. In `frontend/` directory, run `npm run e2e-prod`

   - NOTE: we don't run logged in tests in prod, since we don't want to accidentally
     add test data to prod database

1. local -> rdev

   1. Manually check your rdev link is working. E.g., https://pr-6062-frontend.rdev.single-cell.czi.technology
   1. In `frontend/` directory, run `RDEV_LINK=YOUR_RDEV_LINK_HERE TEST_ACCOUNT_PASS=PUT_PASSWORD_HERE npm run e2e-rdev`.

   - NOTE: Replace `YOUR_RDEV_LINK_HERE` with your rdev link. E.g., `https://pr-6062-frontend.rdev.single-cell.czi.technology`
   - NOTE: Replace `PUT_PASSWORD_HERE` with `corpora/backend/rdev/auth0-secret.test_account_password` in AWS Secret Manager

1. Running Auth-related Tests Locally

   For cases such as testing collection revision tests, the tests require seeded published collections to make the revision from. Tests can’t create and publish a new collection since in order to publish a newly created collection it needs at least one dataset, and our tests cannot upload a dataset to DropBox.

   1. Connect your local frontend app to the deployed `dev` API to fetch the data. This is done in the `configs.js` file.
   1. In `playwright.config.ts`, modify the `MANUAL_COOKIE` variable by supplying a valid cookie `value`.
      - This cookie is retrieved by logging into the `dev` environment and copying the cookie value.
      - **Make sure that you do NOT commit the cookie value into the repo.**
   1. Run your authenticated tests by adding the `USE_COOKIE=true` flag
      - Example: `HEADFUL=true SKIP_LOGIN=true USE_COOKIE=true npm run e2e -- -- tests/features/collection/revision.test.ts`
      - It is recommended to add `.only()` for the tests that you're interested in so that the whole test suite isn't run.

## How to debug Github action failed e2e tests

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
      1. In the test directory, run `npm run e2e-trace trace.zip`
         ![test results](https://user-images.githubusercontent.com/6309723/196766317-95afbb3c-2890-42af-a42e-0b1f7702d73c.png)
         ![Trace Viewer](https://user-images.githubusercontent.com/6309723/196768996-899f1086-13e8-4c3d-804a-eb050bfa6c71.gif)
1. If the information above is not enough to help you pinpoint the root cause of the test failure, we can use [Playwright Inspector](https://playwright.dev/docs/debug) to rerun the [test locally](#how-to-run-e2e-tests-locally) with `--debug` flag and debug the test steps live

   1. Find the test you want to debug.
      The `test-results` directory name is part of the test name, so you can use that to find the actual test case.

      For example, `test-results` directory name `features-collection-collection-Collection-Depl-5df78-te-a-collection-with-a-DOI-in-an-invalid-format-chromium` points us to `frontend/tests/features/collection/collection.test.ts` and test case `doesn't create a collection with a DOI in an invalid format`

   1. Once you find the test, add `.only` to the `test()` call, so Playwright will only run that test in the file
      ![Add .only](https://user-images.githubusercontent.com/6309723/196773185-8553273c-0cfd-4861-ac57-e2e21abfb76a.png)

   1. Use `await page.pause()` anywhere in the test code to pause Playwright execution

      ![page.pause()](https://user-images.githubusercontent.com/6309723/196776657-ac957fd6-935f-44a1-9188-490e85ad3e74.png)

   1. Since we're debugging Dev environment, we need to use the `local -> dev` command as shown in the [How](#how) section above.

   1. In `single-cell-data-portal/frontend/`, run `npm run e2e-dev -- --debug`.

      NOTE: To run in one specific browser, run `npm run e2e-dev -- --project chromium` (chromium|firefox|edge)

   1. Playwright will now spin up a browser window and a debugging console for you to inspect the test live!

      1. You will find the play button in the console, so you can play/pause the test as you see fit
      1. For more information, visit [Playwright Inspector](https://playwright.dev/docs/debug)

      ![Playwright Inspector](https://user-images.githubusercontent.com/6309723/196775284-81bd4853-b4c6-45c4-827b-22e97b97ee06.png)

## Troubleshooting

### preSetup Access Token

If you see an issue that looks like this
```
  ✘  1 [preSetup] › tests/common/playwright.global.preSetup.ts:11:7 › global preSetup › Get access token (0ms)
Error: ENOENT: no such file or directory, open '/var/folders/l5/ygnys3jj7n9f12p826j9448c0000gq/T/playwright-transform-cache-503/f7/context_f76b0cdd8a23a6ac47b3cccbc8fed32460f74ab0.js'

   at playwright.config.ts:11

   9 | import fs from "fs";
  10 | import { LOGIN_STATE_FILENAME } from "tests/common/constants";
> 11 | import { COMMON_PLAYWRIGHT_CONTEXT } from "tests/common/context";
     | ^
  12 | import { getFeatureFlags } from "tests/common/featureFlags";
  13 | import { SKIP_LOGIN } from "tests/common/constants";
  14 | import { shouldUseRdevToken } from "tests/utils/helpers";
```
You can resolve it by just removing the folder
```
rm -rf /var/folders/l5/ygnys3jj7n9f12p826j9448c0000gq/T/playwright-transform-cache-503
```

### AWS Region Missing in Setup

If you see an issue within the Setup that looks like this
```
  ✘  2 [setup] › tests/common/playwright.global.setup.ts:10:7 › global setup › login (13.3s)
🔐🪵 Logging in...
Error: Region is missing
```

You can resolve this by updating your AWS config file to include the following snippet. On a Linux/Mac, this should be in `~/.aws/config`.
```
[default]
role_arn       = arn:aws:iam::***:role/poweruser // Ask someone else on the team for what this should be
source_profile = single-cell-dev-poweruser
region         = us-west-2
```

---

# Test/Env Differential

Where tests are skipped vs. run in different environments.

<table >
<colgroup>
<col />
<col />
<col />
<col />
<col />
<col />
<col />
<col />
</colgroup>
<thead>
  <tr>
    <th></th>
    <th colspan="7">TEST_ENV (defaults to local)</th>
  </tr>
  <tr>
    <th></th>
    <th>local</th>
    <th>localProd</th>
    <th>staging</th>
    <th>prod</th>
    <th>rdev</th>
    <th>happy</th>
    <th>dev</th>
  </tr>
</thead>
<tbody>
  <tr>
    <td>
      <ul>
         <li>"Collection Revision"</li>
      </ul>
   </td>
    <td>skip</td>
    <td>skip</td>
    <td>describe</td>
    <td>skip</td>
    <td>skip</td>
    <td>skip</td>
    <td>describe</td>
  </tr>
  <tr>
    <td>
      <ul>
         <li>`./features/wheresMyGene.test.ts`</li>
         <li>"creates and deletes a collection"</li>
         <li>"dataset order"</li>
         <li>"publish a collection"</li>
         <li>"invalid DOIs"</li>
      </ul>

   </td>
    <td>skip<br></td>
    <td>skip<br></td>
    <td>describe</td>
    <td>describe<br></td>
    <td>skip</td>
    <td>skip</td>
    <td>describe</td>
  </tr>
</tbody>
</table>

# Best Practices

## Define test coverage

1.  Happy Paths:

    Business critical features must at least have the most commonly exercised happy path scenarios covered, so we ensure the majority of the users won't be affected by an uncommon bug

1.  Incremental coverage:

    As uncommon bugs are discovered, create tickets to track and write tests to ensure those edge cases won't happen again. This way we gradually build up the test coverage of our application for quality improvement

## Build meaningful, efficient, and focused tests

1.  Break workflows down with smaller tests:

    Smaller and focused tests increase the team's ability to pinpoint root cause faster, tests also run faster, and thus decrease the bug fix time. For examples:

    1. Testing for a collection revision, seed published collections for the test account, so the revision tests don't waste time creating collections before testing revisions and risk running into bugs related to collection creation, which should be covered by collection creation tests, not revision tests.

    1. Loggedin in tests should avoid exercising the login flow and reuse the authenticated browser state from the login test itself

1.  Create meaningful assertions to verify expected behaviors:

    1. Before writing any tests, think about what major evidences and behaviors that need to happen to prove that the workflow is working. E.g., certain HTML elements exist, certain content exists, URL contains certain query params, etc..

    1. Step into your users' shoes:

       Tests should reflect workflows that our users actually perform. In other words, your test steps should mimic how a user interacts with the application, so we exercise code paths and features our users will use

1.  Anticipate and avoid flakiness:

    End to end tests could go wrong due to many factors outside of your control, such as temporary network blips, spikes in traffic, JS listeners not set up in time to observe user interactions, etc.. Our job as test authors is to anticipate such common conditions and write robust tests that take those conditions into account, so we avoid false alerts as much as possible

    1. Mindset:

       Whenever you write tests, always ask yourself **what could go wrong in your tests**, especially the flaky conditions that are out of your control. Remember Murphy's Law, if anything can go wrong, it will. So identify and handle those factors

    1. Check your test assumptions:

       In your test setup, are you picking the first item from a dropdown? If so, your test assumes that the first item will always be the same one, which is likely untrue as data and/or default sort order changes. As a result, your subsequent assertions must NOT rely on hard coded expected state that tie to the first item today, instead, parse the expected state at run time, so your test is resilient against state changes

       For example, in a Gene Expression test, we select the first tissue from the dropdown list, and instead of hard coding the expected number of cell types from the tissue, we parse the UI to find the expected number of cell types at run time. This way, when BE added/removed cell types associated to the tissue, our test would still pass

    1. Avoid using hard coded wait:

       This is a very common code smell that should be avoided as much as possible.

       For example, it's very easy to write a test to wait for 3 seconds after navigating to a page, hoping 3 seconds would be enough for the page to be loaded and ready for testing. However, due to the flaky conditions listed above, 3 seconds could either be too long or too short, depending on the machine and the network traffic at the time. So it's more reliable to wait for certain conditions and/or events to happen over a hard coded wait. Such as waiting for page load event to happen, or once certain content has been rendered on the page before the actual test begins

    1. Use utility function `tryUntil()`: There will be times when your test assertions will happen before the application is in the expected state due to the flaky factors above. So we have a utility function `tryUntil()` that allows you to retry defined actions/assertions until your expected condition is met. For examples, you can retry clicking on a button until the modal element exists, or retry asserting a certain element exists before throwing an error. There are examples available in the tests that you can look for them by global searching for `tryUntil`

    1. Use `data-testid` HTML attribute for target elements: Since our application's HTML structure changes over time, it's unreliable to select an element based on properties and structures that could easily change over time, such as css classes and element structures (e.g., first div child of a parent). The more reliable way is to add a `data-testid` attribute to your test target, so when the element and/or its context changes, we can still reliably target the element

## Maintain application in a predictable state

1.  Create idempotent tests that restore the application state to the state before the tests are run. This way we eliminate side effects and byproducts of said tests that could cause unexpected errors in tests.

    For example, My Collections page tests used to generate a large amount of private collections that were never cleaned up, which led to My Collections page tests to start timing out due to the page loading very slowly. The solution was to delete the transient collections after each test, so we only have an expected amount of data when each test runs

1.  Write atomic tests:

    Atomic tests mean that each test case should have zero dependencies on the test cases before it, so we eliminate the possibility of side effects from other tests polluting the test environment of your current test.

    This helps avoid flakiness that are hard to reproduce and debug, and also allows us to auto retry a test case that fails temporarily due to flaky conditions

## Keep code DRY and readable

1.  Avoid nesting tests with `beforeEach()` and `afterEach()` when you can extract reusable functions to call explicitly inside each test case. This helps improve readability and maintenance, since all the setup and teardown are isolated within a test case. [(Read more)](https://kentcdodds.com/blog/avoid-nesting-when-youre-testing)

1.  Avoid unnecessary steps that add complexity to tests

    When your tests exercise code paths and features that are not related to your test case, you risk adding noises in test reports, slow down test runs, slow down your debugging speed, and increase points of failure and cost of maintenance

    For examples, performing test steps that users don't do or writing assertions that are irrelevant to your test case.

1.  Write test descriptions that are readable to your team and your future self

    This way it's easier to understand the scope and intention of the test in the test report and helps your team decide whether to add new test cases or modify the existing one when new features are added and/or refactors happen.

    As a bonus point, a well written test description helps code reviewers provide valuable feedback, such as catching test steps that don't assert enough for the intended test case, or suggesting ways to help improve test accuracy and efficiency

## Resources

1. [Datadog testing guide](https://www.datadoghq.com/resources/frontend-monitoring-best-practices/)

   1. The first two sections should be enough! **Best practices for creating
      end-to-end tests** and **Best practices for
      maintaining end-to-end tests**

1. [Avoid nesting when you're testing](https://kentcdodds.com/blog/avoid-nesting-when-youre-testing)
