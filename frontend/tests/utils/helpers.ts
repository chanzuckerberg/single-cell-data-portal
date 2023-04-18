import {
  GetSecretValueCommand,
  SecretsManagerClient,
} from "@aws-sdk/client-secrets-manager";
import { ElementHandle, expect, Page } from "@playwright/test";
import { TEST_ENV } from "tests/common/constants";
import { LOGIN_STATE_FILENAME, TEST_URL } from "../common/constants";

/**
 * (thuang): From oauth/users.json
 * `frontend/tests/common/constants.ts` `API_URL` needs to be `"https://backend.corporanet.local:5000"`
 * If `API_URL` is dev, staging, prod, or rdev, use their corresponding OIDC credentials
 */
const LOCAL_CONTAINER_TEST_USERNAME = "User1";
const LOCAL_CONTAINER_TEST_PASSWORD = "pwd";

const DEPLOYED_TEST_USERNAME = "user@example.com";

let endpoint;

try {
  endpoint = new URL(process.env["BOTO_ENDPOINT_URL"] as string);
} catch (e) {
  console.log("BOTO_ENDPOINT_URL not assigned, assuming running on deployment");
  endpoint = "";
}

const client = new SecretsManagerClient({
  endpoint: String(endpoint),
});

const deployment_stage = process.env.DEPLOYMENT_STAGE || "test";

const secretValueRequest = {
  SecretId: `corpora/backend/${deployment_stage}/auth0-secret`,
};

const command = new GetSecretValueCommand(secretValueRequest);

export const TIMEOUT_MS = 3 * 1000;

// (seve): We use TEST_ENV to describe the environment that playwright is running against. Sometimes the FE tests are run against a local instance of the app which points at a deployed instance of the backend.

//(thuang): BE API doesn't work in local happy
const TEST_ENVS_DEV_STAGING_PROD = ["dev", "staging", "prod"];

export const isDevStagingProd = TEST_ENVS_DEV_STAGING_PROD.includes(TEST_ENV);

// Skip tests unless environment is dev or staging; used by tests that require a deployed environment but also modify
// environment data (e.g. creating collections, which should be avoided in prod).
const TEST_ENVS_DEV_STAGING = ["dev", "staging"];

export const isDevStaging = TEST_ENVS_DEV_STAGING.includes(TEST_ENV);

export async function goToPage(
  url: string = TEST_URL,
  page: Page
): Promise<void> {
  await page.goto(url);
}

export async function login(page: Page): Promise<void> {
  await goToPage(undefined, page);

  const { username, password } = await getTestUsernameAndPassword();

  expect(username).toBeDefined();

  await page.getByText("Log In").click();

  await page.fill('[name="Username"], [name="email"]', username);
  await page.fill('[name="Password"], [name="password"]', password);

  await page.click('[value="login"], [name="submit"]');

  /**
   * (thuang): This is needed to ensure the browser has been navigated back to
   * our app
   */
  await tryUntil(
    () => {
      expect(page.url()).toContain(TEST_URL);
    },
    { page }
  );

  console.log("setting storage state...");

  await page.context().storageState({ path: LOGIN_STATE_FILENAME });
}

interface TryUntilConfigs {
  maxRetry?: number;
  page: Page;
}

export async function tryUntil(
  assert: () => void,
  { maxRetry = 50, page }: TryUntilConfigs
): Promise<void> {
  const WAIT_FOR_MS = 200;

  let retry = 0;

  let savedError: Error = new Error();

  while (retry < maxRetry) {
    try {
      await assert();

      break;
    } catch (error) {
      retry += 1;
      savedError = error as Error;
      await page.waitForTimeout(WAIT_FOR_MS);
    }
  }

  if (retry === maxRetry) {
    savedError.message += " tryUntil() failed";
    throw savedError;
  }
}

export async function getInnerText(
  selector: string,
  page: Page
): Promise<string> {
  await tryUntil(() => page.waitForSelector(selector), { page });

  const element = (await page.$(selector)) as ElementHandle<
    SVGElement | HTMLElement
  >;

  return element.innerText();
}

async function getTestUsernameAndPassword() {
  try {
    const secret = JSON.parse(
      (await client.send(command)).SecretString || "null"
    );

    const { test_account_username, test_account_password } = secret;

    const deployedTestUsername =
      test_account_username || DEPLOYED_TEST_USERNAME;

    const userNames = {
      dev: deployedTestUsername,
      happy: LOCAL_CONTAINER_TEST_USERNAME,
      local: LOCAL_CONTAINER_TEST_USERNAME,
      localProd: LOCAL_CONTAINER_TEST_USERNAME,
      prod: deployedTestUsername,
      rdev: deployedTestUsername,
      staging: deployedTestUsername,
    };

    return {
      password: test_account_password || LOCAL_CONTAINER_TEST_PASSWORD || "",
      username: userNames[TEST_ENV] || deployedTestUsername,
    };
  } catch (error) {
    console.error(error);
    return {
      password: null,
      username: null,
    };
  }
}

export async function waitForElementToBeRemoved(page: Page, selector: string) {
  await tryUntil(
    async () => {
      const element = await page.$(selector);
      await expect(element).toBeNull();
    },
    { page }
  );
}
