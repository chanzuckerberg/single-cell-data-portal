import {
  GetSecretValueCommand,
  SecretsManagerClient,
} from "@aws-sdk/client-secrets-manager";
import { ElementHandle, expect, Locator, Page } from "@playwright/test";
import { TEST_ENV } from "tests/common/constants";
import { LOGIN_STATE_FILENAME, TEST_URL } from "../common/constants";
import {
  CELL_TYPE_LABELS_ID,
  ERROR_NO_TESTID_OR_LOCATOR,
  GENE_LABELS_ID,
} from "../common/constants";
import { TISSUE_NAME_LABEL_CLASS_NAME } from "src/views/WheresMyGeneV2/components/HeatMap/components/YAxisChart/constants";

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
  endpoint = String(new URL(process.env["BOTO_ENDPOINT_URL"] as string));
} catch (e) {
  console.log("BOTO_ENDPOINT_URL not assigned, assuming running on deployment");
}

const client = new SecretsManagerClient({
  // (thuang): Do NOT pass empty string to `endpoint` or it will throw `TypeError: Invalid URL`
  endpoint,
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

const GO_TO_PAGE_TIMEOUT_MS = 2 * 60 * 1000;

export async function goToPage(
  url: string = TEST_URL,
  page: Page
): Promise<void> {
  await tryUntil(
    async () => {
      await Promise.all([
        page.waitForLoadState("networkidle"),
        page.goto(url, { timeout: GO_TO_PAGE_TIMEOUT_MS }),
      ]);
    },
    { page }
  );
}

export async function login(page: Page): Promise<void> {
  console.log("Logging in...");

  await goToPage(undefined, page);

  const { username, password } = await getTestUsernameAndPassword();

  expect(typeof username).toBe("string");
  expect(typeof password).toBe("string");

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

  console.log(`Login success!`);
}

export async function scrollToPageBottom(page: Page): Promise<void> {
  return page.evaluate(() =>
    window.scrollTo(0, document.documentElement.scrollHeight)
  );
}
interface TryUntilConfigs {
  maxRetry?: number;
  page: Page;
  silent?: boolean;
}

export async function tryUntil(
  assert: () => void,
  { maxRetry = 50, page, silent = false }: TryUntilConfigs
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

      if (!silent) {
        console.log("⚠️  tryUntil error-----------------START");
        console.log(savedError.message);
        console.log("⚠️  tryUntil error-----------------END");
      }

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

export async function selectNthOption(page: Page, number: number) {
  // (thuang): Since the first option is now active, we need to offset by 1
  const step = number - 1;

  for (let i = 0; i < step; i++) {
    await page.keyboard.press("ArrowDown");
  }

  await page.keyboard.press("Enter");
  await page.keyboard.press("Escape");
}

export async function waitForElement(page: Page, testId: string) {
  await tryUntil(
    async () => {
      await expect(page.getByTestId(testId)).not.toHaveCount(0);
    },
    { page }
  );
}

export async function getButtonAndClick(page: Page, testID: string) {
  await tryUntil(
    async () => {
      await page.getByTestId(testID).click();
    },
    { page }
  );
}

// for when there are multiple buttons with the same testID
export async function getFirstButtonAndClick(page: Page, testID: string) {
  await tryUntil(
    async () => {
      const buttons = await page.getByTestId(testID).elementHandles();
      await buttons[0].click();
    },
    { page }
  );
}

export async function getCellTypeFmgButtonAndClick(
  page: Page,
  cellType: string
) {
  await waitForElement(page, CELL_TYPE_LABELS_ID);

  await tryUntil(
    async () => {
      await page
        .getByRole("img", {
          name: "display marker genes for " + cellType,
        })
        .click();
    },
    { page }
  );
}

export async function clickDropdownOptionByName({
  page,
  testId,
  name,
}: {
  page: Page;
  testId: string;
  name: string;
}) {
  await page.getByTestId(testId).click();
  await page.getByRole("option").filter({ hasText: name }).click();
  await page.keyboard.press("Escape");
}

export async function getGeneNames(page: Page) {
  return getNames({ page, locator: page.locator(GENE_LABELS_ID) });
}

export async function getCellTypeNames(page: Page) {
  return getNames({ page, testId: CELL_TYPE_LABELS_ID });
}

// (alec) use this instead of locator.count() to make sure that the element is actually present
// when counting
export async function countLocator(locator: Locator) {
  return (await locator.elementHandles()).length;
}

export async function getNames({
  page,
  testId,
  locator,
}: {
  page: Page;
  testId?: string;
  locator?: Locator;
}): Promise<string[]> {
  let labelsLocator: Locator;
  if (testId) {
    labelsLocator = page.getByTestId(testId);
  } else if (locator) {
    labelsLocator = locator;
  } else {
    throw Error(ERROR_NO_TESTID_OR_LOCATOR);
  }

  if ((await labelsLocator.count()) === 0) return [];
  await tryUntil(
    async () => {
      const names = await labelsLocator.allTextContents();
      expect(typeof names[0]).toBe("string");
    },
    { page }
  );

  return labelsLocator.allTextContents();
}

export async function clickUntilOptionsShowUp({
  page,
  testId,
  locator,
}: {
  page: Page;
  testId?: string;
  locator?: Locator;
}) {
  // either testId or locator must be defined, not both
  // locator is used when the element cannot be found using just the test Id from the page
  await tryUntil(
    async () => {
      if (testId) {
        await page.getByTestId(testId).click();
      } else if (locator) {
        await locator.click();
      } else {
        throw Error(ERROR_NO_TESTID_OR_LOCATOR);
      }
      await page.getByRole("tooltip").getByRole("option").elementHandles();
    },
    { page }
  );
}

export async function clickUntilDownloadModalShowsUp({
  page,
  testId,
  locator,
}: {
  page: Page;
  testId?: string;
  locator?: Locator;
}) {
  await tryUntil(
    async () => {
      if (testId) {
        await page.getByTestId(testId).click();
      } else if (locator) {
        await locator.click();
      } else {
        throw Error(ERROR_NO_TESTID_OR_LOCATOR);
      }
      await page.locator(".bp4-dialog").elementHandle();
    },
    { page }
  );
}

export async function clickUntilSidebarShowsUp({
  page,
  testId,
  locator,
}: {
  page: Page;
  testId?: string;
  locator?: Locator;
}) {
  await tryUntil(
    async () => {
      if (testId) {
        await page.getByTestId(testId).click();
      } else if (locator) {
        await locator.click();
      } else {
        throw Error(ERROR_NO_TESTID_OR_LOCATOR);
      }
      await page.locator(".bp4-drawer-header").elementHandle();
    },
    { page }
  );
}

// (thuang): This only works when a dropdown is open
export async function selectFirstOption(page: Page) {
  await selectFirstNOptions(1, page);
}

export async function selectFirstNOptions(count: number, page: Page) {
  for (let i = 0; i < count; i++) {
    await page.keyboard.press("ArrowDown");
    await page.keyboard.press("Enter");
  }

  await page.keyboard.press("Escape");
}

export async function takeSnapshotOfMetaTags(name: string, page: Page) {
  const allMetaTags = await page.locator("meta").all();

  await tryUntil(
    async () => {
      const allMetaTagsHTML = (
        await Promise.all(
          allMetaTags.map(async (metaTag) =>
            metaTag.evaluate((node) => node.outerHTML)
          )
        )
      ).sort();

      expect(JSON.stringify(allMetaTagsHTML)).toMatchSnapshot({
        name: name + "-seoMetaTags.txt",
      });
    },
    { page }
  );
}

export async function expandTissue(page: Page, tissueName: string) {
  await tryUntil(
    async () => {
      const beforeCellTypeNames = await getCellTypeNames(page);
      await page
        .getByTestId(`cell-type-labels-${tissueName}`)
        .getByTestId(TISSUE_NAME_LABEL_CLASS_NAME)
        .click();
      const afterCellTypeNames = await getCellTypeNames(page);
      expect(afterCellTypeNames.length).toBeGreaterThan(
        beforeCellTypeNames.length
      );
    },
    { page }
  );
}

export async function waitForLoadingSpinnerToResolve(page: Page) {
  await page.getByText("Loading").first().waitFor({ state: "hidden" });
}
