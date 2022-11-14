import { ElementHandle, expect, Page } from "@playwright/test";
import { TEST_ENV } from "tests/common/constants";
import { TEST_PASSWORD, TEST_URL, TEST_USERNAME } from "../common/constants";
import { getText } from "./selectors";

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

  expect(process.env.TEST_ACCOUNT_PASS).toBeDefined();

  const cookies = await (await page.context()).cookies();

  if (cookies.length) return;

  const username = TEST_USERNAME;

  await page.click(getText("Log In"));

  await page.fill('[name="Username"], [name="email"]', username);
  await page.fill('[name="Password"], [name="password"]', TEST_PASSWORD);

  await page.click('[value="login"], [name="submit"]');

  await tryUntil(
    () => {
      expect(page.url()).toContain(TEST_URL);
    },
    { page }
  );
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
