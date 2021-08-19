import { ElementHandle } from "playwright";
import { TEST_ENV } from "tests/common/constants";
import { TEST_PASSWORD, TEST_URL, TEST_USERNAME } from "../common/constants";
import { getText } from "./selectors";

export const TIMEOUT_MS = 3 * 1000;

export const describeIfDeployed =
  TEST_ENV.includes("local") || TEST_ENV === "prod" ? describe.skip : describe;

export async function goToPage(url: string = TEST_URL): Promise<void> {
  await page.goto(url);
}

export async function login(): Promise<void> {
  goToPage();

  expect(process.env.TEST_ACCOUNT_PASS).toBeDefined();

  try {
    await expect(page).toHaveSelector(getText("My Collections"), {
      timeout: TIMEOUT_MS,
    });
  } catch (error) {
    const username = TEST_USERNAME;

    await page.click(getText("Log In"));

    await page.fill('[name="Username"], [name="email"]', username);
    await page.fill('[name="Password"], [name="password"]', TEST_PASSWORD);

    await page.click('[value="login"], [name="submit"]');

    await tryUntil(() => {
      expect(page.url()).toContain(TEST_URL);
    });
  }
}

export async function tryUntil(
  assert: () => void,
  maxRetry = 50
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
      savedError = error;
      await page.waitForTimeout(WAIT_FOR_MS);
    }
  }

  if (retry === maxRetry) {
    savedError.message += " tryUntil() failed";
    throw savedError;
  }
}

export async function getInnerText(selector: string): Promise<string> {
  await tryUntil(() => page.waitForSelector(selector));

  const element = (await page.$(selector)) as ElementHandle<
    SVGElement | HTMLElement
  >;

  return element.innerText();
}
