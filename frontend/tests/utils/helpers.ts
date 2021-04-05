import { TEST_PASSWORD, TEST_URL, TEST_USERNAME } from "../common/constants";
import { getText } from "./selectors";

export const TIMEOUT_MS = 3 * 1000;

export async function goToPage(url: string = TEST_URL) {
  await page.goto(url);
}

export async function login() {
  goToPage();
  expect(process.env.TEST_ACCOUNT_PASS).toBeDefined();
  try {
    await expect(page).toHaveSelector(getText("My Collections"), {
      timeout: TIMEOUT_MS,
    });
  } catch (error) {
    const url = page.url();
    const username = TEST_USERNAME;
    await page.click(getText("Log In"));
    await page.fill('[name="Username"], [name="email"]', username);
    await page.fill('[name="Password"], [name="password"]', TEST_PASSWORD);
    await Promise.all([
      page.waitForNavigation({ waitUntil: "networkidle" }),
      page.click('[value="login"], [name="submit"]'),
    ]);

    expect(page.url()).toContain(url);
  }
}

export async function tryUntil(assert: () => void) {
  const MAX_RETRY = 50;
  const WAIT_FOR_MS = 200;

  let retry = 0;

  let savedError: Error = new Error();

  while (retry < MAX_RETRY) {
    try {
      await assert();

      break;
    } catch (error) {
      retry += 1;
      savedError = error;
      await page.waitForTimeout(WAIT_FOR_MS);
    }
  }

  if (retry === MAX_RETRY) {
    savedError.message += " tryUntil() failed";
    throw savedError;
  }
}
