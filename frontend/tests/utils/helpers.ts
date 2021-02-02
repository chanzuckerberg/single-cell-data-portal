import { TEST_URL } from "tests/common/constants";
import { getText } from "./selectors";

export async function goToPage(url: string = TEST_URL) {
  await page.goto(url);
}

export async function login(url: string = TEST_URL) {
  await goToPage(url);
  await page.click(getText("Log In"));
  await page.waitForNavigation({ url, waitUntil: "networkidle" });
}

export async function tryUntil(assert: () => void) {
  const MAX_RETRY = 50;
  const WAIT_FOR_MS = 200;

  let retry = 0;

  while (retry < MAX_RETRY) {
    try {
      await assert();

      break;
    } catch (error) {
      retry += 1;

      await page.waitForTimeout(WAIT_FOR_MS);
    }
  }

  if (retry === MAX_RETRY) {
    throw Error("tryUntil() assertion failed!");
  }
}
