import { TEST_URL } from "../common/constants";

export const TIMEOUT_MS = 3 * 1000;

export async function goToPage(url: string = TEST_URL) {
  await page.goto(url);
}

export async function login() {
  return;
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
