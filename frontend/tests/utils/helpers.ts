export async function goToPage(url: string) {
  await page.goto(url);
}

export async function tryUntil(assert: () => void) {
  const MAX_RETRY = 10;
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
