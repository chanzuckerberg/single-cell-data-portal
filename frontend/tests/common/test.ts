import { test as base, expect } from "@playwright/test";
import { TOO_MANY_REQUESTS_ERROR_MESSAGE_PREFIX } from "src/common/networkGuard";

/**
 * (thuang): Use this `test` object instead of the one from playwright/test
 */
export const test = base.extend<{ saveLogs: void }>({
  page: async ({ page }, use) => {
    const errors: string[] = [];

    page.on("console", (message) => {
      if (message.type() !== "error") return;

      const text = message.text();

      if (text.includes(TOO_MANY_REQUESTS_ERROR_MESSAGE_PREFIX)) {
        errors.push(text);
        console.error(text);
        throw Error(text);
      }
    });

    await use(page);

    expect(errors).toHaveLength(0);
  },
});
