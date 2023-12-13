import { expect } from "@playwright/test";
import { ROUTES } from "src/common/constants/routes";
import { goToPage, tryUntil } from "tests/utils/helpers";
import { TEST_URL } from "../../common/constants";
import { test } from "tests/common/test";

const { describe } = test;

describe("api/deployed_version", () => {
  test("Returns commit SHA", async ({ page }) => {
    await goToPage(`${TEST_URL}${ROUTES.DEPLOYED_VERSION}`, page);

    await tryUntil(
      async () => {
        const commitSha = await page.textContent("body");
        const regex = /(?<="Data Portal":")\w+/;
        const match = regex.exec(commitSha || "");
        expect(match).toBeTruthy();
      },
      { page }
    );
  });
});
