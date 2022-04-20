import { goToPage } from "tests/utils/helpers";
import { getTestID } from "../utils/selectors";

const TIMEOUT_MS = 3 * 1000;

describe("Homepage", () => {
  it("renders the expected elements", async () => {
    await goToPage();
    await expect(page).toHaveSelector(getTestID("logo"));
    await expect(page).toHaveSelector(getTestID("collection-link"));
    await expect(page).not.toHaveSelector(getTestID("visibility-tag"), {
      timeout: TIMEOUT_MS,
    });

    await Promise.all([
      page.waitForNavigation(),
      page.click(getTestID("collection-link")),
    ]);

    expect(page.url()).toContain("collections");
  });
});
