import { goToPage } from "tests/utils/helpers";
import { getTestID } from "../utils/selectors";

describe("Homepage", () => {
  it("renders the expected elements", async () => {
    await goToPage();

    await expect(page).toHaveSelector(getTestID("collections-header"));
    await expect(page).toHaveSelector(getTestID("logo"));
    await expect(page).toHaveSelector(getTestID("collection-link"));
    await expect(page).not.toHaveSelector(getTestID("visibility-tag"));
    await page.click(getTestID("collection-link"));
    await expect(page.url()).toContain("collections");
  });
});
