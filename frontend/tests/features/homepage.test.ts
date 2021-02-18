import { goToPage } from "tests/utils/helpers";
import { getTestClass, getTestID } from "../utils/selectors";

describe("Homepage", () => {
  it("renders the expected elements", async () => {
    await goToPage();

    await expect(page).toHaveSelector(getTestID("collections-header"));
    await expect(page).toHaveSelector(getTestID("logo"));
    await expect(page).toHaveSelector(getTestClass("collection-row"));
    await expect(page).not.toHaveSelector(getTestClass("visibility-tag"));
    await page.click(getTestClass("collection-row"));
    await expect(page.url()).toContain("collections");
  });
});
