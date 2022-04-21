import { goToPage } from "tests/utils/helpers";
import { getTestID } from "../utils/selectors";

describe("Homepage", () => {
  it("renders the expected elements", async () => {
    await goToPage();
    expect(page.url()).toContain("datasets");
    await expect(page).toHaveSelector(getTestID("logo"));
    await expect(page).toHaveSelector(getTestID("collection-link"));

    await Promise.all([
      page.waitForNavigation(),
      page.click(getTestID("collection-link")),
    ]);

    expect(page.url()).toContain("collections");
  });
});
