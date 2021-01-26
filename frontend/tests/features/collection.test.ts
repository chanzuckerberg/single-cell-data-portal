import { TEST_ENV, TEST_URL } from "tests/common/constants";
import { goToPage } from "tests/utils/helpers";
import { getText } from "tests/utils/selectors";

const describeIfCalledByDevEnv = TEST_ENV === "dev" ? describe : describe.skip;

describeIfCalledByDevEnv("CRUD Collection", async () => {
  it("creates and deletes a collection", async () => {
    await goToPage(TEST_URL);
    // LOG IN
    await page.click(getText("Log In"));
    await page.waitForNavigation({ url: TEST_URL, waitUntil: "networkidle" });
    // Create collection
    await page.click(getText("Create Collection"));
    await page.fill("#name", "TEST COLLECTION");
    await page.fill("#description", "TEST DESCRIPTION");
    await page.fill("#contact-name", "TEST");
    await page.fill("#contact-email", "TEST@example.coms");
    await page.click("input[type=checkbox]");
    const [request, response] = await Promise.all([
      page.waitForEvent("request"),
      page.waitForEvent("response"),
      await page.click(getText("Create")),
    ]);

    const cookie = request.headers().cookie;

    const collectionID = (await response.json()) as object;

    await expect(page).toHaveSelector("h3");
    // CLICK DELETE COLLECTION BUTTON HERE
    // TEMP: FETCH ENDPOINT TO DELETE COLLECTION
    const deleteResponse = await fetch(
      `https://api.cellxgene.dev.single-cell.czi.technology/dp/v1/collections/${collectionID}`,
      { headers: { cookie }, method: "DELETE" }
    );

    expect(deleteResponse.status).toBe(202);
    await page.reload();
    await expect(page).not.toHaveSelector("h3");
  });
});
