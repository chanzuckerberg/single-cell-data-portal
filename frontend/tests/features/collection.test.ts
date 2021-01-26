import { ROUTES } from "src/common/constants/routes";
import { apiTemplateToUrl } from "src/common/utils/apiTemplateToUrl";
import { API_URL } from "src/configs/configs";
import { TEST_ENV, TEST_URL } from "tests/common/constants";
import { goToPage } from "tests/utils/helpers";
import { getText } from "tests/utils/selectors";

const describeIfCalledByDevEnv = TEST_ENV === "dev" ? describe : describe.skip;

/*  eslint-disable sort-keys */
const TEST_COLLECTION = {
  name: "TEST COLLECTION",
  description: "TEST DESCRIPTION",
  contactName: "TEST NAME",
  contactEmail: "TEST@example.com",
};
/* eslint-enable sort-keys */

describeIfCalledByDevEnv("Collection", async () => {
  it("creates and deletes a collection", async () => {
    await goToPage(TEST_URL);
    // LOG IN
    await page.click(getText("Log In"));
    await page.waitForNavigation({ url: TEST_URL, waitUntil: "networkidle" });
    // Create collection
    await page.click(getText("Create Collection"));
    await page.fill("#name", TEST_COLLECTION.name);
    await page.fill("#description", TEST_COLLECTION.description);
    await page.fill("#contact-name", TEST_COLLECTION.contactName);
    await page.fill("#contact-email", TEST_COLLECTION.contactEmail);
    await page.click(
      getText("I agree to cellxgene's data submission policies.")
    );
    const [request, response] = await Promise.all([
      page.waitForEvent("request"),
      page.waitForEvent("response"),
      await page.click(getText("Create")),
    ]);

    const cookie = request.headers().cookie;

    const { collectionID } = (await response.json()) as {
      collectionID: string;
    };

    await expect(page).toHaveText(TEST_COLLECTION.name);

    // Try delete
    await page.click("Delete");
    await page.click("Delete Collection");

    await goToPage(
      apiTemplateToUrl(API_URL + ROUTES.PRIVATE_COLLECTION, {
        id: collectionID,
      })
    );
    await expect(page).not.toHaveText(TEST_COLLECTION.name);
  });
});
