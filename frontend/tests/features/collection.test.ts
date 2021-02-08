import { ROUTES } from "src/common/constants/routes";
import { apiTemplateToUrl } from "src/common/utils/apiTemplateToUrl";
import { API_URL } from "src/configs/configs";
import { TEST_ENV } from "tests/common/constants";
import { goToPage, login } from "tests/utils/helpers";
import { getTestTag, getText } from "tests/utils/selectors";

const describeIfDeployed = TEST_ENV.includes("local")
  ? describe.skip
  : describe;

const TEST_COLLECTION = {
  contactEmail: "TEST@example.com",
  contactName: "TEST NAME",
  description: "TEST DESCRIPTION",
  name: "TEST COLLECTION",
};

describe("Collection", async () => {
  describeIfDeployed("Logged In Tests", () => {
    it("creates and deletes a collection", async () => {
      await login();

      const collectionId = await createCollection();

      // Try delete
      await page.click("Delete");
      await page.click("Delete Collection");

      await goToPage(
        apiTemplateToUrl(API_URL + ROUTES.PRIVATE_COLLECTION, {
          id: collectionId,
        })
      );
      await expect(page).not.toHaveText(TEST_COLLECTION.name);
    });

    describe("Publish a collection", () => {
      describe("when no dataset", () => {
        it("shows disabled publish button", async () => {
          await login();

          await createCollection();

          const publishButton = await page.$(
            getTestTag("publish-collection-button")
          );

          expect(publishButton?.getAttribute("disabled")).toBe("");
        });
      });
    });
  });
});

async function createCollection(): Promise<string> {
  await page.click(getText("Create Collection"));
  await page.fill("#name", TEST_COLLECTION.name);
  await page.fill("#description", TEST_COLLECTION.description);
  await page.fill("#contact-name", TEST_COLLECTION.contactName);
  await page.fill("#contact-email", TEST_COLLECTION.contactEmail);
  await page.click(getText("I agree to cellxgene's data submission policies."));

  const [, response] = await Promise.all([
    page.waitForEvent("request"),
    page.waitForEvent("response"),
    await page.click(getText("Create")),
  ]);

  const { collectionId } = (await response.json()) as {
    collectionId: string;
  };

  await expect(page).toHaveText(TEST_COLLECTION.name);

  return collectionId;
}
