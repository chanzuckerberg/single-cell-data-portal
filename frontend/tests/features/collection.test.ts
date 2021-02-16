import { ROUTES } from "src/common/constants/routes";
import { apiTemplateToUrl } from "src/common/utils/apiTemplateToUrl";
import { API_URL } from "src/configs/configs";
import { BLUEPRINT_SAFE_TYPE_OPTIONS, TEST_ENV } from "tests/common/constants";
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
      await goToPage();
      await login();

      const collectionId = await createCollection();

      // Try delete
      await page.click(getText("Delete"));
      await page.click(getText("Delete Collection"));

      await goToPage(
        apiTemplateToUrl(API_URL + ROUTES.PRIVATE_COLLECTION, {
          id: collectionId,
        })
      );
      await expect(page).not.toHaveSelector(getText(TEST_COLLECTION.name));
    });

    describe("Publish a collection", () => {
      describe("when no dataset", () => {
        it("shows disabled publish button", async () => {
          await goToPage();
          await login();

          await createCollection();

          const publishButton = await page.$(
            getTestTag("publish-collection-button")
          );

          expect(await publishButton?.getAttribute("disabled")).toBe("");
        });
      });
    });
  });
});

async function createCollection(): Promise<string> {
  await page.click(getText("Create Collection"));
  await page.type("#name", TEST_COLLECTION.name, BLUEPRINT_SAFE_TYPE_OPTIONS);
  await page.type(
    "#description",
    TEST_COLLECTION.description,
    BLUEPRINT_SAFE_TYPE_OPTIONS
  );
  await page.type(
    "#contact-name",
    TEST_COLLECTION.contactName,
    BLUEPRINT_SAFE_TYPE_OPTIONS
  );
  await page.type(
    "#contact-email",
    TEST_COLLECTION.contactEmail,
    BLUEPRINT_SAFE_TYPE_OPTIONS
  );
  await page.click(getText("I agree to cellxgene's data submission policies."));
  const [response] = await Promise.all([
    page.waitForEvent("response"),
    page.click(getTestTag("create button")),
  ]);

  const { collectionId } = (await response.json()) as {
    collectionId: string;
  };

  await expect(page).toHaveSelector(getText(TEST_COLLECTION.name));

  return collectionId;
}
