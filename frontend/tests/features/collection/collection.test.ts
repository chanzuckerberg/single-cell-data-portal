import { ROUTES } from "src/common/constants/routes";
import { Collection } from "src/common/entities";
import { sortByCellCountDescending } from "src/components/Collection/components/CollectionDatasetsGrid/components/DatasetsGrid";
import { BLUEPRINT_SAFE_TYPE_OPTIONS, TEST_URL } from "tests/common/constants";
import {
  describeIfDeployed,
  goToPage,
  login,
  tryUntil,
} from "tests/utils/helpers";
import { getTestID, getText } from "tests/utils/selectors";
import datasets from "../fixtures/datasets";

const TEST_COLLECTION = {
  contactEmail: "TEST@example.com",
  contactName: "TEST NAME",
  description: "TEST DESCRIPTION",
  name: "TEST COLLECTION",
};

describe("Collection", () => {
  describeIfDeployed("Logged In Tests", () => {
    it("creates and deletes a collection", async () => {
      const timestamp = Date.now();
      await login();

      const collectionName = "TEST_COLLECTION" + timestamp;

      await createCollection({ name: collectionName });

      // Try delete
      await page.click(getTestID("collection-more-button"));
      await page.click(getText("Delete Collection"));

      await Promise.all([
        page.waitForNavigation({ waitUntil: "load" }),
        page.click(".bp3-alert-footer >> text=Delete Collection"),
      ]);

      await tryUntil(async () => {
        await expect(page).not.toHaveSelector(getText(collectionName));
      }, 50);
    });

    describe("dataset order", () => {
      let lastValue = 1_000_000_000;
      sortByCellCountDescending(datasets).forEach((dataset) => {
        expect(dataset.cell_count).toBeLessThanOrEqual(lastValue);
        lastValue = dataset.cell_count ?? 0;
      });
    });

    describe("Publish a collection", () => {
      describe("when no dataset", () => {
        it("shows disabled publish button", async () => {
          await login();

          await createCollection();

          await tryUntil(async () => {
            const publishButton = await page.$(
              getTestID("publish-collection-button")
            );

            expect(await publishButton?.getAttribute("disabled")).toBe("");
          }, 100);
        });
      });
    });
  });
});

async function createCollection(
  collection?: Partial<Collection>
): Promise<string> {
  await goToPage(`${TEST_URL}${ROUTES.MY_COLLECTIONS}`);

  await page.click(getText("Create Collection"));

  const testCollection = { ...TEST_COLLECTION, ...collection };

  await page.type("#name", testCollection.name, BLUEPRINT_SAFE_TYPE_OPTIONS);
  await page.type(
    "#description",
    testCollection.description,
    BLUEPRINT_SAFE_TYPE_OPTIONS
  );
  await page.type(
    "#contact-name",
    testCollection.contactName,
    BLUEPRINT_SAFE_TYPE_OPTIONS
  );
  await page.type(
    "#contact-email",
    testCollection.contactEmail,
    BLUEPRINT_SAFE_TYPE_OPTIONS
  );

  const [response] = await Promise.all([
    page.waitForEvent("response"),
    page.click(getTestID("create-button")),
  ]);

  const { collection_uuid } = (await response.json()) as {
    collection_uuid: string;
  };

  await expect(page).toHaveSelector(getText(testCollection.name));

  return collection_uuid;
}
