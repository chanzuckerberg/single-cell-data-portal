import { ROUTES } from "src/common/constants/routes";
import { Collection } from "src/common/entities";
import { sortByCellCountDescending } from "src/components/Collection/components/CollectionDatasetsGrid/components/DatasetsGrid";
import { INVALID_DOI_ERROR_MESSAGE } from "src/components/CreateCollectionModal/components/Content";
import { BLUEPRINT_SAFE_TYPE_OPTIONS, TEST_URL } from "tests/common/constants";
import {
  describeIfDeployed,
  describeIfDevStaging,
  goToPage,
  login,
  tryUntil,
} from "tests/utils/helpers";
import { getTestID, getText } from "tests/utils/selectors";
import datasets from "../fixtures/datasets";

const TEST_COLLECTION: CollectionFormInput = {
  contact_email: "TEST@example.com",
  contact_name: "TEST NAME",
  description: "TEST DESCRIPTION",
  name: "TEST COLLECTION",
};

/**
 * HTML element ID of DOI input field.
 */
const ELEMENT_ID_INPUT_DOI = "#DOI";

/**
 * Subset of collection fields required for creating/editing a collection.
 */
type CollectionFormInput = Pick<
  Collection,
  "contact_email" | "contact_name" | "description" | "name"
>;

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
        page.click(".bp4-alert-footer >> text=Delete Collection"),
      ]);

      await tryUntil(async () => {
        await expect(page).not.toHaveSelector(getText(collectionName));
      }, 50);
    });

    test("dataset order", () => {
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

  describeIfDevStaging("Deployed Env Tests", () => {
    describe("invalid DOIs", () => {
      it("doesn't create a collection with a DOI in an invalid format", async () => {
        const timestamp = Date.now();

        await login();
        await showCreateForm();
        await populateRequiredInputs({
          ...TEST_COLLECTION,
          name: "TEST_COLLECTION" + timestamp,
        });

        // Specify a DOI that is in an invalid format.
        await populatePublicationDOI("10.1016/j.2022.104097");

        // Attempt submit, confirm error message is displayed.
        const [response] = await submitCreateForm();
        expect(response.status()).toEqual(400);
        await expect(page).toHaveSelector(getText(INVALID_DOI_ERROR_MESSAGE));
      });

      it("doesn't create a collection with an invalid DOI", async () => {
        const timestamp = Date.now();

        await login();
        await showCreateForm();
        await populateRequiredInputs({
          ...TEST_COLLECTION,
          name: "TEST_COLLECTION" + timestamp,
        });

        // Specify a DOI that is a valid format but is not on Crossref.
        await populatePublicationDOI("x");

        // Attempt submit, confirm error message is displayed.
        const [response] = await submitCreateForm();
        expect(response.status()).toEqual(400);
        await expect(page).toHaveSelector(getText(INVALID_DOI_ERROR_MESSAGE));
      });
    });
  });
});

async function createCollection(
  collection?: Partial<Collection>
): Promise<string> {
  await showCreateForm();

  const testCollection = { ...TEST_COLLECTION, ...collection };

  await populateRequiredInputs(testCollection);

  const [response] = await submitCreateForm();

  const { collection_id } = (await response.json()) as {
    collection_id: string;
  };

  await expect(page).toHaveSelector(getText(testCollection.name));

  return collection_id;
}

/**
 * Display the collection form modal.
 */
async function showCreateForm() {
  await goToPage(`${TEST_URL}${ROUTES.MY_COLLECTIONS}`);
  await page.click(getText("Create Collection"));
}

/**
 * Submit create collection form.
 * @returns Form submit response.
 */
async function submitCreateForm() {
  return await Promise.all([
    page.waitForEvent("response"),
    page.click(getTestID("create-button")),
  ]);
}

/**
 * Specify a publication DOI on the collection form.
 * @param value - Value to enter in the DOI input field.
 */
async function populatePublicationDOI(value: string) {
  await page.click(getText("Add Link"));
  await page.click(getText("Publication DOI"));
  await expect(page).toHaveSelector(
    getText(
      "A summary citation linked to this DOI will be automatically added to this collection."
    )
  );
  await page.type(ELEMENT_ID_INPUT_DOI, value, BLUEPRINT_SAFE_TYPE_OPTIONS);
}

/**
 * Populate required values on collection form.
 * @param testCollection - Collection form input values.
 */
async function populateRequiredInputs(testCollection: CollectionFormInput) {
  await page.type("#name", testCollection.name, BLUEPRINT_SAFE_TYPE_OPTIONS);
  await page.type(
    "#description",
    testCollection.description,
    BLUEPRINT_SAFE_TYPE_OPTIONS
  );
  await page.type(
    "#contact-name",
    testCollection.contact_name,
    BLUEPRINT_SAFE_TYPE_OPTIONS
  );
  await page.type(
    "#contact-email",
    testCollection.contact_email,
    BLUEPRINT_SAFE_TYPE_OPTIONS
  );
}
