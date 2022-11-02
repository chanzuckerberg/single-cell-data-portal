import { expect, Page, test } from "@playwright/test";
import HTTP_STATUS_CODE from "src/common/constants/HTTP_STATUS_CODE";
import { ROUTES } from "src/common/constants/routes";
import { Collection } from "src/common/entities";
import { sortByCellCountDescending } from "src/components/Collection/components/CollectionDatasetsGrid/components/DatasetsGrid/common/util";
import { INVALID_DOI_ERROR_MESSAGE } from "src/components/CreateCollectionModal/components/Content/common/constants";
import { API_URL } from "src/configs/configs";
import { BLUEPRINT_SAFE_TYPE_OPTIONS, TEST_URL } from "tests/common/constants";
import {
  goToPage,
  isDevStagingProd,
  login,
  tryUntil,
} from "tests/utils/helpers";
import { getTestID, getText } from "tests/utils/selectors";
import datasets from "../fixtures/datasets";

const { describe, skip } = test;

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
  describe("Logged In Tests", () => {
    skip(
      !isDevStagingProd,
      "Currently push-test runs against dev BE, so login doesn't work for local containers"
    );

    test("creates and deletes a collection", async ({ page }) => {
      const timestamp = Date.now();
      await login(page);

      const collectionName = "TEST_COLLECTION" + timestamp;

      await createCollection({ collection: { name: collectionName }, page });

      // Try delete
      await page.click(getTestID("collection-more-button"));
      await page.click(getText("Delete Collection"));

      await Promise.all([
        page.waitForNavigation({ waitUntil: "load" }),
        page.click(".bp4-alert-footer >> text=Delete Collection"),
      ]);

      await tryUntil(
        async () => {
          await expect(page).not.toHaveSelector(getText(collectionName));
        },
        { page }
      );
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
        test("shows disabled publish button", async ({ page }) => {
          await login(page);

          await createCollection({ page });

          await tryUntil(
            async () => {
              const publishButton = await page.$(
                getTestID("publish-collection-button")
              );

              expect(await publishButton?.getAttribute("disabled")).toBe("");
            },
            { maxRetry: 100, page }
          );
        });
      });
    });
  });

  describe("Deployed Env Tests", () => {
    skip(!isDevStagingProd, "BE DOI endpoints only work in dev, staging, prod");

    describe("invalid DOIs", () => {
      test("doesn't create a collection with a DOI in an invalid format", async ({
        page,
      }) => {
        const timestamp = Date.now();

        await login(page);
        await showCreateForm(page);
        await populateRequiredInputs(
          {
            ...TEST_COLLECTION,
            name: "TEST_COLLECTION" + timestamp,
          },
          page
        );

        // Specify a DOI that is in an invalid format.
        await populatePublicationDOI("INVALID_FORMAT", page);

        // Attempt submit, confirm error message is displayed.
        const [response] = await submitCreateFormInvalid(page);
        expect(response.status()).toEqual(400);
        await expect(page).toHaveSelector(getText(INVALID_DOI_ERROR_MESSAGE));
      });

      test("doesn't create a collection with an invalid DOI", async ({
        page,
      }) => {
        const timestamp = Date.now();

        await login(page);
        await showCreateForm(page);
        await populateRequiredInputs(
          {
            ...TEST_COLLECTION,
            name: "TEST_COLLECTION" + timestamp,
          },
          page
        );

        // Specify a DOI that is a valid format but is not on Crossref.
        const VALID_FORMAT_BUT_NON_EXISTENT_DOI = "10.1016/j.2022.104097";
        await populatePublicationDOI(VALID_FORMAT_BUT_NON_EXISTENT_DOI, page);

        // Attempt submit, confirm error message is displayed.
        const [response] = await submitCreateFormInvalid(page);
        expect(response.status()).toEqual(400);
        await expect(page).toHaveSelector(getText(INVALID_DOI_ERROR_MESSAGE));
      });
    });
  });
});

async function createCollection({
  collection,
  page,
}: {
  collection?: Partial<Collection>;
  page: Page;
}): Promise<string> {
  await showCreateForm(page);

  const testCollection = { ...TEST_COLLECTION, ...collection };

  await populateRequiredInputs(testCollection, page);

  const [response] = await submitCreateForm(page);

  const { collection_id } = (await response.json()) as {
    collection_id: string;
  };

  await expect(page).toHaveSelector(getText(testCollection.name));

  return collection_id;
}

/**
 * Display the collection form modal.
 */
async function showCreateForm(page: Page) {
  await goToPage(`${TEST_URL}${ROUTES.MY_COLLECTIONS}`, page);
  await page.click(getText("Create Collection"));
}

const collectionEndpoint = `${API_URL}/dp/v1/collections`;

/**
 * Submit create invalid collection form.
 * @returns Form submit invalid parameters error response (400).
 */
async function submitCreateFormInvalid(page: Page) {
  return await Promise.all([
    page.waitForResponse(
      (response) =>
        response.url() === collectionEndpoint &&
        response.status() === HTTP_STATUS_CODE.BAD_REQUEST
    ),
    page.click(getTestID("create-button")),
  ]);
}

/**
 * Submit create collection form.
 * @returns Form submit response.
 */
async function submitCreateForm(page: Page) {
  return await Promise.all([
    page.waitForResponse(
      (response) =>
        response.url() === collectionEndpoint &&
        response.status() === HTTP_STATUS_CODE.OK
    ),
    page.click(getTestID("create-button")),
  ]);
}

/**
 * Specify a publication DOI on the collection form.
 * @param value - Value to enter in the DOI input field.
 */
async function populatePublicationDOI(value: string, page: Page) {
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
async function populateRequiredInputs(
  testCollection: CollectionFormInput,
  page: Page
) {
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
