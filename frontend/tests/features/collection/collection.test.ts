import { expect, Page, test } from "@playwright/test";
import HTTP_STATUS_CODE from "src/common/constants/HTTP_STATUS_CODE";
import { ROUTES } from "src/common/constants/routes";
import { Collection } from "src/common/entities";
import { INVALID_DOI_ERROR_MESSAGE } from "src/components/CreateCollectionModal/components/Content/common/constants";
import { TEST_URL } from "tests/common/constants";
import { goToPage, isDevStagingRdev, tryUntil } from "tests/utils/helpers";
import { getTestID } from "tests/utils/selectors";

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
  describe("Logged In Tests @loggedIn", () => {
    skip(
      !isDevStagingRdev,
      `Currently push-test runs against dev BE, so login doesn't work for local containers.
       We also only run the tests in dev, rdev, and staging to avoid polluting prod.
      `
    );

    test("creates and deletes a collection", async ({ page }) => {
      const timestamp = Date.now();
      const collectionName = "TEST_COLLECTION" + timestamp;

      await createCollection({ collection: { name: collectionName }, page });

      // Try delete
      await page.getByTestId("collection-more-button").click();
      await page.getByText("Delete Collection").click();

      // Wait for the page to navigate away after confirming deletion!
      const currentURL = page.url();

      await Promise.all([
        page.waitForURL((url) => url.toString() !== currentURL),
        page.click(".bp5-alert-footer >> text=Delete Collection"),
      ]);

      await tryUntil(
        async () => {
          const isVisible = await page
            .getByText(collectionName, { exact: true })
            .isVisible()
            .catch(() => false);
          expect(isVisible).toBeFalsy();
        },
        { page }
      );
    });

    describe("Publish a collection", () => {
      describe("when no dataset", () => {
        test("shows disabled publish button", async ({ page }) => {
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

  describe("Deployed Env Tests @loggedIn", () => {
    skip(
      !isDevStagingRdev,
      `BE DOI endpoints only work in dev, staging, and prod. And we also only run
      the tests in dev, rdev, and staging to avoid polluting prod.
      `
    );

    describe("invalid DOIs", () => {
      test("doesn't create a collection with a DOI in an invalid format", async ({
        page,
      }) => {
        const timestamp = Date.now();

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

        // Attempt submit and confirm 400 response
        const [response] = await submitCreateFormInvalid(page);
        expect(response.status()).toEqual(400);

        // Wait for the specific error message to be visible
        await expect(page.getByText(INVALID_DOI_ERROR_MESSAGE)).toBeVisible({
          timeout: 6000,
        });
      });

      test("doesn't create a collection with an invalid DOI", async ({
        page,
      }) => {
        const timestamp = Date.now();

        await showCreateForm(page);
        await populateRequiredInputs(
          {
            ...TEST_COLLECTION,
            name: "TEST_COLLECTION" + timestamp,
          },
          page
        );

        // Specify a DOI that has a valid format but does not exist on Crossref
        const VALID_FORMAT_BUT_NON_EXISTENT_DOI = "10.1016/j.2022.104097";
        await populatePublicationDOI(VALID_FORMAT_BUT_NON_EXISTENT_DOI, page);

        // Submit and confirm error
        const [response] = await submitCreateFormInvalid(page);
        expect(response.status()).toEqual(400);

        // Wait for the error message to be visible
        await expect(page.getByText(INVALID_DOI_ERROR_MESSAGE)).toBeVisible({
          timeout: 6000,
        });
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

  expect(page.getByText(testCollection.name)).toBeTruthy();

  return collection_id;
}

/**
 * Display the collection form modal.
 */
async function showCreateForm(page: Page) {
  await goToPage(`${TEST_URL}${ROUTES.COLLECTIONS}`, page);
  await page.getByText("Create Collection").click();
}

const collectionEndpointPath = `/dp/v1/collections`;

/**
 * Submit create invalid collection form.
 * @returns Form submit invalid parameters error response (400).
 */
async function submitCreateFormInvalid(page: Page) {
  return await Promise.all([
    page.waitForResponse(
      (response) =>
        response.url().includes(collectionEndpointPath) &&
        response.status() === HTTP_STATUS_CODE.BAD_REQUEST
    ),
    page.getByTestId("create-button").click(),
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
        response.url().includes(collectionEndpointPath) &&
        response.status() === HTTP_STATUS_CODE.OK
    ),
    page.getByTestId("create-button").click(),
  ]);
}

/**
 * Specify a publication DOI on the collection form.
 * @param value - Value to enter in the DOI input field.
 */
async function populatePublicationDOI(value: string, page: Page) {
  await page.getByText("Add Link").click();
  await page.getByText("Publication DOI").click();

  await expect(
    page.getByText(
      "A summary citation linked to this DOI will be automatically added to this collection."
    )
  ).toBeVisible();

  await page.locator(ELEMENT_ID_INPUT_DOI).fill(value);
}

/**
 * Populate required values on collection form.
 * @param testCollection - Collection form input values.
 */
async function populateRequiredInputs(
  testCollection: CollectionFormInput,
  page: Page
) {
  await page.locator("#name").fill(testCollection.name);
  await page.locator("#name").press("Tab");

  await page.locator("#description").fill(testCollection.description);
  await page.locator("#description").press("Tab");

  await page.locator("#contact-name").fill(testCollection.contact_name);
  await page.locator("#contact-name").press("Tab");

  await page.locator("#contact-email").fill(testCollection.contact_email);
  await page.locator("#contact-email").press("Tab");
}
