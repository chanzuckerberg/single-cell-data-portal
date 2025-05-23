import { expect, Locator, Page, test } from "@playwright/test";
import { ROUTES } from "src/common/constants/routes";
import { Collection } from "src/common/entities";
import {
  DATASET_EDIT_FORM,
  DATASET_EDIT_SAVE,
} from "src/components/Collection/components/CollectionDatasetsGrid/components/Row/DatsetRow/components/EditDataset/constants";
import { DROPDOWN_EDIT_DATASET } from "src/components/Collection/components/CollectionDatasetsGrid/components/Row/DatsetRow/components/MoreDropdown/components/Menu/constants";
import { DATASET_MORE_BUTTON } from "src/components/Collection/components/CollectionDatasetsGrid/components/Row/DatsetRow/components/MoreDropdown/constants";
import { DATASET_TITLE } from "src/components/Datasets/components/Grid/components/DatasetNameCell/constants";
import { TEST_URL } from "tests/common/constants";
import {
  getInnerText,
  goToPage,
  isDevStagingRdev,
  tryUntil,
  waitForElement,
} from "tests/utils/helpers";
import { getTestID } from "tests/utils/selectors";

const { describe, skip } = test;

const COLLECTION_ACTIONS_ID = "collection-actions";
const COLLECTION_NAME_ID = "collection-link";
const COLLECTION_STATUS_BANNER_ID = "revision-status";
const COLLECTION_VIEW_REVISION = "view-collection-revision";
const COLLECTION_ROW_WRITE_PUBLISHED_ID = "collection-row-write-published";
const COLLECTION_ROW_WRITE_REVISION_ID =
  "collection-row-write-published-revision";
const COLLECTIONS_LINK_ID = "collections-link";
const STATUS_TAG_ID = "status-tag";
const TEST_ID_DROPDOWN_REORDER_DATASETS = "dropdown-reorder-datasets";
const TEST_ID_PUBLISH_COLLECTION = "publish-collection-button";
const TEST_ID_REORDER_DATASETS_CANCEL = "datasets-reorder-cancel";
const TEST_ID_REORDER_DATASETS_SAVE = "datasets-reorder-save";

describe("Collection Revision @loggedIn", () => {
  skip(
    !isDevStagingRdev,
    "We only seed published collections for revision test in dev, rdev, and staging"
  );

  test("enables publish if datasets updated", async ({ page }) => {
    const testId = buildCollectionRowLocator(COLLECTION_ROW_WRITE_PUBLISHED_ID);
    await startRevision(page, testId);

    const url = await page.url();
    const collectionId = url.split("/").pop();

    // Fake unique datasets in get collection response.
    await page.route(
      `**/dp/v1/collections/${collectionId}`,
      async (route, request) => {
        if (request.method() === "GET") {
          const originalResponse = await route.fetch();
          const json: Collection = await originalResponse.json();

          // Modify dataset names to simulate a difference
          json.datasets?.forEach((d) => {
            d.name = `fake-dataset-${Math.random()}`;
          });

          await route.fulfill({
            status: 200,
            contentType: "application/json",
            body: JSON.stringify(json),
          });
        } else {
          await route.continue();
        }
      },
      { times: 1 }
    );

    // Reload the page to force data refetch
    await page.reload();
    await page.waitForURL(url);

    // Wait until the publish button becomes enabled
    await tryUntil(
      async () => {
        const publishButton = page.getByTestId(TEST_ID_PUBLISH_COLLECTION);
        await expect(publishButton).toBeEnabled();
      },
      { page }
    );

    await deleteRevision(page);
  });

  test("starts a revision", async ({ page }) => {
    const testId = buildCollectionRowLocator(COLLECTION_ROW_WRITE_PUBLISHED_ID);
    const collectionName = await startRevision(page, testId);

    const publishButton = page.getByTestId(TEST_ID_PUBLISH_COLLECTION);
    await expect(publishButton).toBeDisabled();

    await page.getByTestId(COLLECTIONS_LINK_ID).click();
    await page.waitForURL("**" + ROUTES.COLLECTIONS);

    const collectionRowContinueLocator = page
      .getByTestId(buildCollectionRowLocator(COLLECTION_ROW_WRITE_REVISION_ID))
      .filter({ hasText: collectionName })
      .filter({ hasText: "Published" })
      .filter({ hasText: "Revision" })
      .first();

    // (thuang): Staging is slow due to the amount of collections we fetch,
    // so upping this to avoid flakiness
    const RETRY_TIMES = 150;

    await tryUntil(
      async () => {
        await expect(collectionRowContinueLocator).toBeVisible();
      },
      { maxRetry: RETRY_TIMES, page }
    );

    const viewRevisionButton = collectionRowContinueLocator.getByTestId(
      COLLECTION_VIEW_REVISION
    );
    await expect(viewRevisionButton).toBeVisible();

    await viewRevisionButton.click();

    await deleteRevision(page);
  });

  test("allows editing", async ({ page }) => {
    const testId = buildCollectionRowLocator(COLLECTION_ROW_WRITE_PUBLISHED_ID);
    await startRevision(page, testId);

    const collectionName = await getInnerText(
      getTestID("collection-name"),
      page
    );
    const collectionDescription = await getInnerText(
      getTestID("collection-description"),
      page
    );

    const COLLECTION_CONTACT_ID = "collection-contact";
    const collectionContactName = await getInnerText(
      getTestID(COLLECTION_CONTACT_ID),
      page
    );

    const rawEmail = await page.getAttribute(
      getTestID(COLLECTION_CONTACT_ID),
      "href"
    );
    const collectionContactEmail = rawEmail?.replace(/^mailto:/, "");

    if (!collectionContactEmail) {
      throw new Error("Expected collectionContactEmail to be defined");
    }

    // Open the edit details modal
    await getCollectionMoreButtonLocator(page).click();
    await page.getByTestId("dropdown-edit-details").click();

    const COLLECTION_CONTENT_ID = "collection-form-content";

    // Ensure edit form is visible and fields are pre-filled correctly
    await expect(page.getByTestId(COLLECTION_CONTENT_ID)).toBeVisible();

    await expect(page.locator("input#name")).toHaveValue(collectionName);
    await expect(page.locator("textarea#description")).toHaveValue(
      collectionDescription
    );
    await expect(page.locator("input#contact-name")).toHaveValue(
      collectionContactName
    );
    await expect(page.locator("input#contact-email")).toHaveValue(
      collectionContactEmail
    );

    // Update the collection name
    const newCollectionName = String(Date.now());
    await page.fill("input#name", newCollectionName);

    // Save changes
    await page.getByTestId("create-button").click();

    // Verify the form closes and the new data appears on the page
    await tryUntil(
      async () => {
        await expect(page.getByTestId(COLLECTION_CONTENT_ID)).not.toBeVisible();

        await expect(
          page.getByText(new RegExp(newCollectionName))
        ).toBeVisible();
        await expect(
          page.getByText(new RegExp(collectionDescription))
        ).toBeVisible();
        await expect(
          page.getByText(new RegExp(collectionContactName))
        ).toBeVisible();
      },
      { page }
    );

    // Verify email is unchanged
    const updatedEmail = (
      await page.getAttribute(getTestID(COLLECTION_CONTACT_ID), "href")
    )?.replace(/^mailto:/, "");

    expect(updatedEmail).toBe(collectionContactEmail);

    // Confirm publish is now enabled
    const publishButton = page.getByTestId(TEST_ID_PUBLISH_COLLECTION);
    await expect(publishButton).toBeEnabled();

    await deleteRevision(page);
  });

  describe("edit dataset", () => {
    test("allows rename dataset", async ({ page }) => {
      // Create and navigate to revision
      const testId = buildCollectionRowLocator(
        COLLECTION_ROW_WRITE_PUBLISHED_ID
      );
      await startRevision(page, testId);

      try {
        // Get first dataset row in collection
        const datasetRows = await locateDatasets(page);
        expect(datasetRows.length).toBeGreaterThanOrEqual(1);
        const datasetRow = datasetRows[0];

        // Show the dataset more menu
        const datasetMoreButton = datasetRow.getByTestId(DATASET_MORE_BUTTON);
        await expect(datasetMoreButton).toBeVisible();
        await datasetMoreButton.click();

        // Show the edit dataset modal
        const editDatasetButton = page.getByTestId(DROPDOWN_EDIT_DATASET);
        await expect(editDatasetButton).toBeVisible();
        await editDatasetButton.click();

        // Confirm modal is visible
        const editForm = page.getByTestId(DATASET_EDIT_FORM);
        await expect(editForm).toBeVisible();

        // Confirm input shows current dataset name
        const datasetName = await getInnerText(
          `css=[data-testid^="${DATASET_TITLE}"]`,
          page
        );
        await expect(page.locator("input#title")).toHaveValue(datasetName);

        // Enter new title
        const newDatasetTitle = String(Date.now());
        await page.fill("input#title", newDatasetTitle);

        // Save
        await page.getByTestId(DATASET_EDIT_SAVE).click();

        // Confirm modal closes and progress bar appears
        await tryUntil(
          async () => {
            await expect(page.getByTestId(DATASET_EDIT_FORM)).toBeHidden();
            await expect(page.getByRole("progressbar")).toBeVisible();
          },
          { page }
        );
      } finally {
        // Always clean up
        await deleteRevision(page);
      }
    });
  });

  describe.skip("reorder datasets", () => {
    const MIN_DATASET_COUNT = 3;

    test("enables reorder datasets", async ({ page }) => {
      // Navigate to a revision with multiple datasets
      await startReorderableRevision(page);

      try {
        // Open collection actions menu
        await getCollectionMoreButtonLocator(page).click();

        // Confirm reorder option is visible
        const reorderOption = page.getByTestId(
          TEST_ID_DROPDOWN_REORDER_DATASETS
        );
        await expect(reorderOption).toBeVisible();
      } finally {
        // Always clean up
        await deleteRevision(page);
      }
    });

    test("cancels reorder", async ({ page }) => {
      // Navigate to revision with multiple datasets
      await startReorderableRevision(page);

      try {
        // Enter reorder mode
        await enterReorderMode(page);

        // Confirm cancel button is visible
        const cancelButton = page.getByTestId(TEST_ID_REORDER_DATASETS_CANCEL);
        await expect(cancelButton).toBeVisible();

        // Click cancel to exit reorder mode
        await cancelReorder(page);

        // Confirm cancel button is no longer visible
        await expect(cancelButton).not.toBeVisible();

        // Confirm publish button is visible (and optionally enabled)
        const publishButton = page.getByTestId(TEST_ID_PUBLISH_COLLECTION);
        await expect(publishButton).toBeVisible();
      } finally {
        // Always clean up
        await deleteRevision(page);
      }
    });

    test("saves reorder", async ({ page }) => {
      await startReorderableRevision(page);

      try {
        await enterReorderMode(page);

        const datasets = await locateDatasets(page);
        expect(datasets.length).toBeGreaterThanOrEqual(MIN_DATASET_COUNT);

        const datasetIds = await listDatasetIds(datasets);

        const fromIndex = 0;
        const toIndex = 1;
        await moveDatasetTo(page, datasets, fromIndex, toIndex);

        const reorderSaveButton = page.getByTestId(
          TEST_ID_REORDER_DATASETS_SAVE
        );
        await expect(reorderSaveButton).toBeVisible();

        await reorderSaveButton.click();

        // Wait for save to disappear
        await tryUntil(
          async () => {
            await expect(reorderSaveButton).not.toBeVisible();
          },
          { page }
        );

        const publishButton = page.getByTestId(TEST_ID_PUBLISH_COLLECTION);
        await tryUntil(
          async () => {
            await expect(publishButton).toBeVisible();
          },
          { page }
        );

        // Re-fetch dataset rows after reorder to reflect updated DOM
        const reorderedDatasets = await locateDatasets(page);
        const reorderedDatasetIds = await listDatasetIds(reorderedDatasets);
        const expectedDatasetIds = reorderDatasetIds(
          datasetIds,
          fromIndex,
          toIndex
        );

        expect(reorderedDatasetIds).toEqual(expectedDatasetIds);
      } finally {
        await deleteRevision(page);
      }
    });

    test("move first dataset to second position", async ({ page }) => {
      await startReorderableRevision(page);

      try {
        await enterReorderMode(page);

        const datasets = await locateDatasets(page);
        expect(datasets.length).toBeGreaterThanOrEqual(MIN_DATASET_COUNT);

        const datasetIds = await listDatasetIds(datasets);

        const fromIndex = 0;
        const toIndex = 1;
        await moveDatasetTo(page, datasets, fromIndex, toIndex);

        // Re-fetch datasets in case DOM was updated after reordering
        const updatedDatasets = await locateDatasets(page);
        const reorderedDatasetIds = await listDatasetIds(updatedDatasets);

        const expectedDatasetIds = reorderDatasetIds(
          datasetIds,
          fromIndex,
          toIndex
        );
        expect(reorderedDatasetIds).toEqual(expectedDatasetIds);
      } finally {
        await cancelReorder(page);
        await deleteRevision(page);
      }
    });

    test("move first dataset to last position", async ({ page }) => {
      // Navigate to a revision with more than one dataset
      await startReorderableRevision(page);

      try {
        // Enable reorder mode
        await enterReorderMode(page);

        // Grab the datasets — we need at least three datasets for this test
        const initialDatasets = await locateDatasets(page);
        const datasetCount = initialDatasets.length;
        expect(datasetCount).toBeGreaterThanOrEqual(MIN_DATASET_COUNT);

        // Get a reference to the current dataset order;
        // we'll use this to assert reorder was successful
        const initialDatasetIds = await listDatasetIds(initialDatasets);

        // Drag the first dataset and drop it at the last position
        const fromIndex = 0;
        const toIndex = datasetCount - 1;
        await moveDatasetTo(page, initialDatasets, fromIndex, toIndex);

        // Re-fetch dataset list to reflect any DOM updates after reorder
        const updatedDatasets = await locateDatasets(page);
        const reorderedDatasetIds = await listDatasetIds(updatedDatasets);

        // Compute the expected order and assert correctness
        const expectedDatasetIds = reorderDatasetIds(
          initialDatasetIds,
          fromIndex,
          toIndex
        );
        expect(reorderedDatasetIds).toEqual(expectedDatasetIds);
      } finally {
        // Exit reorder mode and clean up revision
        await cancelReorder(page);
        await deleteRevision(page);
      }
    });

    test("move last dataset to first position", async ({ page }) => {
      // Navigate to a revision with more than one dataset
      await startReorderableRevision(page);

      try {
        // Enable reorder mode
        await enterReorderMode(page);

        // Grab the datasets — we need at least three datasets for this test
        const initialDatasets = await locateDatasets(page);
        const datasetCount = initialDatasets.length;
        expect(datasetCount).toBeGreaterThanOrEqual(MIN_DATASET_COUNT);

        // Get a reference to the current dataset order;
        // we'll use this to assert reorder is successful
        const initialDatasetIds = await listDatasetIds(initialDatasets);

        // Drag the last dataset and drop it at the first position
        const fromIndex = datasetCount - 1;
        const toIndex = 0;
        await moveDatasetTo(page, initialDatasets, fromIndex, toIndex);

        // Re-fetch datasets to avoid using stale DOM references
        const updatedDatasets = await locateDatasets(page);
        const reorderedDatasetIds = await listDatasetIds(updatedDatasets);

        // Generate expected order and confirm the reorder was successful
        const expectedDatasetIds = reorderDatasetIds(
          initialDatasetIds,
          fromIndex,
          toIndex
        );
        expect(reorderedDatasetIds).toEqual(expectedDatasetIds);
      } finally {
        // Exit reorder mode and clean up the revision
        await cancelReorder(page);
        await deleteRevision(page);
      }
    });

    test("move last dataset to second position", async ({ page }) => {
      // Navigate to a revision with more than one dataset
      await startReorderableRevision(page);

      try {
        // Enable reorder mode
        await enterReorderMode(page);

        // Grab the datasets — we need at least three datasets for this test
        const initialDatasets = await locateDatasets(page);
        const datasetCount = initialDatasets.length;
        expect(datasetCount).toBeGreaterThanOrEqual(MIN_DATASET_COUNT);

        // Get the current dataset order; we'll use this to assert the reorder result
        const initialDatasetIds = await listDatasetIds(initialDatasets);

        // Drag the last dataset and drop it at the second position (index 1)
        const fromIndex = datasetCount - 1;
        const toIndex = 1;
        await moveDatasetTo(page, initialDatasets, fromIndex, toIndex);

        // Re-fetch dataset rows to avoid stale references
        const updatedDatasets = await locateDatasets(page);
        const reorderedDatasetIds = await listDatasetIds(updatedDatasets);

        // Generate the expected order and assert correctness
        const expectedDatasetIds = reorderDatasetIds(
          initialDatasetIds,
          fromIndex,
          toIndex
        );
        expect(reorderedDatasetIds).toEqual(expectedDatasetIds);
      } finally {
        // Exit reorder mode and clean up the revision
        await cancelReorder(page);
        await deleteRevision(page);
      }
    });
  });
});

/**
 * (thuang): Wait for 1 min instead of the default 3 minutes, so we fail faster
 */
const WAIT_FOR_MIN_USABLE_COLLECTION_TIMEOUT_MS = 3 * 60 * 1000; // 1 minute

async function startRevision(
  page: Page,
  collectionRowTestId: RegExp
): Promise<string> {
  const MIN_USABLE_COLLECTION_COUNT = 4;
  // (thuang): If we can't find at least 4 usable collections, we'll delete a revision
  await tryUntil(
    async () => {
      await goToPage(TEST_URL + ROUTES.COLLECTIONS, page);

      try {
        await expect(page.getByTestId(collectionRowTestId)).toBeTruthy();

        const collectionRows = await page.getByTestId(collectionRowTestId);

        expect(await collectionRows.count()).toBeGreaterThan(
          MIN_USABLE_COLLECTION_COUNT - 1
        );
      } catch {
        const revisionRowTestId = buildCollectionRowLocator(
          COLLECTION_ROW_WRITE_REVISION_ID
        );
        const hasAnyRevision =
          (await page.getByTestId(revisionRowTestId).count()) > 0;

        if (hasAnyRevision) {
          await page
            .getByTestId(revisionRowTestId)
            .first()
            .getByTestId(COLLECTION_VIEW_REVISION)
            .click();
          await deleteRevision(page);
        }

        throw new Error("No available collection");
      }
    },
    { page, timeoutMs: WAIT_FOR_MIN_USABLE_COLLECTION_TIMEOUT_MS }
  );

  // (thuang): We randomly select a collection row to start a revision
  const collectionRows = await page.getByTestId(collectionRowTestId);
  const collectionRowCount = await collectionRows.count();
  const randomIndex = Math.floor(Math.random() * collectionRowCount);
  const collectionRow = collectionRows.nth(randomIndex);

  expect(collectionRow).not.toBe(null);

  const statusTag = await collectionRow?.getByTestId(STATUS_TAG_ID);

  // We expect to find one status tag, with the text "Published".
  expect(statusTag).not.toBe(null);
  await expect(statusTag).toHaveCount(1);
  await expect(statusTag.getByText("Published")).toBeVisible();

  // Grab the collection link.
  const collectionLink = await collectionRow?.getByTestId(COLLECTION_NAME_ID);
  // Grab the collection name.
  const collectionName = await collectionLink.innerText();
  // Go to collection page.
  await collectionLink?.click();

  await tryUntil(
    async () =>
      await expect(page.getByTestId(COLLECTION_ACTIONS_ID)).toBeVisible(),
    { page }
  );

  // Start revision of collection
  await page
    .getByTestId(COLLECTION_ACTIONS_ID)
    .locator("button")
    .filter({ hasText: "Start Revision" })
    .click();

  const PRIVATE_REVISION_TEXT =
    "This is a private revision of a published collection. Open Published Collection";

  await tryUntil(
    async () => {
      await expect(page.getByTestId(COLLECTION_STATUS_BANNER_ID)).toBeVisible();

      const collectionStatusBanner = page.getByTestId(
        COLLECTION_STATUS_BANNER_ID
      );

      await expect(
        collectionStatusBanner.getByText(PRIVATE_REVISION_TEXT)
      ).toBeVisible();
    },

    { page }
  );

  return collectionName;
}

async function deleteRevision(page: Page) {
  const DROPDOWN_CANCEL_ID = "dropdown-cancel-revision";

  // Grab the published collection route from the collection status banner.
  const collectionRoute = await page
    .getByTestId(COLLECTION_STATUS_BANNER_ID)
    .locator("a")
    .getAttribute("href");

  /**
   * (thuang): Sometimes the dropdown is already open, so we need to check if it's
   * visible before clicking on the "More" button
   */
  if (!(await page.$(getTestID(DROPDOWN_CANCEL_ID)))) {
    await getCollectionMoreButtonLocator(page).click();
    await tryUntil(
      async () => {
        await expect(page.getByTestId(DROPDOWN_CANCEL_ID)).toBeVisible();
      },
      { page }
    );
  }

  await page.getByTestId(DROPDOWN_CANCEL_ID).click();

  await page
    .locator(".bp5-alert-footer")
    .locator("button")
    .filter({ hasText: "Cancel Revision" })
    .click();

  // The page routes to published collection when the revision is successfully deleted.
  await page.waitForURL(`${TEST_URL}${collectionRoute}`);
}

function getCollectionMoreButtonLocator(page: Page) {
  return page
    .getByTestId("collection-actions")
    .getByTestId("collection-more-button");
}

/**
 * Build a locator from the given stem, matches on strings starting with the stem followed
 * by -d where d is any number.
 * @param testIdStem - Base string to build locator from.
 * @returns Regular expression for finding collections or revisions.
 */
function buildCollectionRowLocator(testIdStem: string): RegExp {
  return new RegExp(`^${testIdStem}-\\d+$`);
}

/**
 * Build a locator from the given stem, specific to finding collections with more than three
 * datasets. Matches on strings starting with the stem followed by -d where d is any number greater
 * than three.
 * @param testIdStem - Base string to build locator from.
 * @returns Regular expression for finding collections or revisions.
 */
function buildReorderableCollectionRowLocator(): RegExp {
  return new RegExp(`^${COLLECTION_ROW_WRITE_PUBLISHED_ID}-([3-9]|\\d{2,})$`);
}

/**
 * Cancel reorder mode by clicking the "cancel reorder" button.
 * @param page - Playwright page model.
 */
async function cancelReorder(page: Page): Promise<void> {
  // Click cancel to return back to edit mode.
  await page.getByTestId(TEST_ID_REORDER_DATASETS_CANCEL).click();
}

/**
 * Enter reorder mode by clicking the "reorder datasets" menu option.
 * @param page - Playwright page model.
 */
async function enterReorderMode(page: Page): Promise<void> {
  await getCollectionMoreButtonLocator(page).click();
  await page.getByTestId(TEST_ID_DROPDOWN_REORDER_DATASETS).click();
}

/**
 * Grab the set of IDs for the given datasets.
 * @param datasets - Array of locators representing dataset rows.
 * @returns Array of dataset IDs.
 */
async function listDatasetIds(datasets: Locator[]): Promise<(string | null)[]> {
  return await Promise.all(
    datasets.map(async (dataset) => await dataset.getAttribute("data-testid"))
  );
}

/**
 * Find all datasets (that is, dataset rows) on the page.
 * @param page - Playwright page model.
 */
async function locateDatasets(page: Page): Promise<Locator[]> {
  return await page.getByTestId(/^dataset-row-/).all();
}

/**
 * Calculate the updated dataset ID by moving the dataset ID from the fromIndex to the toIndex.
 * @param datasetIds - Array of dataset IDs in the pre-reorder order.
 * @param fromIndex - Index of dataset to move.
 * @param toIndex - Index of where to move dataset to.
 * @returns Array of dataset IDs updated to match new order.
 */
function reorderDatasetIds(
  datasetIds: (string | null)[],
  fromIndex: number,
  toIndex: number
): (string | null)[] {
  const reorderedIds = [...datasetIds];
  reorderedIds.splice(toIndex, 0, ...reorderedIds.splice(fromIndex, 1));
  return reorderedIds;
}

/**
 * Move dataset at from index to the given to index.
 * @param page - Playwright page model.
 * @param datasets - Array of locators representing dataset rows.
 * @param fromIndex - Index of dataset to move.
 * @param toIndex - Index of where to move dataset to.
 */
async function moveDatasetTo(
  page: Page,
  datasets: Locator[],
  fromIndex: number,
  toIndex: number
): Promise<void> {
  await datasets[fromIndex].hover();
  await page.mouse.down();

  await datasets[toIndex].hover();
  await page.mouse.up();
}

/**
 * Navigate to a revision that has at least three datasets that can be reordered by
 * the test user.
 * @param page - Playwright page model.
 */
async function startReorderableRevision(page: Page): Promise<void> {
  // Start a revision for a collection that has at least three datasets.
  const testId = buildReorderableCollectionRowLocator();
  await startRevision(page, testId);

  // Enable reorder feature flag; to be removed with #6493.
  const url = await page.url();
  await goToPage(`${url}?reorder=true`, page);
  await waitForElement(page, COLLECTION_ACTIONS_ID);
}
