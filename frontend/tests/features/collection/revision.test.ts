import { expect, Locator, Page, test } from "@playwright/test";
import { ROUTES } from "src/common/constants/routes";
import { Collection } from "src/common/entities";
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
      `*/**/dp/v1/collections/${collectionId}`,
      async (route, request) => {
        // Handle GET collection requests.
        if (request.method() === "GET") {
          const response = await route.fetch();
          const json: Collection = await response.json();

          // Modify name of each dataset to fake diff between revision and
          // published counterpart.
          json.datasets?.forEach((d) => {
            d.name = `${Math.random()}`;
          });

          await route.fulfill({ response, json });
        } else {
          // We're not expecting POST (create revision) or DELETE (delete
          // revision) requests at this point; handling these cases for
          // completion but they will result in a failing test.
          await route.continue();
        }
      },
      { times: 1 }
    );

    // Reload page required to force re-fetch of collection and published
    // counterpart (via the useCollection hook) rather than using the (React
    // Query) cached values.
    await page.reload();
    await page.waitForURL(url);

    // We have faked changes in datasets; publish button should be enabled.
    await tryUntil(
      async () => {
        const publishButton = await page.$(
          getTestID("publish-collection-button")
        );
        await expect(publishButton).toBeEnabled();
      },
      { page }
    );

    await deleteRevision(page);
  });

  /**
   * TODO(#5666): Enable this test once #5666 is resolved
   * https://app.zenhub.com/workspaces/single-cell-5e2a191dad828d52cc78b028/issues/gh/chanzuckerberg/single-cell-data-portal/5666
   */
  test("starts a revision", async ({ page }) => {
    const testId = buildCollectionRowLocator(COLLECTION_ROW_WRITE_PUBLISHED_ID);
    const collectionName = await startRevision(page, testId);

    const publishButton = await page.$(getTestID("publish-collection-button"));

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
    // so upping this for avoid flakiness
    const RETRY_TIMES = 100;

    await tryUntil(
      async () => {
        await expect(collectionRowContinueLocator).toBeVisible();
      },
      { maxRetry: RETRY_TIMES, page }
    );

    const viewRevisionButton = await collectionRowContinueLocator.getByTestId(
      COLLECTION_VIEW_REVISION
    );

    await expect(viewRevisionButton).toBeVisible();

    await viewRevisionButton?.click();

    await deleteRevision(page);
  });

  /**
   * TODO(#5666): Enable this test once #5666 is resolved
   * https://app.zenhub.com/workspaces/single-cell-5e2a191dad828d52cc78b028/issues/gh/chanzuckerberg/single-cell-data-portal/5666
   */

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

    const collectionContactEmail = (
      await page.getAttribute(getTestID(COLLECTION_CONTACT_ID), "href")
    )?.replace(/^mailto:/, "");

    await getCollectionMoreButtonLocator(page).click();

    await page.getByTestId("dropdown-edit-details").click();

    const COLLECTION_CONTENT_ID = "collection-form-content";

    // Assert that the edit form is visible with the right input values
    await page.waitForSelector(getTestID(COLLECTION_CONTENT_ID));

    expect(await page.inputValue("input#name")).toBe(collectionName);
    expect(await page.inputValue("textarea#description")).toBe(
      collectionDescription
    );
    expect(await page.inputValue("input#contact-name")).toBe(
      collectionContactName
    );
    expect(await page.inputValue("input#contact-email")).toBe(
      collectionContactEmail
    );

    const newCollectionName = String(Date.now());

    await page.fill("input#name", newCollectionName);

    await page.getByTestId("create-button").click();

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

    expect(
      (
        await page.getAttribute(getTestID(COLLECTION_CONTACT_ID), "href")
      )?.replace(/^mailto:/, "")
    ).toBe(collectionContactEmail);

    const publishButton = await page.$(getTestID("publish-collection-button"));

    await expect(publishButton).toBeEnabled();

    await deleteRevision(page);
  });

  describe("reorder datasets", () => {
    const MIN_DATASET_COUNT = 3;

    test("enables reorder datasets", async ({ page }) => {
      // Navigate to revision with more than one dataset.
      await startReorderableRevision(page);

      // Confirm reorder menu option is available.
      await getCollectionMoreButtonLocator(page).click();
      await expect(
        page.getByTestId(TEST_ID_DROPDOWN_REORDER_DATASETS)
      ).toBeVisible();

      // Tear down.
      await deleteRevision(page);
    });

    test("cancels reorder", async ({ page }) => {
      // Navigate to revision with more than one dataset.
      await startReorderableRevision(page);

      // Enable reorder.
      await enterReorderMode(page);

      // Confirm cancel is available.
      await expect(
        page.getByTestId(TEST_ID_REORDER_DATASETS_CANCEL)
      ).toBeVisible();

      // Click cancel to return back to edit mode.
      await cancelReorder(page);

      // Confirm cancel button is no longer visible.
      await expect(
        page.getByTestId(TEST_ID_REORDER_DATASETS_CANCEL)
      ).not.toBeVisible();

      // Confirm publish button is visible.
      await expect(page.getByTestId(TEST_ID_PUBLISH_COLLECTION)).toBeVisible();

      // Tear dfown.
      await deleteRevision(page);
    });

    test("saves reorder", async ({ page }) => {
      // Navigate to revision with more than one dataset.
      await startReorderableRevision(page);

      // Enable reorder.
      await enterReorderMode(page);

      // Grab the datasets - we need at least three datasets for this test.
      const datasets = await locateDatasets(page);
      expect(datasets.length).toBeGreaterThanOrEqual(MIN_DATASET_COUNT);

      // Get a reference to the current dataset order; we'll use this to assert reorder is successful.
      const datasetIds = await listDatasetIds(datasets);

      // Drag first dataset and drop at second position.
      const fromIndex = 0;
      const toIndex = 1;
      await moveDatasetTo(page, datasets, fromIndex, toIndex);

      // Confirm save is available.
      await expect(
        page.getByTestId(TEST_ID_REORDER_DATASETS_SAVE)
      ).toBeVisible();

      // Click save to save order and to return back to edit mode.
      await page.getByTestId(TEST_ID_REORDER_DATASETS_SAVE).click();

      // Confirm save button is no longer visible.
      await tryUntil(
        async () => {
          await expect(
            page.getByTestId(TEST_ID_REORDER_DATASETS_SAVE)
          ).not.toBeVisible();
        },
        { page }
      );

      // Confirm publish button is visible.
      await tryUntil(
        async () => {
          await expect(
            page.getByTestId(TEST_ID_PUBLISH_COLLECTION)
          ).toBeVisible();
        },
        { page }
      );

      // Confirm post-drop order of datasets has been retained.
      const reorderedDatasetIds = await listDatasetIds(datasets);
      const expectedDatasetIds = reorderDatasetIds(
        datasetIds,
        fromIndex,
        toIndex
      );
      expect(reorderedDatasetIds).toEqual(expectedDatasetIds);

      // Tear down.
      await deleteRevision(page);
    });

    test("move first dataset to second position", async ({ page }) => {
      // Navigate to revision with more than one dataset.
      await startReorderableRevision(page);

      // Enable reorder.
      await enterReorderMode(page);

      // Grab the datasets - we need at least three datasets for this test.
      const datasets = await locateDatasets(page);
      expect(datasets.length).toBeGreaterThanOrEqual(MIN_DATASET_COUNT);

      // Get a reference to the current dataset order; we'll use this to assert reorder is successful.
      const datasetIds = await listDatasetIds(datasets);

      // Drag first dataset and drop at second position.
      const fromIndex = 0;
      const toIndex = 1;
      await moveDatasetTo(page, datasets, fromIndex, toIndex);

      // Confirm post-drop order of datasets has been retained.
      const reorderedDatasetIds = await listDatasetIds(datasets);
      const expectedDatasetIds = reorderDatasetIds(
        datasetIds,
        fromIndex,
        toIndex
      );
      expect(reorderedDatasetIds).toEqual(expectedDatasetIds);

      // Tear down.
      await cancelReorder(page);
      await deleteRevision(page);
    });

    test("move first dataset to last position", async ({ page }) => {
      // Navigate to revision with more than one dataset.
      await startReorderableRevision(page);

      // Enable reorder.
      await enterReorderMode(page);

      // Grab the datasets - we need at least three datasets for this test.
      const datasets = await locateDatasets(page);
      const datasetCount = datasets.length;
      expect(datasetCount).toBeGreaterThanOrEqual(MIN_DATASET_COUNT);

      // Get a reference to the current dataset order; we'll use this to assert reorder is successful.
      const datasetIds = await listDatasetIds(datasets);

      // Drag first dataset and drop at last position.
      const fromIndex = 0;
      const toIndex = datasetCount - 1;
      await moveDatasetTo(page, datasets, fromIndex, toIndex);

      // Confirm post-drop order of datasets has been retained.
      const reorderedDatasetIds = await listDatasetIds(datasets);
      const expectedDatasetIds = reorderDatasetIds(
        datasetIds,
        fromIndex,
        toIndex
      );
      expect(reorderedDatasetIds).toEqual(expectedDatasetIds);

      // Tear down.
      await cancelReorder(page);
      await deleteRevision(page);
    });

    test("move last dataset to first position", async ({ page }) => {
      // Navigate to revision with more than one dataset.
      await startReorderableRevision(page);

      // Enable reorder.
      await enterReorderMode(page);

      // Grab the datasets - we need at least three datasets for this test.
      const datasets = await locateDatasets(page);
      const datasetCount = datasets.length;
      expect(datasetCount).toBeGreaterThanOrEqual(MIN_DATASET_COUNT);

      // Get a reference to the current dataset order; we'll use this to assert reorder is successful.
      const datasetIds = await listDatasetIds(datasets);

      // Drag last dataset and drop at fist position.
      const fromIndex = datasetCount - 1;
      const toIndex = 0;
      await moveDatasetTo(page, datasets, fromIndex, toIndex);

      // Confirm post-drop order of datasets has been retained.
      const reorderedDatasetIds = await listDatasetIds(datasets);
      const expectedDatasetIds = reorderDatasetIds(
        datasetIds,
        fromIndex,
        toIndex
      );
      expect(reorderedDatasetIds).toEqual(expectedDatasetIds);

      // Tear down.
      await cancelReorder(page);
      await deleteRevision(page);
    });

    test.skip("move last dataset to second position", async ({ page }) => {
      // Navigate to revision with more than one dataset.
      await startReorderableRevision(page);

      // Enable reorder.
      await enterReorderMode(page);

      // Grab the datasets - we need at least three datasets for this test.
      const datasets = await locateDatasets(page);
      const datasetCount = datasets.length;
      expect(datasetCount).toBeGreaterThanOrEqual(MIN_DATASET_COUNT);

      // Get a reference to the current dataset order; we'll use this to assert reorder is successful.
      const datasetIds = await listDatasetIds(datasets);

      // Drag last dataset and drop at second position.
      const fromIndex = datasetCount - 1;
      const toIndex = 1;
      await moveDatasetTo(page, datasets, fromIndex, toIndex);

      // Confirm post-drop order of datasets has been retained.
      const reorderedDatasetIds = await listDatasetIds(datasets);
      const expectedDatasetIds = reorderDatasetIds(
        datasetIds,
        fromIndex,
        toIndex
      );
      expect(reorderedDatasetIds).toEqual(expectedDatasetIds);

      // Tear down.
      await cancelReorder(page);
      await deleteRevision(page);
    });
  });
});

/**
 * (thuang): Wait for 1 min instead of the default 3 minutes, so we fail faster
 */
const WAIT_FOR_MIN_USABLE_COLLECTION_TIMEOUT_MS = 1 * 60 * 1000; // 1 minute

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

        await tryUntil(
          async () => {
            const collectionRows = await page.getByTestId(collectionRowTestId);
            expect(await collectionRows.count()).toBeGreaterThan(
              MIN_USABLE_COLLECTION_COUNT - 1
            );
          },
          { page }
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

  // Simulate drag over rows between from row and to row to trigger drag and drop calculations correctly.
  // For example, if we're dragging from 0 to 3, we want to drag over 1 then 2. If we're dragging from 3 to
  // 0, we want to drag over 2 then 1.
  const dragOverIndices =
    fromIndex > toIndex
      ? Array.from(
          { length: fromIndex - toIndex - 1 },
          (_, i) => fromIndex - i - 1
        )
      : Array.from(
          { length: toIndex - fromIndex - 1 },
          (_, i) => i + fromIndex + 1
        );
  await Promise.all(
    dragOverIndices.map(async (i) => {
      return await datasets[i].hover();
    })
  );

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
