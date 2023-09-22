import { expect, Page, test } from "@playwright/test";
import { ROUTES } from "src/common/constants/routes";
import { TEST_URL } from "tests/common/constants";
import {
  getInnerText,
  goToPage,
  isDevStaging,
  tryUntil,
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

describe("Collection Revision @loggedIn", () => {
  skip(
    !isDevStaging,
    "We only seed published collections for revision test in dev and staging"
  );

  /**
   * TODO(#5666): Enable this test once #5666 is resolved
   * https://app.zenhub.com/workspaces/single-cell-5e2a191dad828d52cc78b028/issues/gh/chanzuckerberg/single-cell-data-portal/5666
   */
  test.skip("starts a revision", async ({ page }) => {
    const collectionName = await startRevision(page);

    const publishButton = await page.$(getTestID("publish-collection-button"));

    await expect(publishButton).toBeDisabled();

    await page.getByTestId(COLLECTIONS_LINK_ID).click();
    await page.waitForURL("**" + ROUTES.COLLECTIONS);

    const collectionRowContinueLocator = page
      .getByTestId(COLLECTION_ROW_WRITE_REVISION_ID)
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
  test.skip("allows editing", async ({ page }) => {
    await startRevision(page);

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
});

async function startRevision(page: Page): Promise<string> {
  const MIN_USABLE_COLLECTION_COUNT = 4;

  // (thuang): If we can't find at least 4 usable collections, we'll delete a revision
  await tryUntil(
    async () => {
      await goToPage(TEST_URL + ROUTES.COLLECTIONS, page);

      try {
        await expect(
          page.getByTestId(COLLECTION_ROW_WRITE_PUBLISHED_ID)
        ).toBeTruthy();

        await tryUntil(
          async () => {
            const collectionRows = await page.getByTestId(
              COLLECTION_ROW_WRITE_PUBLISHED_ID
            );
            expect(await collectionRows.count()).toBeGreaterThan(
              MIN_USABLE_COLLECTION_COUNT - 1
            );
          },
          { page }
        );
      } catch {
        await page
          .getByTestId(COLLECTION_ROW_WRITE_REVISION_ID)
          .first()
          .getByTestId(COLLECTION_VIEW_REVISION)
          .click();
        await deleteRevision(page);
        throw new Error("No available collection");
      }
    },
    { page }
  );

  // (thuang): We randomly select a collection row to start a revision
  const collectionRows = await page.getByTestId(
    COLLECTION_ROW_WRITE_PUBLISHED_ID
  );
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
    .locator(".bp4-alert-footer")
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
