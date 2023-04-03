import { ElementHandle, expect, Page, test } from "@playwright/test";
import { ROUTES } from "src/common/constants/routes";
import { TEST_URL } from "tests/common/constants";
import {
  getInnerText,
  goToPage,
  isDevStaging,
  tryUntil,
} from "tests/utils/helpers";
import { getTestID, getText } from "tests/utils/selectors";

const { describe, skip } = test;

const COLLECTION_ROW_ID = "collection-row";

describe("Collection Revision", () => {
  skip(
    !isDevStaging,
    "We only seed published collections for revision test in dev and staging"
  );

  test("starts a revision", async ({ page }) => {
    const collectionName = await startRevision(page);

    const publishButton = await page.$(getTestID("publish-collection-button"));

    await expect(publishButton).toBeDisabled();

    await page.getByText("My Collections").click();

    const COLLECTION_ROW_SELECTOR = `*${page.getByTestId(
      COLLECTION_ROW_ID
    )} >> text="${collectionName}"`;

    let collectionRowContinue: ElementHandle | null = null;

    // (thuang): Staging is slow due to the amount of collections we fetch,
    // so upping this for avoid flakiness
    const RETRY_TIMES = 100;

    await tryUntil(
      async () => {
        const collections = await page.$$(COLLECTION_ROW_SELECTOR);

        for (const collection of collections) {
          const hasPublished = Boolean(
            await collection.$(getText("Published"))
          );
          const hasRevisionPending = Boolean(
            await collection.$(getTestID("Revision Pending"))
          );

          if (hasPublished && hasRevisionPending) {
            collectionRowContinue = collection;
            break;
          }
        }

        if (!collectionRowContinue) {
          throw new Error("Collection not found");
        }
      },
      { maxRetry: RETRY_TIMES, page }
    );

    expect(collectionRowContinue).not.toBe(null);

    const actionButtonContinue = await (
      collectionRowContinue as unknown as ElementHandle<HTMLElement>
    ).$(getTestID("revision-action-button"));

    await expect(actionButtonContinue).toMatchText("Continue");

    await actionButtonContinue?.click();

    await deleteRevision(page);
  });

  test("allows editing", async ({ page }) => {
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

    await page.getByTestId("collection-more-button").click();

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
        expect(page.getByTestId(COLLECTION_CONTENT_ID)).toBeFalsy();
      },
      { page }
    );

    await expect(page).toMatchText(new RegExp(newCollectionName));
    await expect(page).toMatchText(new RegExp(collectionDescription));
    await expect(page).toMatchText(new RegExp(collectionContactName));

    const REVISION_STATUS_TEXT =
      "This collection has changed since you last published it";

    await expect(page).toMatchText(new RegExp(REVISION_STATUS_TEXT));

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
  await goToPage(TEST_URL + ROUTES.MY_COLLECTIONS, page);
  // (thuang): Wait for collections to load to prevent race condition
  await page.waitForLoadState("networkidle");

  const MIN_USABLE_COLLECTION_COUNT = 4;

  // (thuang): If we can't find at least 4 usable collections, we'll delete a revision
  await tryUntil(
    async () => {
      try {
        await expect(page.getByText("Start Revision")).toBeTruthy();

        await tryUntil(
          async () => {
            const collectionRows = await page.locator(
              getText("Start Revision")
            );
            expect(await collectionRows.count()).toBeGreaterThan(
              MIN_USABLE_COLLECTION_COUNT - 1
            );
          },
          { page }
        );
      } catch {
        await page.getByText("Continue").click();
        await deleteRevision(page);
        throw new Error("No available collection");
      }
    },
    { page }
  );

  /**
   * (thuang): NOTE: the `*` at the beginning of the string captures the
   * `COLLECTION_ROW_ID` selector element
   * @see https://playwright.dev/docs/selectors#intermediate-matches
   */
  const COLLECTION_ROW_SELECTOR_START = `*${page.getByTestId(
    COLLECTION_ROW_ID
  )} >> text="Start Revision"`;

  // (thuang): We randomly select a collection row to start a revision
  const collectionRows = await page.$$(COLLECTION_ROW_SELECTOR_START);
  const collectionRowCount = await collectionRows.length;

  const collectionRow =
    collectionRows[Math.floor(Math.random() * collectionRowCount)];

  expect(collectionRow).not.toBe(null);

  const VISIBILITY_TAG_ID = "visibility-tag";

  const visibilityTag = await collectionRow?.$(getTestID(VISIBILITY_TAG_ID));

  expect(visibilityTag).not.toBe(null);

  await expect(visibilityTag).toMatchText("Published");

  const collectionName = await collectionRow?.$eval(
    getTestID("collection-link"),
    (element: HTMLElement) => element.textContent
  );

  const actionButton = await collectionRow?.$(
    getTestID("revision-action-button")
  );

  await actionButton?.click();

  const PRIVATE_REVISION_TEXT =
    "This is a private revision of a public collection.";

  await tryUntil(
    async () => await expect(page).toHaveSelector(getTestID("revision-status")),
    { page }
  );

  const revisionStatusElement = await page.$(getTestID("revision-status"));

  await expect(revisionStatusElement).toMatchText(PRIVATE_REVISION_TEXT);

  return collectionName as string;
}

async function deleteRevision(page: Page) {
  const DROPDOWN_CANCEL_ID = "dropdown-cancel-revision";

  if (!(await page.$(getTestID(DROPDOWN_CANCEL_ID)))) {
    await page.getByTestId("collection-more-button").click();
  }

  await page.getByTestId(DROPDOWN_CANCEL_ID).click();

  await page.locator(".bp4-alert-footer >> text=Cancel Revision").click();

  await page.waitForURL(TEST_URL + ROUTES.MY_COLLECTIONS);
}
