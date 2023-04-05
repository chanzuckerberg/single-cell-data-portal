/**
 * Test suite for select filter-related utils.
 */
import { expect, test } from "@playwright/test";
import { TEST_URL } from "tests/common/constants";

const { describe } = test;
const CHEVRON_LEFT = '[data-icon="chevron-left"]';
const FIRST_OPTION = '[data-option-index="0"]';
const DATASET_FILTER = "dataset-filter";

describe("Left side bar", () => {
  test.beforeEach(async ({ page }) => {
    // navigate to url
    await page.goto(TEST_URL + "/gene-expression", { timeout: 60000 });
  });

  test("Left side bar collapse and expand", async ({ page }) => {
    // click cheevron left to collaspe the left tab
    await page.locator(CHEVRON_LEFT).click();

    //verify the left tab is collapsed

    expect(await page.getByTestId(DATASET_FILTER).isVisible()).toBeFalsy();
  });
  test("Should be able select and de-select options for datasetfilter", async ({
    page,
  }) => {
    // click the dataset filter
    await page
      .getByTestId(DATASET_FILTER)

      .click({ position: { x: 0, y: 0 } });

    // select the first dataset
    await page.locator(FIRST_OPTION).click();
    await page.getByTestId(DATASET_FILTER).click({ position: { x: 0, y: 0 } });

    const filter_label = '[data-testid="dataset-filter"] [role="button"]';
    //expect the  selected filter to be visible
    await expect(page.locator(filter_label)).toBeVisible();

    //click the cancel button
    await page.getByTestId("ClearIcon").click();
    const visibility = await page.locator(filter_label).isVisible();
    expect(visibility).toBeFalsy();
  });
  test("Should be able select and de-select options for disease filter", async ({
    page,
  }) => {
    // click the dataset filter
    await page
      .getByTestId("disease-filter")
      .click({ position: { x: 0, y: 0 } });

    // select the first dataset
    await page.locator(FIRST_OPTION).click();

    //click filter again to close pop up
    await page
      .getByTestId("disease-filter")
      .click({ position: { x: 0, y: 0 } });
    const filter_label = '[data-testid="disease-filter"] [role="button"]';
    //expect the  selected filter to be visible
    await expect(page.locator(filter_label)).toBeVisible();

    //click the cancel button
    await page.getByTestId("ClearIcon").click();
    const visibility = await page.locator(filter_label).isVisible();
    expect(visibility).toBeFalsy();
  });
  test("Should be able select and de-select options for Self-Reported Ethnicity filter", async ({
    page,
  }) => {
    // click the dataset filter
    await page
      .getByTestId("self-reported-ethnicity-filter")
      .click({ position: { x: 0, y: 0 } });

    // select the first dataset
    await page.locator(FIRST_OPTION).click();

    //click filter again to close pop up
    await page
      .getByTestId("self-reported-ethnicity-filter")
      .click({ position: { x: 0, y: 0 } });
    const filter_label =
      '[data-testid="self-reported-ethnicity-filter"] [role="button"]';
    //expect the  selected filter to be visible
    await expect(page.locator(filter_label)).toBeVisible();

    //click the cancel button
    await page.getByTestId("ClearIcon").click();
    const visibility = await page.locator(filter_label).isVisible();
    expect(visibility).toBeFalsy();
  });
  test("Should be able select and de-select options for Self-Reported sex-filter", async ({
    page,
  }) => {
    // click the dataset filter
    await page.getByTestId("sex-filter").click({ position: { x: 0, y: 0 } });

    // select the first dataset
    await page.locator(FIRST_OPTION).click();

    //click filter again to close pop up
    await page.getByTestId("sex-filter").click({ position: { x: 0, y: 0 } });
    const filter_label = '[data-testid="sex-filter"] [role="button"]';
    //expect the  selected filter to be visible
    await expect(page.locator(filter_label)).toBeVisible();

    //click the cancel button
    await page.getByTestId("ClearIcon").click();
    const visibility = await page.locator(filter_label).isVisible();
    expect(visibility).toBeFalsy();
  });
});
