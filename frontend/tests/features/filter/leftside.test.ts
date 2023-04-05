/**
 * Test suite for select filter-related utils.
 */
import { expect, test } from "@playwright/test";
import { TEST_URL } from "tests/common/constants";

const { describe } = test;
const chevron_left = '[data-icon="chevron-left"]';
const first_option = '[data-option-index="0"]';

describe("Left side bar", () => {
  test.beforeEach(async ({ page }) => {
    // navigate to url
    await page.goto(TEST_URL + "/gene-expression", { timeout: 60000 });
  });

  test("Left side bar collapse and expand", async ({ page }) => {
    // click cheevron left to collaspe the left tab
    await page.locator(chevron_left).click();

    //verify the left tab is collapsed
    // eslint-disable-next-line sonarjs/no-duplicate-string
    expect(await page.getByTestId("dataset-filter").isVisible()).toBeFalsy();
  });
  test("Should be able select and de-select options for datasetfilter", async ({
    page,
  }) => {
    // click the dataset filter
    await page
      .getByTestId("dataset-filter")
      .click({ position: { x: 0, y: 0 } });

    // select the first dataset
    await page.locator(first_option).click();
    await page
      .getByTestId("dataset-filter")
      .click({ position: { x: 0, y: 0 } });

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
    await page.locator(first_option).click();

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
    await page.locator(first_option).click();

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
    await page.locator(first_option).click();

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
