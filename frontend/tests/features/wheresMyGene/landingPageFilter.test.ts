/**
 * Test suite for select filter-related utils.
 */
import { expect, test } from "@playwright/test";
import { goToWMG } from "tests/common/utils";
import { isDevStagingProd } from "tests/utils/helpers";
import { getTestID } from "tests/utils/selectors";
const CHEVRON_LEFT = '[data-icon="chevron-left"]';
const FIRST_OPTION = '[data-option-index="0"]';
const DATASET_FILTER = "dataset-filter";
const DISEASE_FILTER = "disease-filter";
const SELF_REPORTED_ETHNICITY_FILTER = "self-reported-ethnicity-filter";
const SEX_FILTER = "sex-filter";
const { describe, skip } = test;

describe("Left side bar", () => {
  skip(!isDevStagingProd, "WMG BE API does not work locally or in rdev");
  test.beforeEach(async ({ page }) => {
    // navigate to url
    await goToWMG(page);
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

    //close popup
    await page.keyboard.press("Escape");

    const filter_label = `${getTestID(DATASET_FILTER)} [role="button"]`;
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
    await page.keyboard.press("Escape");

    const filter_label = `${getTestID(DISEASE_FILTER)} [role="button"]`;
    //expect the  selected filter to be visible\

    await expect(page.locator(filter_label)).toBeVisible();

    //click the cancel button
    await page.getByTestId("ClearIcon").click();
    const visibility = await page.locator(filter_label).isVisible();
    expect(visibility).toBeFalsy();
  });
  test.only("Should be able select and de-select options for Self-Reported Ethnicity filter", async ({
    page,
  }) => {
    // click the dataset filter
    await page
      .getByTestId("self-reported-ethnicity-filter")
      .click({ position: { x: 0, y: 0 } });

    // select the first dataset
    await page.locator(FIRST_OPTION).click();

    //close pop up
    await page.keyboard.press("Escape");
    const filter_label = `${getTestID(
      SELF_REPORTED_ETHNICITY_FILTER
    )} [role="button"]`;
    //expect the  selected filter to be visible
    await expect(page.locator(filter_label)).toBeVisible();

    //click the cancel button
    await page.getByTestId("ClearIcon").click();
    const visibility = await page.locator(filter_label).isVisible();
    expect(visibility).toBeFalsy();
  });
  test.only("Should be able select and de-select options for Self-Reported sex-filter", async ({
    page,
  }) => {
    // click the dataset filter
    await page.getByTestId("sex-filter").click({ position: { x: 0, y: 0 } });

    // select the first dataset
    await page.locator(FIRST_OPTION).click();

    // close pop up
    await page.keyboard.press("Escape");
    const filter_label = `${getTestID(SEX_FILTER)} [role="button"]`;
    //expect the  selected filter to be visible
    await expect(page.locator(filter_label)).toBeVisible();

    //click the cancel button
    await page.getByTestId("ClearIcon").click();
    const visibility = await page.locator(filter_label).isVisible();
    expect(visibility).toBeFalsy();
  });
});
