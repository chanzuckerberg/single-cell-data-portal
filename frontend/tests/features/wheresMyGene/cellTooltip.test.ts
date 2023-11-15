import { expect } from "@playwright/test";
import { expandTissue, getCellTypeNames, goToPage } from "tests/utils/helpers";
import { TEST_URL } from "tests/common/constants";
import { ROUTES } from "src/common/constants/routes";
import {
  CELL_TYPE_NAME_LABEL_CLASS_NAME,
  CELL_TYPE_ROW_CLASS_NAME,
} from "src/views/WheresMyGeneV2/components/HeatMap/components/YAxisChart/constants";
import { test } from "tests/common/test";
import assert from "assert";
import {
  NO_MARKER_GENES_DESCRIPTION,
  NO_MARKER_GENES_FOR_BLOOD_DESCRIPTION,
  TOO_FEW_CELLS_NO_MARKER_GENES_DESCRIPTION,
} from "src/views/WheresMyGeneV2/components/CellInfoSideBar/constants";
const { describe } = test;

const NO_MARKER_GENES_DESCRIPTION_ID = "no-marker-genes-description";

describe("cell tooltip", () => {
  test(`Should verify cell tooltip hover`, async ({ page }) => {
    await goToPage(`${TEST_URL}${ROUTES.WHERE_IS_MY_GENE}`, page);

    const cellName = CELL_TYPE_NAME_LABEL_CLASS_NAME;

    // (thuang): Expand blood tissue to find truncated cell names
    await expandTissue(page, "blood");

    const cellsAfter = await getCellTypeNames(page);

    // get the number of cell names displayed

    let checkedCells = 0;

    // no need to check all truncated cell names
    for (let i = 0; i < cellsAfter.length && checkedCells < 2; i++) {
      // skip cells that are not truncated
      if (!cellsAfter[i].includes("...")) continue;

      // get the text displayed when user hovers over truncated cell names
      const hoverText =
        (await page
          .getByTestId(CELL_TYPE_ROW_CLASS_NAME)
          .nth(i)
          .getByTestId("cell-type-full-name")
          .textContent()) || "no text";

      await page.getByTestId(cellName).nth(i).hover();

      //expect the hover text to be displayed
      await expect(page.getByText(hoverText)).toBeVisible();

      checkedCells++;
    }
  });
  test(`Should verify blood cells have no marker genes`, async ({ page }) => {
    await goToPage(`${TEST_URL}${ROUTES.WHERE_IS_MY_GENE}`, page);

    // Expand blood tissue
    await expandTissue(page, "blood");

    // Click stem cell info icon
    await page.getByTestId("cell-type-info-button-blood-stem cell").click();

    // Verify copy is what we expect
    const noMarkerGenesDescription = (await page
      .getByTestId(NO_MARKER_GENES_DESCRIPTION_ID)
      .textContent()) as string;
    assert.strictEqual(
      noMarkerGenesDescription.trim(),
      NO_MARKER_GENES_FOR_BLOOD_DESCRIPTION
    );
  });
  test(`Should verify cell types with < 25 cells have no marker genes`, async ({
    page,
  }) => {
    await goToPage(`${TEST_URL}${ROUTES.WHERE_IS_MY_GENE}`, page);

    // Expand blood tissue
    await expandTissue(page, "adipose-tissue");

    // Click naive B cell info icon
    await page
      .getByTestId("cell-type-info-button-adipose tissue-naive B cell")
      .click();

    // Verify copy is what we expect
    const noMarkerGenesDescription = (await page
      .getByTestId(NO_MARKER_GENES_DESCRIPTION_ID)
      .textContent()) as string;
    assert.strictEqual(
      noMarkerGenesDescription.trim(),
      TOO_FEW_CELLS_NO_MARKER_GENES_DESCRIPTION
    );
  });
  test(`Should verify copy for cell types with no marker genes`, async ({
    page,
  }) => {
    await goToPage(`${TEST_URL}${ROUTES.WHERE_IS_MY_GENE}`, page);

    // Expand blood tissue
    await expandTissue(page, "yolk-sac");

    // Click yolk sac somatic cell info icon
    await page
      .getByTestId("cell-type-info-button-yolk sac-somatic cell")
      .click();

    // Verify copy is what we expect
    const noMarkerGenesDescription = (await page
      .getByTestId(NO_MARKER_GENES_DESCRIPTION_ID)
      .textContent()) as string;
    assert.strictEqual(
      noMarkerGenesDescription.trim(),
      NO_MARKER_GENES_DESCRIPTION
    );
  });
});
