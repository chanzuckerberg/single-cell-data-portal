/**
 * Smoke test suite that will be run in Travis CI
 * Tests included in this file are expected to be relatively stable and test core features
 */

/* eslint-disable no-await-in-loop -- await in loop is needed to emulate sequential user actions  */

import { Classes } from "@blueprintjs/core";
import { appUrlBase, DATASET, DATASET_TRUNCATE } from "./config";

import {
  clickOn,
  goToPage,
  waitByClass,
  waitByID,
  getTestId,
  getTestClass,
  getAllByClass,
  clickOnUntil,
  getOneElementInnerHTML,
  getElementCoordinates,
  tryUntil,
} from "./puppeteerUtils";

import {
  calcDragCoordinates,
  calcTransformDragCoordinates,
  drag,
  scroll,
  expandCategory,
  subset,
  createGeneset,
  deleteGeneset,
  assertGenesetExists,
  assertGenesetDoesNotExist,
  getCellSetCount,
  editGenesetName,
  assertGeneExistsInGeneset,
  removeGene,
  assertGeneDoesNotExist,
  expandGene,
  colorByGeneset,
  assertColorLegendLabel,
  colorByGene,
  clip,
  getAllCategoriesAndCounts,
  selectCategory,
  addGeneToSetAndExpand,
  requestGeneInfo,
  assertGeneInfoCardExists,
  assertGeneInfoCardIsMinimized,
  minimizeGeneInfo,
  removeGeneInfo,
  addGeneToSearch,
  assertGeneInfoDoesNotExist,
  keyboardUndo,
  keyboardRedo,
} from "./cellxgeneActions";

import { datasets } from "./data";

import { scaleMax } from "../../src/util/camera";

const BLUEPRINT_SKELETON_CLASS_NAME = Classes.SKELETON;

// geneset CRUD
const genesetToDeleteName = "geneset_to_delete";
const meanExpressionBrushGenesetName = "second_gene_set";

// const GENES_TO_ADD = ["PF4","PPBP","GNG11","SDPR","NRGN","SPARC","RGS18","TPM4","GP9","GPX1","CD9","TUBB1","ITGA2B"]
// initial text, the text we type in, the result
const editableGenesetName = "geneset_to_edit";
const editText = "_111";
const newGenesetName = "geneset_to_edit_111";

// add gene to set
const geneToAddToSet = "RER1";
const setToAddGeneTo = "fill_this_geneset";

// remove gene from set
const geneToRemove = "SIK1";
const setToRemoveFrom = "empty_this_geneset";

// brush a gene
const geneToBrushAndColorBy = "SIK1";
const brushThisGeneGeneset = "brush_this_gene";
const geneBrushedCellCount = "109";
const subsetGeneBrushedCellCount = "96";

// open gene info card
const geneToRequestInfo = "SIK1";

const genesetDescriptionID =
  "geneset-description-tooltip-fourth_gene_set: fourth description";
const genesetDescriptionString = "fourth_gene_set: fourth description";
const genesetToCheckForDescription = "fourth_gene_set";

const data = datasets[DATASET];
const dataTruncate = datasets[DATASET_TRUNCATE];

const defaultBaseUrl = "d";
const pageUrl = appUrlBase.includes("localhost")
  ? [appUrlBase, defaultBaseUrl, DATASET].join("/")
  : appUrlBase;

const pageUrlTruncate = [appUrlBase, defaultBaseUrl, DATASET_TRUNCATE].join(
  "/",
);

describe("did launch", () => {
  test("page launched", async () => {
    await goToPage(pageUrl);

    const element = await getOneElementInnerHTML(getTestId("header"));

    expect(element).toMatchSnapshot();
  });
});

describe("breadcrumbs loads", () => {
  test("dataset and collection from breadcrumbs appears", async () => {
    await goToPage(pageUrl);

    const datasetElement = await getOneElementInnerHTML(
      getTestId("bc-Dataset"),
    );
    const collectionsElement = await getOneElementInnerHTML(
      getTestId("bc-Collection"),
    );
    expect(datasetElement).toMatchSnapshot();
    expect(collectionsElement).toMatchSnapshot();
  });

  test("datasets from breadcrumbs appears on clicking collections", async () => {
    await goToPage(pageUrl);

    await clickOn(`bc-Dataset`);
    await waitByID("dataset-menu-item-Sed eu nisi condimentum");
    const element = await getOneElementInnerHTML(
      getTestId("dataset-menu-item-Sed eu nisi condimentum"),
    );
    expect(element).toMatchSnapshot();
  });
});

describe("metadata loads", () => {
  test("categories and values from dataset appear", async () => {
    await goToPage(pageUrl);

    for (const label of Object.keys(data.categorical)) {
      const element = await getOneElementInnerHTML(
        getTestId(`category-${label}`),
      );

      expect(element).toMatchSnapshot();

      await clickOn(`${label}:category-expand`);

      const categories = await getAllCategoriesAndCounts(label);

      expect(Object.keys(categories)).toMatchObject(
        // @ts-expect-error ts-migrate(7053) FIXME: Element implicitly has an 'any' type because expre... Remove this comment to see the full error message
        Object.keys(data.categorical[label]),
      );

      expect(Object.values(categories)).toMatchObject(
        // @ts-expect-error ts-migrate(7053) FIXME: Element implicitly has an 'any' type because expre... Remove this comment to see the full error message
        Object.values(data.categorical[label]),
      );
    }
  });

  test("categories and values from dataset appear and properly truncate if applicable", async () => {
    await goToPage(pageUrlTruncate);

    for (const label of Object.keys(dataTruncate.categorical)) {
      const element = await getOneElementInnerHTML(
        getTestId(`category-${label}`),
      );

      expect(element).toMatchSnapshot();

      await clickOn(`${label}:category-expand`);

      const categories = await getAllCategoriesAndCounts(label);

      expect(Object.keys(categories)).toMatchObject(
        // @ts-expect-error ts-migrate(7053) FIXME: Element implicitly has an 'any' type because expre... Remove this comment to see the full error message
        Object.keys(dataTruncate.categorical[label]),
      );

      expect(Object.values(categories)).toMatchObject(
        // @ts-expect-error ts-migrate(7053) FIXME: Element implicitly has an 'any' type because expre... Remove this comment to see the full error message
        Object.values(dataTruncate.categorical[label]),
      );
    }
  });

  test("continuous data appears", async () => {
    await goToPage(pageUrl);

    for (const label of Object.keys(data.continuous)) {
      await waitByID(`histogram-${label}`);
    }
  });
});

describe("cell selection", () => {
  test("selects all cells cellset 1", async () => {
    await goToPage(pageUrl);

    const cellCount = await getCellSetCount(1);
    expect(cellCount).toBe(data.dataframe.nObs);
  });

  test("selects all cells cellset 2", async () => {
    await goToPage(pageUrl);

    const cellCount = await getCellSetCount(2);
    expect(cellCount).toBe(data.dataframe.nObs);
  });

  test("selects cells via lasso", async () => {
    await goToPage(pageUrl);

    for (const cellset of data.cellsets.lasso) {
      const cellset1 = await calcDragCoordinates(
        "layout-graph",
        cellset["coordinates-as-percent"],
      );

      await drag("layout-graph", cellset1.start, cellset1.end, true);
      const cellCount = await getCellSetCount(1);
      expect(cellCount).toBe(cellset.count);
    }
  });

  test("selects cells via categorical", async () => {
    await goToPage(pageUrl);

    for (const cellset of data.cellsets.categorical) {
      await clickOn(`${cellset.metadata}:category-expand`);
      await clickOn(`${cellset.metadata}:category-select`);

      for (const value of cellset.values) {
        await clickOn(`categorical-value-select-${cellset.metadata}-${value}`);
      }

      const cellCount = await getCellSetCount(1);

      expect(cellCount).toBe(cellset.count);
    }
  });

  test("selects cells via continuous", async () => {
    await goToPage(pageUrl);

    for (const cellset of data.cellsets.continuous) {
      const histBrushableAreaId = `histogram-${cellset.metadata}-plot-brushable-area`;

      const coords = await calcDragCoordinates(
        histBrushableAreaId,
        cellset["coordinates-as-percent"],
      );

      await drag(histBrushableAreaId, coords.start, coords.end);

      const cellCount = await getCellSetCount(1);

      expect(cellCount).toBe(cellset.count);
    }
  });
});

describe("subset", () => {
  test("subset - cell count matches", async () => {
    await goToPage(pageUrl);

    for (const select of data.subset.cellset1) {
      if (select.kind === "categorical") {
        await selectCategory(select.metadata, select.values, true);
      }
    }

    await clickOn("subset-button");

    for (const label of Object.keys(data.subset.categorical)) {
      const categories = await getAllCategoriesAndCounts(label);

      expect(Object.keys(categories)).toMatchObject(
        // @ts-expect-error ts-migrate(7053) FIXME: Element implicitly has an 'any' type because expre... Remove this comment to see the full error message
        Object.keys(data.subset.categorical[label]),
      );

      expect(Object.values(categories)).toMatchObject(
        // @ts-expect-error ts-migrate(7053) FIXME: Element implicitly has an 'any' type because expre... Remove this comment to see the full error message
        Object.values(data.subset.categorical[label]),
      );
    }
  });

  test("lasso after subset", async () => {
    await goToPage(pageUrl);

    for (const select of data.subset.cellset1) {
      if (select.kind === "categorical") {
        await selectCategory(select.metadata, select.values, true);
      }
    }

    await clickOn("subset-button");

    const lassoSelection = await calcDragCoordinates(
      "layout-graph",
      data.subset.lasso["coordinates-as-percent"],
    );

    await drag("layout-graph", lassoSelection.start, lassoSelection.end, true);

    const cellCount = await getCellSetCount(1);
    expect(cellCount).toBe(data.subset.lasso.count);
  });
});

describe("clipping", () => {
  test("clip continuous", async () => {
    await goToPage(pageUrl);

    await clip(data.clip.min, data.clip.max);
    const histBrushableAreaId = `histogram-${data.clip.metadata}-plot-brushable-area`;
    const coords = await calcDragCoordinates(
      histBrushableAreaId,
      data.clip["coordinates-as-percent"],
    );
    await drag(histBrushableAreaId, coords.start, coords.end);
    const cellCount = await getCellSetCount(1);
    expect(cellCount).toBe(data.clip.count);

    // ensure categorical data appears properly
    for (const label of Object.keys(data.categorical)) {
      const element = await getOneElementInnerHTML(
        getTestId(`category-${label}`),
      );

      expect(element).toMatchSnapshot();

      await clickOn(`${label}:category-expand`);

      const categories = await getAllCategoriesAndCounts(label);

      expect(Object.keys(categories)).toMatchObject(
        // @ts-expect-error ts-migrate(7053) FIXME: Element implicitly has an 'any' type because expre... Remove this comment to see the full error message
        Object.keys(data.categorical[label]),
      );

      expect(Object.values(categories)).toMatchObject(
        // @ts-expect-error ts-migrate(7053) FIXME: Element implicitly has an 'any' type because expre... Remove this comment to see the full error message
        Object.values(data.categorical[label]),
      );
    }
  });
});

// interact with UI elements just that they do not break
describe("ui elements don't error", () => {
  test("color by", async () => {
    await goToPage(pageUrl);

    const allLabels = [
      ...Object.keys(data.categorical),
      ...Object.keys(data.continuous),
    ];

    for (const label of allLabels) {
      await clickOn(`colorby-${label}`);
    }
  });

  test("pan and zoom", async () => {
    await goToPage(pageUrl);

    await clickOn("mode-pan-zoom");
    const panCoords = await calcDragCoordinates(
      "layout-graph",
      data.pan["coordinates-as-percent"],
    );

    await drag("layout-graph", panCoords.start, panCoords.end, false);

    await page.evaluate("window.scrollBy(0, 1000);");
  });
});

describe("centroid labels", () => {
  test("labels are created", async () => {
    await goToPage(pageUrl);

    const labels = Object.keys(data.categorical);

    await clickOn(`colorby-${labels[0]}`);
    await clickOn("centroid-label-toggle");

    // Toggle colorby for each category and check to see if labels are generated
    for (let i = 0, { length } = labels; i < length; i += 1) {
      const label = labels[i];
      // first label is already enabled
      if (i !== 0) await clickOn(`colorby-${label}`);
      const generatedLabels = await getAllByClass("centroid-label");
      // Number of labels generated should be equal to size of the object
      expect(generatedLabels).toHaveLength(
        // @ts-expect-error ts-migrate(7053) FIXME: Element implicitly has an 'any' type because expre... Remove this comment to see the full error message
        Object.keys(data.categorical[label]).length,
      );
    }
  });
});

describe("graph overlay", () => {
  test("transform centroids correctly", async () => {
    await goToPage(pageUrl);

    const category = Object.keys(data.categorical)[0];

    await clickOn(`colorby-${category}`);
    await clickOn("centroid-label-toggle");
    await clickOn("mode-pan-zoom");

    const panCoords = await calcTransformDragCoordinates(
      "layout-graph",
      data.pan["coordinates-as-percent"],
    );

    // @ts-expect-error ts-migrate(7053) FIXME: Element implicitly has an 'any' type because expre... Remove this comment to see the full error message
    const categoryValue = Object.keys(data.categorical[category])[0];
    const initialCoordinates = await getElementCoordinates(
      `${categoryValue}-centroid-label`,
    );

    await tryUntil(async () => {
      await drag("layout-graph", panCoords.start, panCoords.end, false);

      const terminalCoordinates = await getElementCoordinates(
        `${categoryValue}-centroid-label`,
      );

      expect(terminalCoordinates[0] - initialCoordinates[0]).toBeCloseTo(
        panCoords.end.x - panCoords.start.x,
      );
      expect(terminalCoordinates[1] - initialCoordinates[1]).toBeCloseTo(
        panCoords.end.y - panCoords.start.y,
      );
    });
  });
});

test("zoom limit is 12x", async () => {
  await goToPage(pageUrl);

  const category = Object.keys(data.categorical)[0];

  await clickOn(`colorby-${category}`);
  await clickOn("centroid-label-toggle");
  await clickOn("mode-pan-zoom");

  // @ts-expect-error ts-migrate(7053) FIXME: Element implicitly has an 'any' type because expre... Remove this comment to see the full error message
  const categoryValue = Object.keys(data.categorical[category])[0];
  const initialCoordinates = await getElementCoordinates(
    `${categoryValue}-centroid-label`,
  );

  await tryUntil(async () => {
    await scroll({
      testId: "layout-graph",
      deltaY: -10000,
      coords: initialCoordinates,
    });
    await page.waitForTimeout(1000);
    const newGraph = await page.waitForSelector(
      "[data-test-id^=graph-wrapper-]",
    );
    const newGraphTestId = await newGraph?.evaluate((el) =>
      el.getAttribute("data-test-id"),
    );
    const newDistance = newGraphTestId?.split("distance=").at(-1);
    expect(parseFloat(newDistance)).toBe(scaleMax);
  });
});

test("pan zoom mode resets lasso selection", async () => {
  await goToPage(pageUrl);

  const panzoomLasso = data.features.panzoom.lasso;

  const lassoSelection = await calcDragCoordinates(
    "layout-graph",
    panzoomLasso["coordinates-as-percent"],
  );

  await drag("layout-graph", lassoSelection.start, lassoSelection.end, true);
  await waitByID("lasso-element", { visible: true });

  const initialCount = await getCellSetCount(1);

  expect(initialCount).toBe(panzoomLasso.count);

  await clickOn("mode-pan-zoom");
  await clickOn("mode-lasso");

  const modeSwitchCount = await getCellSetCount(1);

  expect(modeSwitchCount).toBe(initialCount);
});

test("lasso moves after pan", async () => {
  await goToPage(pageUrl);

  const panzoomLasso = data.features.panzoom.lasso;
  const coordinatesAsPercent = panzoomLasso["coordinates-as-percent"];

  const lassoSelection = await calcDragCoordinates(
    "layout-graph",
    coordinatesAsPercent,
  );

  await drag("layout-graph", lassoSelection.start, lassoSelection.end, true);
  await waitByID("lasso-element", { visible: true });

  const initialCount = await getCellSetCount(1);

  expect(initialCount).toBe(panzoomLasso.count);

  await clickOn("mode-pan-zoom");

  const panCoords = await calcDragCoordinates(
    "layout-graph",
    coordinatesAsPercent,
  );

  await drag("layout-graph", panCoords.start, panCoords.end, false);
  await clickOn("mode-lasso");

  const panCount = await getCellSetCount(2);

  expect(panCount).toBe(initialCount);
});

/* eslint-enable no-await-in-loop -- await in loop is needed to emulate sequential user actions */

/*
Tests included below are specific to annotation features
*/

// eslint-disable-next-line @typescript-eslint/no-explicit-any --- FIXME: disabled temporarily on migrate to TS.
async function setup(config: any) {
  await goToPage(pageUrl);
  // page
  //  .on('console', message =>
  //    console.log(`${message.type().substr(0, 3).toUpperCase()} ${message.text()}`))
  if (config.withSubset) {
    await subset({ x1: 0.1, y1: 0.1, x2: 0.8, y2: 0.8 });
  }
}

describe.each([
  { withSubset: true, tag: "subset" },
  { withSubset: false, tag: "whole" },
])("geneSET crud operations and interactions", (config) => {
  test("brush on geneset mean", async () => {
    await setup(config);
    await createGeneset(meanExpressionBrushGenesetName);
    await addGeneToSetAndExpand(meanExpressionBrushGenesetName, "SIK1");

    const histBrushableAreaId = `histogram-${meanExpressionBrushGenesetName}-plot-brushable-area`;

    const coords = await calcDragCoordinates(histBrushableAreaId, {
      x1: 0.1,
      y1: 0.5,
      x2: 0.9,
      y2: 0.5,
    });
    await drag(histBrushableAreaId, coords.start, coords.end);
    await clickOn(`cellset-button-1`);
    const cellCount = await getCellSetCount(1);
    if (config.withSubset) {
      expect(cellCount).toBe("113");
    } else {
      expect(cellCount).toBe("131");
    }
  });
  test("color by mean expression", async () => {
    await setup(config);
    await createGeneset(meanExpressionBrushGenesetName);
    await addGeneToSetAndExpand(meanExpressionBrushGenesetName, "SIK1");

    await colorByGeneset(meanExpressionBrushGenesetName);
    await assertColorLegendLabel(meanExpressionBrushGenesetName);
  });
  test("diffexp", async () => {
    if (config.withSubset) return;

    await setup(config);

    // set the two cell sets to b cells vs nk cells
    await expandCategory(`louvain`);
    await clickOn(`louvain:category-select`);
    await clickOn(`categorical-value-select-louvain-B cells`);
    await clickOn(`cellset-button-1`);
    await clickOn(`categorical-value-select-louvain-B cells`);
    await clickOn(`categorical-value-select-louvain-NK cells`);
    await clickOn(`cellset-button-2`);

    // run diffexp
    await clickOn(`diffexp-button`);
    await waitByClass("pop-1-geneset-expand");
    await expect(page).toClick(getTestClass("pop-1-geneset-expand"));

    await waitUntilNoSkeletonDetected();

    let genesHTML = await getOneElementInnerHTML(
      getTestClass("gene-set-genes"),
    );

    expect(genesHTML).toMatchSnapshot();

    // (thuang): We need to assert Pop2 geneset is expanded, because sometimes
    // the click is so fast that it's not registered
    await tryUntil(async () => {
      await expect(page).toClick(getTestClass("pop-1-geneset-expand"));
      await expect(page).toClick(getTestClass("pop-2-geneset-expand"));

      await waitUntilNoSkeletonDetected();

      const geneset = await page.$(getTestClass("geneset"));
      expect(geneset).toBeTruthy();

      await waitByClass("geneset");
      // (thuang): Assumes Pop2 geneset has NKG7 gene
      await waitByID("NKG7:gene-label");
    });

    genesHTML = await getOneElementInnerHTML(getTestClass("gene-set-genes"));

    expect(genesHTML).toMatchSnapshot();

    async function waitUntilNoSkeletonDetected() {
      await tryUntil(async () => {
        const skeleton = await page.$(`.${BLUEPRINT_SKELETON_CLASS_NAME}`);
        expect(skeleton).toBeFalsy();
      });
    }
  });
  test("create a new geneset and undo/redo", async () => {
    if (config.withSubset) return;

    await setup(config);

    const genesetName = `test-geneset-foo-123`;
    await assertGenesetDoesNotExist(genesetName);
    await createGeneset(genesetName);
    /* note: as of June 2021, the aria label is in the truncate component which clones the element */
    await assertGenesetExists(genesetName);
    await keyboardUndo();
    await assertGenesetDoesNotExist(genesetName);
    await keyboardRedo();
    await assertGenesetExists(genesetName);
  });
  test("edit geneset name and undo/redo", async () => {
    await setup(config);
    await createGeneset(editableGenesetName);
    await editGenesetName(editableGenesetName, editText);
    await assertGenesetExists(newGenesetName);
    await keyboardUndo();
    await assertGenesetExists(editableGenesetName);
    await keyboardRedo();
    await assertGenesetExists(newGenesetName);
  });
  test("delete a geneset and undo/redo", async () => {
    if (config.withSubset) return;

    await setup(config);
    await createGeneset(genesetToDeleteName);
    await deleteGeneset(genesetToDeleteName);
    await keyboardUndo();
    await assertGenesetExists(genesetToDeleteName);
    await keyboardRedo();
    await assertGenesetDoesNotExist(genesetToDeleteName);
  });
  test("geneset description", async () => {
    if (config.withSubset) return;

    await setup(config);
    await createGeneset(genesetToCheckForDescription);
    await clickOnUntil(
      `${genesetToCheckForDescription}:geneset-expand`,
      async () => {
        expect(page).toMatchElement(getTestId(genesetDescriptionID), {
          text: genesetDescriptionString,
        });
      },
    );
  });
});

describe.each([
  { withSubset: true, tag: "subset" },
  { withSubset: false, tag: "whole" },
])("GENE crud operations and interactions", (config) => {
  test("add a gene to geneset and undo/redo", async () => {
    await setup(config);
    await createGeneset(setToAddGeneTo);
    await addGeneToSetAndExpand(setToAddGeneTo, geneToAddToSet);
    await assertGeneExistsInGeneset(geneToAddToSet);
    await keyboardUndo();
    await assertGeneDoesNotExist(geneToAddToSet);
    await keyboardRedo();
    await assertGeneExistsInGeneset(geneToAddToSet);
  });
  test("expand gene and brush", async () => {
    await setup(config);
    await createGeneset(brushThisGeneGeneset);
    await addGeneToSetAndExpand(brushThisGeneGeneset, geneToBrushAndColorBy);
    await expandGene(geneToBrushAndColorBy);
    const histBrushableAreaId = `histogram-${geneToBrushAndColorBy}-plot-brushable-area`;

    const coords = await calcDragCoordinates(histBrushableAreaId, {
      x1: 0.25,
      y1: 0.5,
      x2: 0.55,
      y2: 0.5,
    });
    await drag(histBrushableAreaId, coords.start, coords.end);
    const cellCount = await getCellSetCount(1);
    if (config.withSubset) {
      expect(cellCount).toBe(subsetGeneBrushedCellCount);
    } else {
      expect(cellCount).toBe(geneBrushedCellCount);
    }
  });
  test("color by gene in geneset", async () => {
    await setup(config);
    await createGeneset(meanExpressionBrushGenesetName);
    await addGeneToSetAndExpand(meanExpressionBrushGenesetName, "SIK1");

    await colorByGene("SIK1");
    await assertColorLegendLabel("SIK1");
  });
  test("delete gene from geneset and undo/redo", async () => {
    // We've already deleted the gene
    if (config.withSubset) return;

    await setup(config);
    await createGeneset(setToRemoveFrom);
    await addGeneToSetAndExpand(setToRemoveFrom, geneToRemove);

    await removeGene(geneToRemove);
    await assertGeneDoesNotExist(geneToRemove);
    await keyboardUndo();
    await assertGeneExistsInGeneset(geneToRemove);
    await keyboardRedo();
    await assertGeneDoesNotExist(geneToRemove);
  });
  test("open gene info card and hide/remove", async () => {
    await setup(config);
    await addGeneToSearch(geneToRequestInfo);
    await requestGeneInfo(geneToRequestInfo);
    await assertGeneInfoCardExists(geneToRequestInfo);
    await minimizeGeneInfo();
    await assertGeneInfoCardIsMinimized(geneToRequestInfo);
    await removeGeneInfo();
    await assertGeneInfoDoesNotExist(geneToRequestInfo);
  });
});
