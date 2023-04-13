import { expect, test } from "@playwright/test";
import { ADD_GENE_BTN } from "tests/utils/constants";
import {
  goToWMG,
  searchAndAddGene,
  searchAndAddTissue,
} from "tests/utils/geneAndTissueUtils";
const { describe } = test;

describe("Rankit value tests", () => {
  test("Should verify MALAT1 is expressed in all values", async ({
    page,
  }) => {
    const GENE = "MALAT1";
    const TISSUE = "blood";
    await goToWMG(page);
    // click +Tissue button
    await searchAndAddTissue(page, TISSUE);

    // add gene
    await searchAndAddGene(page, GENE);
    
  });

});
