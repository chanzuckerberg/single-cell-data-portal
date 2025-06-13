/**
 * Test suite for filter-related queries.
 */

// App dependencies
import { expect } from "@playwright/test";
import {
  calculateMonthsSincePublication,
  calculateRecency,
  CollectionResponse,
  createPublicationDateValues,
  createTaggedTissueOntology,
  TissueOntology,
  processMultiValueOntologies,
} from "src/common/queries/filter";
import { PUBLICATION_DATE_VALUES } from "src/components/common/Filter/common/constants";
import { TISSUE_TYPE } from "src/components/common/Filter/common/entities";
import { test } from "tests/common/test";

const { describe } = test;

const ONTOLOGY_TERM_ID_ASIAN = "HANCESTRO:0008";
const ONTOLOGY_TERM_ID_DUTCH = "HANCESTRO:0320";
const ONTOLOGY_TERM_ID_CUBAN = "HANCESTRO:0405";
const ONTOLOGY_TERM_ID_BRAIN = "UBERON:0000955";

const LABEL_ASIAN = "Asian";
const LABEL_DUTCH = "Dutch";
const LABEL_CUBAN = "Cuban";

describe("filter", () => {
  describe("Calculate Months Since Publication", () => {
    test("calculates one month since publication", () => {
      // Publication date - 12/2021
      // Today's date - 1/2022
      const monthsSincePublication = calculateMonthsSincePublication(
        1,
        2022,
        12,
        2021
      );
      expect(monthsSincePublication).toEqual(1);
    });
    test("calculates six months since publication", () => {
      // Publication date - 12/2021
      // Today's date - 1/2022
      const monthsSincePublication = calculateMonthsSincePublication(
        1,
        2022,
        7,
        2021
      );
      expect(monthsSincePublication).toEqual(6);
    });
    test("calculates thirteen months since publication", () => {
      // Publication date - 12/2021
      // Today's date - 1/2022
      const monthsSincePublication = calculateMonthsSincePublication(
        1,
        2022,
        12,
        2020
      );
      expect(monthsSincePublication).toEqual(13);
    });
  });
  test("calculates 23 months since publication", () => {
    // Publication date - 12/2021
    // Today's date - 1/2022
    const monthsSincePublication = calculateMonthsSincePublication(
      12,
      2022,
      1,
      2021
    );
    expect(monthsSincePublication).toEqual(23);
  });
  describe("Calculate Date Bins", () => {
    test("calculates bins for 1 month since publication", () => {
      const dateBins = createPublicationDateValues(1);
      // Expecting all date ranges (that is, 1, 3, 6, 12, 24 and 36).
      expect(dateBins.length).toEqual(6);
      validateDateValue(dateBins, PUBLICATION_DATE_VALUES);
    });
    test("calculates bins for 2 months since publication", () => {
      const dateBins = createPublicationDateValues(2);
      // Expecting 3, 6, 12, 24 and 36.
      expect(dateBins.length).toEqual(5);
      validateDateValue(dateBins, PUBLICATION_DATE_VALUES.slice(1));
    });
    test("calculates bins for 4 months since publication", () => {
      const dateBins = createPublicationDateValues(4);
      // Expecting 6, 12, 24 and 36.
      expect(dateBins.length).toEqual(4);
      validateDateValue(dateBins, PUBLICATION_DATE_VALUES.slice(2));
    });
    test("calculates bins for 7 months since publication", () => {
      const dateBins = createPublicationDateValues(7);
      // Expecting 12, 24 and 36.
      expect(dateBins.length).toEqual(3);
      validateDateValue(dateBins, PUBLICATION_DATE_VALUES.slice(3));
    });
    test("calculates bins for 32 months since publication", () => {
      const dateBins = createPublicationDateValues(32);
      // Expecting 36.
      expect(dateBins.length).toEqual(1);
      validateDateValue(dateBins, PUBLICATION_DATE_VALUES.slice(5));
    });
  });
  describe("Calculate Recency", () => {
    test("calculates recency for collection with publisher metadata", () => {
      const publishedAt = 1646092800;
      const collection = {
        publisher_metadata: {
          published_at: publishedAt,
          published_day: 11,
          published_month: 1,
          published_year: 2022,
        },
      } as CollectionResponse;
      const recency = calculateRecency(
        collection,
        collection.publisher_metadata
      );
      expect(recency).toEqual(publishedAt);
    });
    test("calculates recency for collection with revised at", () => {
      const revisedAt = 1644527777.095609; // JS month
      const collection = {
        // No publisher metadata
        revised_at: revisedAt,
      } as CollectionResponse;
      const recency = calculateRecency(
        collection,
        collection.publisher_metadata
      );
      expect(recency).toEqual(revisedAt);
    });
    test("calculates recency for collection with published at", () => {
      const publishedAt = 1644526776.095609;
      const collection = {
        // No publisher metadata or revised_at
        published_at: publishedAt,
      } as CollectionResponse;
      const recency = calculateRecency(
        collection,
        collection.publisher_metadata
      );
      expect(recency).toEqual(publishedAt);
    });
  });
  describe("Process tissue type", () => {
    test("handles 3.x.x format", () => {
      // TODO remove this test with #6266.
      const tissue = {
        label: "brain",
        ontology_term_id: ONTOLOGY_TERM_ID_BRAIN,
      } as TissueOntology; // Force 3.x.x format.

      const processedTissue = createTaggedTissueOntology(tissue);
      expect(processedTissue.label).toEqual(tissue.label);
      expect(processedTissue.ontology_term_id).toEqual(tissue.ontology_term_id);
    });
    test("handles organoid", () => {
      const tissue = {
        label: "brain",
        ontology_term_id: ONTOLOGY_TERM_ID_BRAIN,
        tissue_type: TISSUE_TYPE.ORGANOID,
      };

      const processedTissue = createTaggedTissueOntology(tissue);
      expect(processedTissue.label).toEqual(
        `${tissue.label} (${TISSUE_TYPE.ORGANOID})`
      );
      expect(processedTissue.ontology_term_id).toEqual(
        `${tissue.ontology_term_id} (${TISSUE_TYPE.ORGANOID})`
      );
    });
    test("handles cell culture", () => {
      const tissue = {
        label: "epithelial cell",
        ontology_term_id: "CL:0000066",
        tissue_type: TISSUE_TYPE.CELL_CULTURE,
      };

      const processedTissue = createTaggedTissueOntology(tissue);
      expect(processedTissue.label).toEqual(
        `${tissue.label} (${TISSUE_TYPE.CELL_CULTURE})`
      );
      expect(processedTissue.ontology_term_id).toEqual(
        `${tissue.ontology_term_id} (${TISSUE_TYPE.CELL_CULTURE})`
      );
    });
    test("handles tissue", () => {
      const tissue = {
        label: "brain",
        ontology_term_id: ONTOLOGY_TERM_ID_BRAIN,
        tissue_type: TISSUE_TYPE.TISSUE,
      };

      const processedTissue = createTaggedTissueOntology(tissue);
      expect(processedTissue.label).toEqual(tissue.label);
      expect(processedTissue.ontology_term_id).toEqual(tissue.ontology_term_id);
    });
  });
  describe("Process Self Reported Ethnicity", () => {
    test("splits multiethnicity", () => {
      const selfReportedEthnicity = {
        label: `${LABEL_ASIAN}||${LABEL_DUTCH}||${LABEL_CUBAN}`,
        ontology_term_id: `${ONTOLOGY_TERM_ID_ASIAN}||${ONTOLOGY_TERM_ID_DUTCH}||${ONTOLOGY_TERM_ID_CUBAN}`,
      };

      const ethnicities = processMultiValueOntologies([selfReportedEthnicity]);

      expect(ethnicities.length).toEqual(3);
      const [asian, dutch, cuban] = ethnicities;
      expect(asian.label).toEqual(LABEL_ASIAN);
      expect(asian.ontology_term_id).toEqual(ONTOLOGY_TERM_ID_ASIAN);
      expect(dutch.label).toEqual(LABEL_DUTCH);
      expect(dutch.ontology_term_id).toEqual(ONTOLOGY_TERM_ID_DUTCH);
      expect(cuban.label).toEqual(LABEL_CUBAN);
      expect(cuban.ontology_term_id).toEqual(ONTOLOGY_TERM_ID_CUBAN);
    });
    test("handles single ethnicity", () => {
      const selfReportedEthnicity = {
        label: `${LABEL_ASIAN}`,
        ontology_term_id: `${ONTOLOGY_TERM_ID_ASIAN}`,
      };
      const ethnicities = processMultiValueOntologies([selfReportedEthnicity]);
      expect(ethnicities.length).toEqual(1);
      const [asian] = ethnicities;
      expect(asian.label).toEqual(LABEL_ASIAN);
      expect(asian.ontology_term_id).toEqual(ONTOLOGY_TERM_ID_ASIAN);
    });
  });
  /**
   * (thuang): Disease examples:
   * [
   *  {'label': "atrial fibrillation, familial, 16 || mitral valve insufficiency",'ontology_term_id': "MONDO:0800349 || MONDO:1030008"},
   *  {'label': "Hodgkin's lymphoma, lymphocytic-histiocytic predominance",'ontology_term_id': "MONDO:0004604"},
   *  {'label': "Hodgkin's lymphoma, lymphocytic-histiocytic predominance || Weil's disease",'ontology_term_id': "MONDO:0004604 || MONDO:0043004"},
   *  {'label': "mitral valve insufficiency",'ontology_term_id': "MONDO:1030008"},
   *  {'label': "normal",'ontology_term_id': "PATO:0000461"},
   * ]
   * @see: https://app.zenhub.com/workspaces/single-cell-5e2a191dad828d52cc78b028/issues/gh/chanzuckerberg/single-cell/748
   */
  describe("Process Disease", () => {
    test("splits multidisease", () => {
      const disease = {
        label:
          "atrial fibrillation, familial, 16 || mitral valve insufficiency || Hodgkin's lymphoma, lymphocytic-histiocytic predominance",
        ontology_term_id: "MONDO:0800349 || MONDO:1030008 || MONDO:0004604",
      };
      const diseases = processMultiValueOntologies([disease]);
      expect(diseases.length).toEqual(3);
      const [disease1, disease2, disease3] = diseases;
      expect(disease1.label).toEqual("atrial fibrillation, familial, 16");
      expect(disease1.ontology_term_id).toEqual("MONDO:0800349");
      expect(disease2.label).toEqual("mitral valve insufficiency");
      expect(disease2.ontology_term_id).toEqual("MONDO:1030008");
      expect(disease3.label).toEqual(
        "Hodgkin's lymphoma, lymphocytic-histiocytic predominance"
      );
      expect(disease3.ontology_term_id).toEqual("MONDO:0004604");
    });
  });
});

/**
 * Check each returned, actual date value and confirm it is in the expected.
 * @param actualDataValues - Array of numbers (in the result) to check.
 * @param expectedDateValues - Array of numbers that are expected in the result.
 */
function validateDateValue(
  actualDataValues: number[],
  expectedDateValues: number[]
) {
  actualDataValues.forEach((dateValue: number, index: number) => {
    expect(dateValue).toEqual(expectedDateValues[index]);
  });
}
