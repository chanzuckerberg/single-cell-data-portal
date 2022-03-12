/**
 * Test suite for category filter hook.
 */

// App dependencies
import { buildNextOntologyCategoryFilters } from "src/common/hooks/useCategoryFilter";
import {
  CATEGORY_CONFIGS_BY_CATEGORY_KEY,
  CATEGORY_KEY,
  OntologyCategoryConfig,
} from "src/components/common/Filter/common/entities";

describe("useCategoryFilter", () => {
  describe("Ontology Category", () => {
    const ONTOLOGY_ID_HUMAN_PRENATAL = "HsapDv:0000045";
    const ONTOLOGY_ID_HUMAN_EMBRYONIC_HUMAN = "HsapDv:0000002";
    const ONTOLOGY_ID_HUMAN_CARNEGIE_CS1 = "HsapDv:0000003";
    const ONTOLOGY_ID_HUMAN_CLEAVAGE_CS2 = "HsapDv:0000004";
    const ONTOLOGY_ID_HUMAN_BLASTULA_CS3_5 = "HsapDv:0000006";
    const ONTOLOGY_ID_HUMAN_GASTRULA_CS6 = "HsapDv:0000010";
    const ONTOLOGY_ID_HUMAN_NEURULA_CS7_8 = "HsapDv:0000012";
    const ONTOLOGY_ID_HUMAN_ORGANOGENESIS_CS9_23 = "HsapDv:0000015";
    const ONTOLOGY_ID_HUMAN_FETAL = "HsapDv:0000037";

    describe("buildNextOntologyCategoryFilters", () => {
      /**
       * Remove selected leaf node
       * - Selected value - Carnegie (CS1) HsapDv:0000003.
       */
      it("de-selects ontology leaf", () => {
        const idToDeselect = ONTOLOGY_ID_HUMAN_CARNEGIE_CS1;
        const filters = [
          {
            id: CATEGORY_KEY.DEVELOPMENT_STAGE_ANCESTORS,
            value: [idToDeselect], // Add ID to de-select as already selected
          },
        ];
        const nextFilters = buildNextOntologyCategoryFilters(
          CATEGORY_KEY.DEVELOPMENT_STAGE_ANCESTORS,
          idToDeselect,
          filters,
          (
            CATEGORY_CONFIGS_BY_CATEGORY_KEY[
              CATEGORY_KEY.DEVELOPMENT_STAGE_ANCESTORS
            ] as OntologyCategoryConfig
          ).ontology
        );
        expect(nextFilters.length).toEqual(0);
      });

      /**
       * Select leaf, nothing currently selected => selects self.
       *
       * - Selected value: Carnegie (CS1) HsapDv:0000003.
       *
       * - Expected selected set: Carnegie (CS1) HsapDv:0000003
       */
      it("selects ontology leaf, no siblings currently selected", () => {
        const idToSelected = ONTOLOGY_ID_HUMAN_CARNEGIE_CS1;
        const nextFilters = buildNextOntologyCategoryFilters(
          CATEGORY_KEY.DEVELOPMENT_STAGE_ANCESTORS,
          idToSelected,
          [], // No filters selected
          (
            CATEGORY_CONFIGS_BY_CATEGORY_KEY[
              CATEGORY_KEY.DEVELOPMENT_STAGE_ANCESTORS
            ] as OntologyCategoryConfig
          ).ontology
        );
        expect(nextFilters.length).toEqual(1);
        expect(nextFilters[0]).toEqual(idToSelected);
      });

      /**
       * Select leaf, some siblings selected => leaf selected.
       *
       * - Current selected values: Blastula (CS3–5), Gastrula (CS6),
       *   Neurula (CS7–8), Organogenesis (CS9–23)
       *
       * - Selected value: Carnegie (CS1)
       *
       * - Expected selected set: Blastula (CS3–5), Gastrula (CS6),
       *   Neurula (CS7–8), Organogenesis (CS9–23) and Carnegie (CS1)
       */
      it("selects ontology leaf, all siblings currently selected", () => {
        const nextFilters = buildNextOntologyCategoryFilters(
          CATEGORY_KEY.DEVELOPMENT_STAGE_ANCESTORS,
          ONTOLOGY_ID_HUMAN_CARNEGIE_CS1,
          [
            {
              id: CATEGORY_KEY.DEVELOPMENT_STAGE_ANCESTORS,
              value: [
                ONTOLOGY_ID_HUMAN_BLASTULA_CS3_5,
                ONTOLOGY_ID_HUMAN_GASTRULA_CS6,
                ONTOLOGY_ID_HUMAN_NEURULA_CS7_8,
                ONTOLOGY_ID_HUMAN_ORGANOGENESIS_CS9_23,
              ],
            },
          ],
          (
            CATEGORY_CONFIGS_BY_CATEGORY_KEY[
              CATEGORY_KEY.DEVELOPMENT_STAGE_ANCESTORS
            ] as OntologyCategoryConfig
          ).ontology
        );
        expect(nextFilters.length).toEqual(5);
        expect(
          nextFilters.includes(ONTOLOGY_ID_HUMAN_CARNEGIE_CS1)
        ).toBeTruthy();
        expect(
          nextFilters.includes(ONTOLOGY_ID_HUMAN_BLASTULA_CS3_5)
        ).toBeTruthy();
        expect(
          nextFilters.includes(ONTOLOGY_ID_HUMAN_GASTRULA_CS6)
        ).toBeTruthy();
        expect(
          nextFilters.includes(ONTOLOGY_ID_HUMAN_NEURULA_CS7_8)
        ).toBeTruthy();
        expect(
          nextFilters.includes(ONTOLOGY_ID_HUMAN_ORGANOGENESIS_CS9_23)
        ).toBeTruthy();
      });

      /**
       * Select leaf, all siblings selected => parent selected.
       *
       * - Current selected values: Cleavage (CS2), Blastula (CS3–5), Gastrula (CS6),
       *   Neurula (CS7–8), Organogenesis (CS9–23)
       *
       * - Selected value: Carnegie (CS1)
       *
       * - Expected selected set: Embryonic human (that is, parent of all selected values)
       */
      it("selects ontology leaf, all siblings currently selected", () => {
        const nextFilters = buildNextOntologyCategoryFilters(
          CATEGORY_KEY.DEVELOPMENT_STAGE_ANCESTORS,
          ONTOLOGY_ID_HUMAN_CARNEGIE_CS1,
          [
            {
              id: CATEGORY_KEY.DEVELOPMENT_STAGE_ANCESTORS,
              value: [
                ONTOLOGY_ID_HUMAN_CLEAVAGE_CS2,
                ONTOLOGY_ID_HUMAN_BLASTULA_CS3_5,
                ONTOLOGY_ID_HUMAN_GASTRULA_CS6,
                ONTOLOGY_ID_HUMAN_NEURULA_CS7_8,
                ONTOLOGY_ID_HUMAN_ORGANOGENESIS_CS9_23,
              ],
            },
          ],
          (
            CATEGORY_CONFIGS_BY_CATEGORY_KEY[
              CATEGORY_KEY.DEVELOPMENT_STAGE_ANCESTORS
            ] as OntologyCategoryConfig
          ).ontology
        );
        expect(nextFilters.length).toEqual(1);
        expect(nextFilters[0]).toEqual(ONTOLOGY_ID_HUMAN_EMBRYONIC_HUMAN);
      });

      /**
       * Select leaf, all siblings and aunt selected => ancestor selected.
       *
       * - Current selected values: Cleavage (CS2), Blastula (CS3–5), Gastrula (CS6),
       *   Neurula (CS7–8), Organogenesis (CS9–23), Fetal (>56 days–birth)
       *
       * - Selected value: Carnegie (CS1)
       *
       * - Expected selected set: Prenatal (conception–birth) (that is, parent of all siblings is selected and aunt
       *   is selected resulting in the grandparent being selected)
       */
      it("selects ontology leaf, all siblings and aunt currently selected", () => {
        const nextFilters = buildNextOntologyCategoryFilters(
          CATEGORY_KEY.DEVELOPMENT_STAGE_ANCESTORS,
          ONTOLOGY_ID_HUMAN_CARNEGIE_CS1,
          [
            {
              id: CATEGORY_KEY.DEVELOPMENT_STAGE_ANCESTORS,
              value: [
                ONTOLOGY_ID_HUMAN_CLEAVAGE_CS2, // Sibling
                ONTOLOGY_ID_HUMAN_BLASTULA_CS3_5, // Sibling
                ONTOLOGY_ID_HUMAN_GASTRULA_CS6, // Sibling
                ONTOLOGY_ID_HUMAN_NEURULA_CS7_8, // Sibling
                ONTOLOGY_ID_HUMAN_ORGANOGENESIS_CS9_23, // Sibling
                ONTOLOGY_ID_HUMAN_FETAL, // Aunt
              ],
            },
          ],
          (
            CATEGORY_CONFIGS_BY_CATEGORY_KEY[
              CATEGORY_KEY.DEVELOPMENT_STAGE_ANCESTORS
            ] as OntologyCategoryConfig
          ).ontology
        );
        expect(nextFilters.length).toEqual(1);
        expect(nextFilters[0]).toEqual(ONTOLOGY_ID_HUMAN_PRENATAL);
      });

      /**
       * Select node, nothing currently selected => selects self.
       *
       * - Current selected values: -
       *
       * - Selected value: Embryonic human
       *
       * - Expected selected set: Embryonic human
       */
      it("selects ontology node, nothing currently selected", () => {
        const idToSelected = ONTOLOGY_ID_HUMAN_EMBRYONIC_HUMAN;
        const nextFilters = buildNextOntologyCategoryFilters(
          CATEGORY_KEY.DEVELOPMENT_STAGE_ANCESTORS,
          idToSelected,
          [], // No filters
          (
            CATEGORY_CONFIGS_BY_CATEGORY_KEY[
              CATEGORY_KEY.DEVELOPMENT_STAGE_ANCESTORS
            ] as OntologyCategoryConfig
          ).ontology
        );
        expect(nextFilters.length).toEqual(1);
        expect(nextFilters[0]).toEqual(idToSelected);
      });

      /**
       * Select node, child currently selected => removes child, selects self.
       *
       * - Current selected values: Carnegie (CS1)
       *
       * - Selected value: Embryonic human
       *
       * - Expected selected set: Embryonic human
       */
      it("selects ontology node, child currently selected", () => {
        const idToSelected = ONTOLOGY_ID_HUMAN_EMBRYONIC_HUMAN;
        const nextFilters = buildNextOntologyCategoryFilters(
          CATEGORY_KEY.DEVELOPMENT_STAGE_ANCESTORS,
          idToSelected,
          [
            {
              id: CATEGORY_KEY.DEVELOPMENT_STAGE_ANCESTORS,
              value: [
                ONTOLOGY_ID_HUMAN_CARNEGIE_CS1, // Child selected
              ],
            },
          ],
          (
            CATEGORY_CONFIGS_BY_CATEGORY_KEY[
              CATEGORY_KEY.DEVELOPMENT_STAGE_ANCESTORS
            ] as OntologyCategoryConfig
          ).ontology
        );
        expect(nextFilters.length).toEqual(1);
        expect(nextFilters[0]).toEqual(idToSelected);
      });

      /**
       * Select node, sibling currently selected => selects parent.
       *
       * - Current selected values: Fetal (>56 days–birth)
       *
       * - Selected value: Embryonic human
       *
       * - Expected selected set: Prenatal human
       */
      it("selects ontology node, sibling currently selected", () => {
        const nextFilters = buildNextOntologyCategoryFilters(
          CATEGORY_KEY.DEVELOPMENT_STAGE_ANCESTORS,
          ONTOLOGY_ID_HUMAN_EMBRYONIC_HUMAN,
          [
            {
              id: CATEGORY_KEY.DEVELOPMENT_STAGE_ANCESTORS,
              value: [
                ONTOLOGY_ID_HUMAN_FETAL, // Sibling selected
              ],
            },
          ],
          (
            CATEGORY_CONFIGS_BY_CATEGORY_KEY[
              CATEGORY_KEY.DEVELOPMENT_STAGE_ANCESTORS
            ] as OntologyCategoryConfig
          ).ontology
        );
        expect(nextFilters.length).toEqual(1);
        expect(nextFilters[0]).toEqual(ONTOLOGY_ID_HUMAN_PRENATAL); // Parent selected
      });

      /**
       * Select node, ancestor currently selected => selects self, aunts.
       *
       * - Current selected values: Prenatal human
       *
       * - Selected value: Neurula (CS7 - 8)
       *
       * - Expected selected set: Neurula (CS7 - 8), Feta (>56 days-birth)
       */
      it("selects ontology node, ancestor currently selected", () => {
        const nextFilters = buildNextOntologyCategoryFilters(
          CATEGORY_KEY.DEVELOPMENT_STAGE_ANCESTORS,
          ONTOLOGY_ID_HUMAN_NEURULA_CS7_8,
          [
            {
              id: CATEGORY_KEY.DEVELOPMENT_STAGE_ANCESTORS,
              value: [
                ONTOLOGY_ID_HUMAN_PRENATAL, // Ancestor selected
              ],
            },
          ],
          (
            CATEGORY_CONFIGS_BY_CATEGORY_KEY[
              CATEGORY_KEY.DEVELOPMENT_STAGE_ANCESTORS
            ] as OntologyCategoryConfig
          ).ontology
        );
        expect(nextFilters.length).toEqual(2);
        expect(
          nextFilters.includes(ONTOLOGY_ID_HUMAN_NEURULA_CS7_8)
        ).toBeTruthy(); // Self selected
        expect(nextFilters.includes(ONTOLOGY_ID_HUMAN_FETAL)).toBeTruthy(); // Aunt selected
      });
    });
  });
});
