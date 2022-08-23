/**
 * Test suite for category filter hook.
 */

// App dependencies
import { buildNextOntologyCategoryFilters } from "src/common/hooks/useCategoryFilter";
import { CATEGORY_FILTER_CONFIGS_BY_ID } from "src/components/common/Filter/common/constants";
import {
  CATEGORY_FILTER_ID,
  CuratedOntologyCategoryFilterConfig,
} from "src/components/common/Filter/common/entities";

describe("useCategoryFilter", () => {
  describe("Ontology Category", () => {
    const ONTOLOGY_ID_HUMAN_PRENATAL = "HsapDv:0000045";
    const ONTOLOGY_ID_HUMAN_EMBRYONIC_HUMAN = "HsapDv:0000002";
    const ONTOLOGY_ID_HUMAN_CARNEGIE_CS1 = "HsapDv:0000003";
    const ONTOLOGY_ID_HUMAN_NEURULA_CS7_8 = "HsapDv:0000012";
    const ONTOLOGY_ID_HUMAN_FETAL = "HsapDv:0000037";
    const ONTOLOGY_ID_HUMAN_IMMATURE = "HsapDv:0000080";
    const ONTOLOGY_ID_HUMAN_NEWBORN = "HsapDv:0000082";
    const ONTOLOGY_ID_HUMAN_INFANT = "HsapDv:0000083";
    const ONTOLOGY_ID_HUMAN_CHILD = "HsapDv:0000081";

    // Set of category value keys returned from the BE, resulting in the following available tree structure for
    // human:
    //
    // Prenatal human - HsapDv:0000045
    //   Embryonic - HsapDv:0000002
    //     Carnegie CS1 - HsapDv:0000003
    //     Neurula - HsapDv:0000012
    //   Fetal - HsapDv:0000037
    // Immature - HsapDv:0000080
    //   Newborn human - HsapDv:0000082
    //   Infant - HsapDv:0000083
    //   Child - HsapDv:0000081
    const CATEGORY_VALUE_KEYS = new Set([
      ONTOLOGY_ID_HUMAN_PRENATAL,
      ONTOLOGY_ID_HUMAN_EMBRYONIC_HUMAN,
      ONTOLOGY_ID_HUMAN_CARNEGIE_CS1,
      ONTOLOGY_ID_HUMAN_NEURULA_CS7_8,
      ONTOLOGY_ID_HUMAN_FETAL,
      ONTOLOGY_ID_HUMAN_IMMATURE,
      ONTOLOGY_ID_HUMAN_CHILD,
      ONTOLOGY_ID_HUMAN_NEWBORN,
      ONTOLOGY_ID_HUMAN_INFANT,
    ]);

    describe("buildNextOntologyCategoryFilters", () => {
      /**
       * Remove selected leaf node, no siblings.
       *
       * - Current selected set: Carnegie CS1, Neurula
       *
       * - Selected value: Neurula
       */
      it(`de-selects leaf ${ONTOLOGY_ID_HUMAN_CARNEGIE_CS1}, no siblings`, () => {
        const idToDeselect = ONTOLOGY_ID_HUMAN_CARNEGIE_CS1;
        const filters = [
          {
            id: CATEGORY_FILTER_ID.DEVELOPMENT_STAGE,
            value: [idToDeselect], // Add ID to de-select as already selected
          },
        ];

        const nextFilters = buildNextOntologyCategoryFilters(
          CATEGORY_FILTER_ID.DEVELOPMENT_STAGE,
          [idToDeselect],
          filters,
          CATEGORY_VALUE_KEYS,
          (
            CATEGORY_FILTER_CONFIGS_BY_ID[
              CATEGORY_FILTER_ID.DEVELOPMENT_STAGE
            ] as CuratedOntologyCategoryFilterConfig
          ).mask
        );
        expect(nextFilters.length).toEqual(0);
      });

      /**
       * Remove selected leaf node with siblings.
       *
       * - Current selected set: Carnegie CS1, Neurula
       *
       * - Selected value: Neurula
       */
      it(`de-selects leaf ${ONTOLOGY_ID_HUMAN_CARNEGIE_CS1}, with siblings`, () => {
        const idToDeselect = ONTOLOGY_ID_HUMAN_CARNEGIE_CS1;
        const filters = [
          {
            id: CATEGORY_FILTER_ID.DEVELOPMENT_STAGE,
            value: [
              ONTOLOGY_ID_HUMAN_EMBRYONIC_HUMAN,
              idToDeselect,
              ONTOLOGY_ID_HUMAN_NEURULA_CS7_8,
            ], // Add ID to de-select as already selected
          },
        ];
        const nextFilters = buildNextOntologyCategoryFilters(
          CATEGORY_FILTER_ID.DEVELOPMENT_STAGE,
          [idToDeselect],
          filters,
          CATEGORY_VALUE_KEYS,
          (
            CATEGORY_FILTER_CONFIGS_BY_ID[
              CATEGORY_FILTER_ID.DEVELOPMENT_STAGE
            ] as CuratedOntologyCategoryFilterConfig
          ).mask
        );
        expect(nextFilters.length).toEqual(1);
        expect(
          nextFilters.includes(ONTOLOGY_ID_HUMAN_NEURULA_CS7_8)
        ).toBeTruthy();
      });

      /**
       * Remove intermediate selected parent node, siblings unselected.
       *
       * - Selected value: Embryonic human
       *
       * - Current selected set: Embryonic human, Carnegie CS1, Neurula
       *
       * - Expected selected set: []
       */
      it(`de-selects node ${ONTOLOGY_ID_HUMAN_EMBRYONIC_HUMAN}, sibling unselected`, () => {
        const idToDeselect = ONTOLOGY_ID_HUMAN_EMBRYONIC_HUMAN;
        const filters = [
          {
            id: CATEGORY_FILTER_ID.DEVELOPMENT_STAGE,
            value: [
              idToDeselect,
              ONTOLOGY_ID_HUMAN_CARNEGIE_CS1,
              ONTOLOGY_ID_HUMAN_NEURULA_CS7_8,
            ], // Add ID to de-select as already selected
          },
        ];
        const nextFilters = buildNextOntologyCategoryFilters(
          CATEGORY_FILTER_ID.DEVELOPMENT_STAGE,
          [idToDeselect],
          filters,
          CATEGORY_VALUE_KEYS,
          (
            CATEGORY_FILTER_CONFIGS_BY_ID[
              CATEGORY_FILTER_ID.DEVELOPMENT_STAGE
            ] as CuratedOntologyCategoryFilterConfig
          ).mask
        );
        expect(nextFilters.length).toEqual(0);
      });

      /**
       * Remove intermediate selected parent node with siblings selected.
       *
       * - Selected value: Embryonic human
       *
       * - Current selected set: Prental, Embryonic human, Carenegie CS1, Neurula, Fetal
       *
       * - Expected selected set: [Fetal]
       */
      it(`de-selects node ${ONTOLOGY_ID_HUMAN_EMBRYONIC_HUMAN}, sibling selected`, () => {
        const idToDeselect = ONTOLOGY_ID_HUMAN_EMBRYONIC_HUMAN;
        const filters = [
          {
            id: CATEGORY_FILTER_ID.DEVELOPMENT_STAGE,
            value: [
              ONTOLOGY_ID_HUMAN_PRENATAL,
              idToDeselect,
              ONTOLOGY_ID_HUMAN_CARNEGIE_CS1,
              ONTOLOGY_ID_HUMAN_NEURULA_CS7_8,
              ONTOLOGY_ID_HUMAN_FETAL,
            ], // Add ID to de-select as already selected
          },
        ];
        const nextFilters = buildNextOntologyCategoryFilters(
          CATEGORY_FILTER_ID.DEVELOPMENT_STAGE,
          [idToDeselect],
          filters,
          CATEGORY_VALUE_KEYS,
          (
            CATEGORY_FILTER_CONFIGS_BY_ID[
              CATEGORY_FILTER_ID.DEVELOPMENT_STAGE
            ] as CuratedOntologyCategoryFilterConfig
          ).mask
        );
        expect(nextFilters.length).toEqual(1);
        expect(nextFilters.includes(ONTOLOGY_ID_HUMAN_FETAL)).toBeTruthy();
      });

      /**
       * Select leaf, nothing currently selected => selects self.
       *
       * - Selected value: Carnegie (CS1) HsapDv:0000003.
       *
       * - Expected selected set: Carnegie (CS1) HsapDv:0000003
       */
      it(`selects leaf ${ONTOLOGY_ID_HUMAN_CARNEGIE_CS1}, no siblings currently selected`, () => {
        const idToSelected = ONTOLOGY_ID_HUMAN_CARNEGIE_CS1;
        const nextFilters = buildNextOntologyCategoryFilters(
          CATEGORY_FILTER_ID.DEVELOPMENT_STAGE,
          [idToSelected],
          [], // No filters selected
          CATEGORY_VALUE_KEYS,
          (
            CATEGORY_FILTER_CONFIGS_BY_ID[
              CATEGORY_FILTER_ID.DEVELOPMENT_STAGE
            ] as CuratedOntologyCategoryFilterConfig
          ).mask
        );
        expect(nextFilters.length).toEqual(1);
        expect(nextFilters[0]).toEqual(idToSelected);
      });

      /**
       * Select leaf, some siblings selected => leaf selected.
       *
       * - Current selected values: Newborn
       *
       * - Selected value: Infant
       *
       * - Expected selected set: Newborn, Infant
       */
      it(`selects leaf ${ONTOLOGY_ID_HUMAN_INFANT}, some siblings currently selected`, () => {
        const nextFilters = buildNextOntologyCategoryFilters(
          CATEGORY_FILTER_ID.DEVELOPMENT_STAGE,
          [ONTOLOGY_ID_HUMAN_INFANT],
          [
            {
              id: CATEGORY_FILTER_ID.DEVELOPMENT_STAGE,
              value: [ONTOLOGY_ID_HUMAN_NEWBORN],
            },
          ],
          CATEGORY_VALUE_KEYS,
          (
            CATEGORY_FILTER_CONFIGS_BY_ID[
              CATEGORY_FILTER_ID.DEVELOPMENT_STAGE
            ] as CuratedOntologyCategoryFilterConfig
          ).mask
        );
        expect(nextFilters.length).toEqual(2);
        expect(nextFilters.includes(ONTOLOGY_ID_HUMAN_NEWBORN)).toBeTruthy();
        expect(nextFilters.includes(ONTOLOGY_ID_HUMAN_INFANT)).toBeTruthy();
      });

      /**
       * Select leaf, all siblings selected => siblings selected, parent selected
       *
       * - Current selected values: Newborn, infant
       *
       * - Selected value: Child
       *
       * - Expected selected set: Immature, Newborn, Infant, Child
       */
      it(`selects leaf ${ONTOLOGY_ID_HUMAN_CHILD}, all siblings currently selected`, () => {
        const nextFilters = buildNextOntologyCategoryFilters(
          CATEGORY_FILTER_ID.DEVELOPMENT_STAGE,
          [ONTOLOGY_ID_HUMAN_CHILD],
          [
            {
              id: CATEGORY_FILTER_ID.DEVELOPMENT_STAGE,
              value: [ONTOLOGY_ID_HUMAN_NEWBORN, ONTOLOGY_ID_HUMAN_INFANT],
            },
          ],
          CATEGORY_VALUE_KEYS,
          (
            CATEGORY_FILTER_CONFIGS_BY_ID[
              CATEGORY_FILTER_ID.DEVELOPMENT_STAGE
            ] as CuratedOntologyCategoryFilterConfig
          ).mask
        );
        expect(nextFilters.length).toEqual(4);
        expect(nextFilters.includes(ONTOLOGY_ID_HUMAN_NEWBORN)).toBeTruthy();
        expect(nextFilters.includes(ONTOLOGY_ID_HUMAN_INFANT)).toBeTruthy();
        expect(nextFilters.includes(ONTOLOGY_ID_HUMAN_CHILD)).toBeTruthy();
        expect(nextFilters.includes(ONTOLOGY_ID_HUMAN_IMMATURE)).toBeTruthy();
      });

      /**
       * Select leaf, all siblings and aunt selected => ancestor selected.
       *
       * - Current selected values: Neurula, Fetal
       *
       * - Selected value: Carnegie (CS1)
       *
       * - Expected selected set: Prenatal, Embryonic, Neurula, Carnegie and Fetal
       */
      it(`selects leaf ${ONTOLOGY_ID_HUMAN_CARNEGIE_CS1}, all siblings and aunt currently selected`, () => {
        const nextFilters = buildNextOntologyCategoryFilters(
          CATEGORY_FILTER_ID.DEVELOPMENT_STAGE,
          [ONTOLOGY_ID_HUMAN_CARNEGIE_CS1],
          [
            {
              id: CATEGORY_FILTER_ID.DEVELOPMENT_STAGE,
              value: [ONTOLOGY_ID_HUMAN_NEURULA_CS7_8, ONTOLOGY_ID_HUMAN_FETAL],
            },
          ],
          CATEGORY_VALUE_KEYS,
          (
            CATEGORY_FILTER_CONFIGS_BY_ID[
              CATEGORY_FILTER_ID.DEVELOPMENT_STAGE
            ] as CuratedOntologyCategoryFilterConfig
          ).mask
        );
        expect(nextFilters.length).toEqual(5);
        expect(nextFilters.includes(ONTOLOGY_ID_HUMAN_PRENATAL)).toBeTruthy();
        expect(
          nextFilters.includes(ONTOLOGY_ID_HUMAN_EMBRYONIC_HUMAN)
        ).toBeTruthy();
        expect(
          nextFilters.includes(ONTOLOGY_ID_HUMAN_CARNEGIE_CS1)
        ).toBeTruthy();
        expect(
          nextFilters.includes(ONTOLOGY_ID_HUMAN_NEURULA_CS7_8)
        ).toBeTruthy();
        expect(nextFilters.includes(ONTOLOGY_ID_HUMAN_FETAL)).toBeTruthy();
      });

      /**
       * Select parent node, nothing currently selected => selects self and children.
       *
       * - Current selected values: -
       *
       * - Selected value: Embryonic human
       *
       * - Expected selected set: Embryonic human, Carnegie, Neurula
       */
      it(`selects node ${ONTOLOGY_ID_HUMAN_EMBRYONIC_HUMAN}, nothing currently selected`, () => {
        const nextFilters = buildNextOntologyCategoryFilters(
          CATEGORY_FILTER_ID.DEVELOPMENT_STAGE,
          [ONTOLOGY_ID_HUMAN_EMBRYONIC_HUMAN],
          [], // No filters
          CATEGORY_VALUE_KEYS,
          (
            CATEGORY_FILTER_CONFIGS_BY_ID[
              CATEGORY_FILTER_ID.DEVELOPMENT_STAGE
            ] as CuratedOntologyCategoryFilterConfig
          ).mask
        );
        expect(nextFilters.length).toEqual(3);
        expect(
          nextFilters.includes(ONTOLOGY_ID_HUMAN_EMBRYONIC_HUMAN)
        ).toBeTruthy();
        expect(
          nextFilters.includes(ONTOLOGY_ID_HUMAN_CARNEGIE_CS1)
        ).toBeTruthy();
        expect(
          nextFilters.includes(ONTOLOGY_ID_HUMAN_NEURULA_CS7_8)
        ).toBeTruthy();
      });

      /**
       * Select node, child currently selected => selects self and all children.
       *
       * - Current selected values: Carnegie (CS1)
       *
       * - Selected value: Embryonic human
       *
       * - Expected selected set: Embryonic human, Carnegie, Neurula
       */
      it(`selects node ${ONTOLOGY_ID_HUMAN_EMBRYONIC_HUMAN}, child currently selected`, () => {
        const nextFilters = buildNextOntologyCategoryFilters(
          CATEGORY_FILTER_ID.DEVELOPMENT_STAGE,
          [ONTOLOGY_ID_HUMAN_EMBRYONIC_HUMAN],
          [
            {
              id: CATEGORY_FILTER_ID.DEVELOPMENT_STAGE,
              value: [ONTOLOGY_ID_HUMAN_CARNEGIE_CS1],
            },
          ],
          CATEGORY_VALUE_KEYS,
          (
            CATEGORY_FILTER_CONFIGS_BY_ID[
              CATEGORY_FILTER_ID.DEVELOPMENT_STAGE
            ] as CuratedOntologyCategoryFilterConfig
          ).mask
        );
        expect(nextFilters.length).toEqual(3);
        expect(
          nextFilters.includes(ONTOLOGY_ID_HUMAN_EMBRYONIC_HUMAN)
        ).toBeTruthy();
        expect(
          nextFilters.includes(ONTOLOGY_ID_HUMAN_CARNEGIE_CS1)
        ).toBeTruthy();
        expect(
          nextFilters.includes(ONTOLOGY_ID_HUMAN_NEURULA_CS7_8)
        ).toBeTruthy();
      });

      /**
       * Select node, sibling currently selected => selects parent and all descendents.
       *
       * - Current selected values: Fetal (>56 days–birth)
       *
       * - Selected value: Embryonic human
       *
       * - Expected selected set: Prenatal human, Embryonic, Carnegie, Neurula, Fetal
       */
      it(`selects node ${ONTOLOGY_ID_HUMAN_EMBRYONIC_HUMAN}, sibling currently selected`, () => {
        const nextFilters = buildNextOntologyCategoryFilters(
          CATEGORY_FILTER_ID.DEVELOPMENT_STAGE,
          [ONTOLOGY_ID_HUMAN_EMBRYONIC_HUMAN],
          [
            {
              id: CATEGORY_FILTER_ID.DEVELOPMENT_STAGE,
              value: [
                ONTOLOGY_ID_HUMAN_FETAL, // Sibling selected
              ],
            },
          ],
          CATEGORY_VALUE_KEYS,
          (
            CATEGORY_FILTER_CONFIGS_BY_ID[
              CATEGORY_FILTER_ID.DEVELOPMENT_STAGE
            ] as CuratedOntologyCategoryFilterConfig
          ).mask
        );
        expect(nextFilters.length).toEqual(5);
        expect(nextFilters.includes(ONTOLOGY_ID_HUMAN_PRENATAL)).toBeTruthy();
        expect(
          nextFilters.includes(ONTOLOGY_ID_HUMAN_EMBRYONIC_HUMAN)
        ).toBeTruthy();
        expect(
          nextFilters.includes(ONTOLOGY_ID_HUMAN_CARNEGIE_CS1)
        ).toBeTruthy();
        expect(
          nextFilters.includes(ONTOLOGY_ID_HUMAN_NEURULA_CS7_8)
        ).toBeTruthy();
        expect(nextFilters.includes(ONTOLOGY_ID_HUMAN_FETAL)).toBeTruthy();
      });

      /**
       * Select node, ancestor currently selected => de-selects self.
       *
       * - Current selected values: Prenatal human, Embryonic, Carnegie, Neurula, Fetal
       *
       * - Selected value: Carnegie
       *
       * - Expected selected set: Neurula, Fetal
       */
      it(`de-selects leaf ${ONTOLOGY_ID_HUMAN_CARNEGIE_CS1}, ancestor currently selected`, () => {
        const nextFilters = buildNextOntologyCategoryFilters(
          CATEGORY_FILTER_ID.DEVELOPMENT_STAGE,
          [ONTOLOGY_ID_HUMAN_CARNEGIE_CS1],
          [
            {
              id: CATEGORY_FILTER_ID.DEVELOPMENT_STAGE,
              value: [
                ONTOLOGY_ID_HUMAN_PRENATAL,
                ONTOLOGY_ID_HUMAN_EMBRYONIC_HUMAN,
                ONTOLOGY_ID_HUMAN_CARNEGIE_CS1,
                ONTOLOGY_ID_HUMAN_NEURULA_CS7_8,
                ONTOLOGY_ID_HUMAN_FETAL,
              ],
            },
          ],
          CATEGORY_VALUE_KEYS,
          (
            CATEGORY_FILTER_CONFIGS_BY_ID[
              CATEGORY_FILTER_ID.DEVELOPMENT_STAGE
            ] as CuratedOntologyCategoryFilterConfig
          ).mask
        );
        expect(nextFilters.length).toEqual(2);
        expect(
          nextFilters.includes(ONTOLOGY_ID_HUMAN_NEURULA_CS7_8)
        ).toBeTruthy();
        expect(nextFilters.includes(ONTOLOGY_ID_HUMAN_FETAL)).toBeTruthy(); // Aunt selected
      });
    });
  });
});
