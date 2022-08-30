/**
 * Test suite for category filter hook.
 */

// App dependencies
import {
  buildNextOntologyCategoryFilters,
  buildUINodesByCategoryValueId,
  MultiPanelUINode,
} from "src/common/hooks/useCategoryFilter";
import { CATEGORY_FILTER_CONFIGS_BY_ID } from "src/components/common/Filter/common/constants";
import {
  CategoryValueId,
  CATEGORY_FILTER_ID,
  CuratedOntologyCategoryFilterConfig,
  OrFilterPrefix,
} from "src/components/common/Filter/common/entities";

describe("useCategoryFilter", () => {
  describe("Multi-Panel Category", () => {
    const TERM_ID_BLADDER_LUMEN = "UBERON:0009958";
    const TERM_ID_BLADDER_ORGAN = "UBERON:0018707";
    const TERM_ID_BLOOD = "UBERON:0000178";
    const TERM_ID_BONE_MARROW = "UBERON:0002371";
    const TERM_ID_HEMATOPOIETIC_SYSTEM = "UBERON:0002390";
    const TERM_ID_KIDNEY = "UBERON:0002113";
    const TERM_ID_SPLEEN = "UBERON:0002106";
    const TERM_ID_THYMUS = "UBERON:0002370";
    const TERM_ID_RENAL_MEDULLA = "UBERON:0000362";
    const TERM_ID_RENAL_SYSTEM = "UBERON:0001008";
    const TERM_ID_UMBILICAL_CORD_BLOOD = "UBERON:0012168";
    const TERM_ID_URETER = "UBERON:0000056";
    const TERM_ID_URETHRA = "UBERON:0000057";
    const TERM_ID_VENOUS_BLOOD = "UBERON:0013756";

    const INFERRED_HEMATOPOIETIC_SYSTEM = `${OrFilterPrefix.INFERRED}:${TERM_ID_HEMATOPOIETIC_SYSTEM}`;
    const INFERRED_RENAL_SYSTEM = `${OrFilterPrefix.INFERRED}:${TERM_ID_RENAL_SYSTEM}`;
    const INFERRED_BLADDER_ORGAN = `${OrFilterPrefix.INFERRED}:${TERM_ID_BLADDER_ORGAN}`;
    const INFERRED_BLOOD = `${OrFilterPrefix.INFERRED}:${TERM_ID_BLOOD}`;
    const INFERRED_BONE_MARROW = `${OrFilterPrefix.INFERRED}:${TERM_ID_BONE_MARROW}`;
    const INFERRED_KIDNEY = `${OrFilterPrefix.INFERRED}:${TERM_ID_KIDNEY}`;
    const INFERRED_SPLEEN = `${OrFilterPrefix.INFERRED}:${TERM_ID_SPLEEN}`;
    const INFERRED_THYYMUS = `${OrFilterPrefix.INFERRED}:${TERM_ID_THYMUS}`;

    const EXPLICIT_BLADDER_LUMEN = `${OrFilterPrefix.EXPLICIT}:${TERM_ID_BLADDER_LUMEN}`;
    const EXPLICIT_BLOOD = `${OrFilterPrefix.EXPLICIT}:${TERM_ID_BLOOD}`;
    const EXPLICIT_BONE_MARROW = `${OrFilterPrefix.EXPLICIT}:${TERM_ID_BONE_MARROW}`;
    const EXPLICIT_KIDNEY = `${OrFilterPrefix.EXPLICIT}:${TERM_ID_KIDNEY}`;
    const EXPLICIT_RENAL_MEDULLA = `${OrFilterPrefix.EXPLICIT}:${TERM_ID_RENAL_MEDULLA}`;
    const EXPLICIT_SPLEEN = `${OrFilterPrefix.EXPLICIT}:${TERM_ID_SPLEEN}`;
    const EXPLICIT_THYMUS = `${OrFilterPrefix.EXPLICIT}:${TERM_ID_THYMUS}`;
    const EXPLICIT_UMBILICAL_CORD_BLOOD = `${OrFilterPrefix.EXPLICIT}:${TERM_ID_UMBILICAL_CORD_BLOOD}`;
    const EXPLICIT_URETER = `${OrFilterPrefix.EXPLICIT}:${TERM_ID_URETER}`;
    const EXPLICIT_URETHRA = `${OrFilterPrefix.EXPLICIT}:${TERM_ID_URETHRA}`;
    const EXPLICIT_VENOUS_BLOOD = `${OrFilterPrefix.EXPLICIT}:${TERM_ID_VENOUS_BLOOD}`;

    const TISSUE_SYSTEMS = [
      INFERRED_HEMATOPOIETIC_SYSTEM,
      INFERRED_RENAL_SYSTEM,
    ];
    const TISSUE_ORGANS = [
      INFERRED_BLADDER_ORGAN,
      INFERRED_BLOOD,
      INFERRED_BONE_MARROW,
      INFERRED_KIDNEY,
      INFERRED_SPLEEN,
      INFERRED_THYYMUS,
    ];
    const TISSUES = [
      EXPLICIT_BLADDER_LUMEN,
      EXPLICIT_BLOOD,
      EXPLICIT_BONE_MARROW,
      EXPLICIT_KIDNEY,
      EXPLICIT_RENAL_MEDULLA,
      EXPLICIT_SPLEEN,
      EXPLICIT_THYMUS,
      EXPLICIT_UMBILICAL_CORD_BLOOD,
      EXPLICIT_URETER,
      EXPLICIT_URETHRA,
      EXPLICIT_VENOUS_BLOOD,
    ];

    const CATEGORY_VALUE_IDS_BY_PANEL = [
      TISSUE_SYSTEMS,
      TISSUE_ORGANS,
      TISSUES,
    ];

    describe("buildUINodesByCategoryValueId", () => {
      let uiNodesByCategoryValueId: Map<CategoryValueId, MultiPanelUINode>;
      beforeAll(() => {
        uiNodesByCategoryValueId = buildUINodesByCategoryValueId(
          CATEGORY_VALUE_IDS_BY_PANEL
        );
      });

      /**
       * Build parent child relationships for hematopoietic system. hematopoietic system has "UI" children in both organ
       * and tissue panels.
       * - Parents: -
       * - Children: blood non-specific, umbilical cord blood, venous blood
       */
      it("builds parent child relationships for hematopoietic system", () => {
        const uiNode = uiNodesByCategoryValueId?.get(
          INFERRED_HEMATOPOIETIC_SYSTEM
        );

        const uiParents = uiNode?.uiParents;
        expect(uiParents?.length).toEqual(0);

        const uiChildren = uiNode?.uiChildren;
        expect(uiChildren?.length).toEqual(4);
        expect(uiChildren?.includes(INFERRED_BLOOD)).toBeTruthy();
        expect(uiChildren?.includes(INFERRED_BONE_MARROW)).toBeTruthy();
        expect(uiChildren?.includes(INFERRED_SPLEEN)).toBeTruthy();
        expect(uiChildren?.includes(INFERRED_THYYMUS)).toBeTruthy();
      });

      /**
       * Build parent child relationships for renal system. renal system only has organs as direct "UI" children.
       * - Parents: -
       * - Children: kidney, urethra, bladder lumen, ureter
       */
      it("builds parent child relationships for renal system", () => {
        const uiNode = uiNodesByCategoryValueId?.get(INFERRED_RENAL_SYSTEM);

        const uiParents = uiNode?.uiParents;
        expect(uiParents?.length).toEqual(0);

        const uiChildren = uiNode?.uiChildren;
        expect(uiChildren?.length).toEqual(4);
        expect(uiChildren?.includes(INFERRED_KIDNEY)).toBeTruthy();
        expect(uiChildren?.includes(EXPLICIT_URETHRA)).toBeTruthy();
        expect(uiChildren?.includes(EXPLICIT_BLADDER_LUMEN)).toBeTruthy();
        expect(uiChildren?.includes(EXPLICIT_URETER)).toBeTruthy();
      });

      /**
       * Build parent child relationships for blood. blood has a combination of a non-specific values as well as is_a
       * or part_of values.
       * - Parents: renal system
       * - Children: blood non-specific, umbilical cord blood, venous blood
       */
      it("builds parent child relationships for blood", () => {
        const uiNode = uiNodesByCategoryValueId?.get(INFERRED_BLOOD);

        const uiParents = uiNode?.uiParents;
        expect(uiParents?.length).toEqual(1);
        expect(uiParents?.includes(INFERRED_HEMATOPOIETIC_SYSTEM)).toBeTruthy();

        const uiChildren = uiNode?.uiChildren;
        expect(uiChildren?.length).toEqual(3);
        expect(uiChildren?.includes(EXPLICIT_BLOOD)).toBeTruthy();
        expect(
          uiChildren?.includes(EXPLICIT_UMBILICAL_CORD_BLOOD)
        ).toBeTruthy();
        expect(uiChildren?.includes(EXPLICIT_VENOUS_BLOOD)).toBeTruthy();
      });

      /**
       * Build parent child relationships for bladder lumen. bladder lumen has "UI" parents in both the organ and
       * system panels.
       * - Parents: renal system, bladder organ
       * - Children: -
       */
      it("builds parent child relationships for bladder lumen", () => {
        const uiNode = uiNodesByCategoryValueId?.get(EXPLICIT_BLADDER_LUMEN);

        const uiParents = uiNode?.uiParents;
        expect(uiParents?.length).toEqual(2);
        expect(uiParents?.includes(INFERRED_RENAL_SYSTEM)).toBeTruthy();
        expect(uiParents?.includes(INFERRED_BLADDER_ORGAN)).toBeTruthy();

        const uiChildren = uiNode?.uiChildren;
        expect(uiChildren?.length).toEqual(0);
      });

      /**
       * Build parent child relationships for ureter. ureter only has systems as direct "UI" parents.
       * - Parents: renal system
       * - Children: -
       */
      it("builds parent child relationships for ureter", () => {
        const uiNode = uiNodesByCategoryValueId?.get(EXPLICIT_URETER);

        const uiParents = uiNode?.uiParents;
        expect(uiParents?.length).toEqual(1);
        expect(uiParents?.includes(INFERRED_RENAL_SYSTEM)).toBeTruthy();

        const uiChildren = uiNode?.uiChildren;
        expect(uiChildren?.length).toEqual(0);
      });
    });
  });

  // TODO(cc) remove skip
  describe.skip("Curated Ontology Category", () => {
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
          idToDeselect,
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
          idToDeselect,
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
          idToDeselect,
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
          idToDeselect,
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
        const idToSelect = ONTOLOGY_ID_HUMAN_CARNEGIE_CS1;
        const nextFilters = buildNextOntologyCategoryFilters(
          CATEGORY_FILTER_ID.DEVELOPMENT_STAGE,
          idToSelect,
          [], // No filters selected
          CATEGORY_VALUE_KEYS,
          (
            CATEGORY_FILTER_CONFIGS_BY_ID[
              CATEGORY_FILTER_ID.DEVELOPMENT_STAGE
            ] as CuratedOntologyCategoryFilterConfig
          ).mask
        );
        expect(nextFilters.length).toEqual(1);
        expect(nextFilters[0]).toEqual(idToSelect);
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
          ONTOLOGY_ID_HUMAN_INFANT,
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
          ONTOLOGY_ID_HUMAN_CHILD,
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
          ONTOLOGY_ID_HUMAN_CARNEGIE_CS1,
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
          ONTOLOGY_ID_HUMAN_EMBRYONIC_HUMAN,
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
          ONTOLOGY_ID_HUMAN_EMBRYONIC_HUMAN,
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
          ONTOLOGY_ID_HUMAN_EMBRYONIC_HUMAN,
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
          ONTOLOGY_ID_HUMAN_CARNEGIE_CS1,
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
