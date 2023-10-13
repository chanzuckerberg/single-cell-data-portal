/**
 * Test suite for category filter hook.
 *
 * Multi-Panel Tests
 * *****************
 *
 * The majority of the multi-panel tests are based on the following tissue-specific hierarchy:
 *
 * -- bladder organ
 * ---- bladder lumen
 *
 * hematopoietic system
 * -- blood
 * ---- blood, non-specific
 * ---- umbilical cord blood
 * ---- venous blood
 * -- bone marrow
 * ---- bone marrow, non-specific
 * -- spleen
 * ---- spleen, non-specific
 * -- thymus
 * ---- thymus, non-specific
 *
 * immune system
 * -- bone marrow
 * ---- bone marrow, non-specific
 * -- lymph node
 * -- spleen
 * ---- spleen, non-specific
 * -- thymus
 * ---- thymus, non-specific
 *
 * renal system
 * ---- ureter
 * ---- urethra
 * ---- bladder lumen
 * -- kidney
 *
 * For tests specific to cell type, the following hierarchy is used:
 *
 * hematopoietic cell
 * -- leukocyte
 * -- T cell
 * ------ alpha-beta T cell
 *
 */

// App dependencies
import { expect } from "@playwright/test";
import { Row } from "react-table";
import { buildNextOntologyCategoryFilters } from "src/common/hooks/useCategoryFilter/common/curatedOntologyUtils";
import {
  buildInferredOntologyTermId,
  buildMultiPanelCategoryView,
  buildUINodesByCategoryValueId,
  keyCategoryValueIdsByPanel,
  listMultiPanelSelectedViews,
  listPartiallySelectedCategoryValues,
  onRemoveMultiPanelCategoryValueTag,
  overrideSelectedParents,
} from "src/common/hooks/useCategoryFilter/common/multiPanelOntologyUtils";
import {
  CATEGORY_FILTER_CONFIGS_BY_ID,
  CELL_TYPE_DESCENDANTS,
  TISSUE_DESCENDANTS,
  TISSUE_SYSTEM_ONTOLOGY_TERM_SET,
} from "src/components/common/Filter/common/constants";
import {
  CategoryValueId,
  CATEGORY_FILTER_ID,
  CuratedOntologyCategoryFilterConfig,
  DatasetRow,
  MultiPanelCategoryFilterUIState,
  MultiPanelOntologyCategoryView,
  MultiPanelOntologyFilterConfig,
  MultiPanelUINode,
  MultiPanelUIState,
  OrFilterPrefix,
  SelectCategoryValue,
  SelectCategoryValueView,
} from "src/components/common/Filter/common/entities";
import { test } from "tests/common/test";

const { describe } = test;

describe("useCategoryFilter", () => {
  describe("Multi-Panel Category", () => {
    const TERM_ID_BLADDER_LUMEN = "UBERON:0009958";
    const TERM_ID_BLADDER_ORGAN = "UBERON:0018707";
    const TERM_ID_BLOOD = "UBERON:0000178";
    const TERM_ID_BONE_MARROW = "UBERON:0002371";
    const TERM_ID_HEMATOPOIETIC_SYSTEM = "UBERON:0002390";
    const TERM_ID_IMMUNE_SYSTEM = "UBERON:0002405";
    const TERM_ID_KIDNEY = "UBERON:0002113";
    const TERM_ID_LYMPH_NODE = "UBERON:0000029";
    const TERM_ID_SPLEEN = "UBERON:0002106";
    const TERM_ID_RENAL_MEDULLA = "UBERON:0000362";
    const TERM_ID_RENAL_SYSTEM = "UBERON:0001008";
    const TERM_ID_THORACIC_LYMPH_NODE = "UBERON:0007644";
    const TERM_ID_THYMUS = "UBERON:0002370";
    const TERM_ID_UMBILICAL_CORD_BLOOD = "UBERON:0012168";
    const TERM_ID_URETER = "UBERON:0000056";
    const TERM_ID_URETHRA = "UBERON:0000057";
    const TERM_ID_VENOUS_BLOOD = "UBERON:0013756";

    const INFERRED_BLADDER_ORGAN = `${OrFilterPrefix.INFERRED}:${TERM_ID_BLADDER_ORGAN}`;
    const INFERRED_BLOOD = `${OrFilterPrefix.INFERRED}:${TERM_ID_BLOOD}`;
    const INFERRED_BONE_MARROW = `${OrFilterPrefix.INFERRED}:${TERM_ID_BONE_MARROW}`;
    const INFERRED_HEMATOPOIETIC_SYSTEM = `${OrFilterPrefix.INFERRED}:${TERM_ID_HEMATOPOIETIC_SYSTEM}`;
    const INFERRED_IMMUNE_SYSTEM = `${OrFilterPrefix.INFERRED}:${TERM_ID_IMMUNE_SYSTEM}`;
    const INFERRED_KIDNEY = `${OrFilterPrefix.INFERRED}:${TERM_ID_KIDNEY}`;
    const INFERRED_LYMPH_NODE = `${OrFilterPrefix.INFERRED}:${TERM_ID_LYMPH_NODE}`;
    const INFERRED_RENAL_SYSTEM = `${OrFilterPrefix.INFERRED}:${TERM_ID_RENAL_SYSTEM}`;
    const INFERRED_SPLEEN = `${OrFilterPrefix.INFERRED}:${TERM_ID_SPLEEN}`;
    const INFERRED_THYMUS = `${OrFilterPrefix.INFERRED}:${TERM_ID_THYMUS}`;

    const EXPLICIT_BLADDER_LUMEN = `${OrFilterPrefix.EXPLICIT}:${TERM_ID_BLADDER_LUMEN}`;
    const EXPLICIT_BLOOD = `${OrFilterPrefix.EXPLICIT}:${TERM_ID_BLOOD}`;
    const EXPLICIT_BONE_MARROW = `${OrFilterPrefix.EXPLICIT}:${TERM_ID_BONE_MARROW}`;
    const EXPLICIT_KIDNEY = `${OrFilterPrefix.EXPLICIT}:${TERM_ID_KIDNEY}`;
    const EXPLICIT_RENAL_MEDULLA = `${OrFilterPrefix.EXPLICIT}:${TERM_ID_RENAL_MEDULLA}`;
    const EXPLICIT_SPLEEN = `${OrFilterPrefix.EXPLICIT}:${TERM_ID_SPLEEN}`;
    const EXPLICIT_THORACIC_LYMPH_NODE = `${OrFilterPrefix.EXPLICIT}:${TERM_ID_THORACIC_LYMPH_NODE}`;
    const EXPLICIT_THYMUS = `${OrFilterPrefix.EXPLICIT}:${TERM_ID_THYMUS}`;
    const EXPLICIT_UMBILICAL_CORD_BLOOD = `${OrFilterPrefix.EXPLICIT}:${TERM_ID_UMBILICAL_CORD_BLOOD}`;
    const EXPLICIT_URETER = `${OrFilterPrefix.EXPLICIT}:${TERM_ID_URETER}`;
    const EXPLICIT_URETHRA = `${OrFilterPrefix.EXPLICIT}:${TERM_ID_URETHRA}`;
    const EXPLICIT_VENOUS_BLOOD = `${OrFilterPrefix.EXPLICIT}:${TERM_ID_VENOUS_BLOOD}`;

    const INFERRED_HEMATOPOIETIC_SYSTEM_CATEGORY_VALUE = {
      categoryValueId: INFERRED_HEMATOPOIETIC_SYSTEM,
      count: 0,
      selected: false,
      selectedPartial: false,
    };

    const INFERRED_IMMUNE_SYSTEM_CATEGORY_VALUE = {
      categoryValueId: INFERRED_IMMUNE_SYSTEM,
      count: 0,
      selected: false,
      selectedPartial: false,
    };

    const INFERRED_RENAL_SYSTEM_CATEGORY_VALUE = {
      categoryValueId: INFERRED_RENAL_SYSTEM,
      count: 0,
      selected: false,
      selectedPartial: false,
    };

    const INFERRED_BLADDER_ORGAN_CATEGORY_VALUE = {
      categoryValueId: INFERRED_BLADDER_ORGAN,
      count: 0,
      selected: false,
      selectedPartial: false,
    };

    const INFERRED_BLOOD_CATEGORY_VALUE = {
      categoryValueId: INFERRED_BLOOD,
      count: 0,
      selected: false,
      selectedPartial: false,
    };

    const INFERRED_BONE_MARROW_CATEGORY_VALUE = {
      categoryValueId: INFERRED_BONE_MARROW,
      count: 0,
      selected: false,
      selectedPartial: false,
    };

    const INFERRED_KIDNEY_CATEGORY_VALUE = {
      categoryValueId: INFERRED_KIDNEY,
      count: 0,
      selected: false,
      selectedPartial: false,
    };

    const INFERRED_LYMPH_NODE_CATEGORY_VALUE = {
      categoryValueId: INFERRED_LYMPH_NODE,
      count: 0,
      selected: false,
      selectedPartial: false,
    };

    const INFERRED_SPLEEN_CATEGORY_VALUE = {
      categoryValueId: INFERRED_SPLEEN,
      count: 0,
      selected: false,
      selectedPartial: false,
    };

    const INFERRED_THYMUS_CATEGORY_VALUE = {
      categoryValueId: INFERRED_THYMUS,
      count: 0,
      selected: false,
      selectedPartial: false,
    };

    const EXPLICIT_BLADDER_LUMEN_CATEGORY_VALUE = {
      categoryValueId: EXPLICIT_BLADDER_LUMEN,
      count: 0,
      selected: false,
      selectedPartial: false,
    };

    const EXPLICIT_BLOOD_CATEGORY_VALUE = {
      categoryValueId: EXPLICIT_BLOOD,
      count: 0,
      selected: false,
      selectedPartial: false,
    };

    const EXPLICIT_BONE_MARROW_CATEGORY_VALUE = {
      categoryValueId: EXPLICIT_BONE_MARROW,
      count: 0,
      selected: false,
      selectedPartial: false,
    };

    const EXPLICIT_KIDNEY_CATEGORY_VALUE = {
      categoryValueId: EXPLICIT_KIDNEY,
      count: 0,
      selected: false,
      selectedPartial: false,
    };

    const EXPLICIT_RENAL_MEDULLA_CATEGORY_VALUE = {
      categoryValueId: EXPLICIT_RENAL_MEDULLA,
      count: 0,
      selected: false,
      selectedPartial: false,
    };

    const EXPLICIT_SPLEEN_CATEGORY_VALUE = {
      categoryValueId: EXPLICIT_SPLEEN,
      count: 0,
      selected: false,
      selectedPartial: false,
    };

    const EXPLICIT_THYMUS_CATEGORY_VALUE = {
      categoryValueId: EXPLICIT_THYMUS,
      count: 0,
      selected: false,
      selectedPartial: false,
    };

    const EXPLICIT_UMBILICAL_CORD_BLOOD_CATEGORY_VALUE = {
      categoryValueId: EXPLICIT_UMBILICAL_CORD_BLOOD,
      count: 0,
      selected: false,
      selectedPartial: false,
    };

    const EXPLICIT_URETER_CATEGORY_VALUE = {
      categoryValueId: EXPLICIT_URETER,
      count: 0,
      selected: false,
      selectedPartial: false,
    };

    const EXPLICIT_URETHRA_CATEGORY_VALUE = {
      categoryValueId: EXPLICIT_URETHRA,
      count: 0,
      selected: false,
      selectedPartial: false,
    };

    const EXPLICIT_VENOUS_BLOOD_CATEGORY_VALUE = {
      categoryValueId: EXPLICIT_VENOUS_BLOOD,
      count: 0,
      selected: false,
      selectedPartial: false,
    };

    const KEYED_CATEGORY_VALUES = new Map<CategoryValueId, SelectCategoryValue>(
      [
        [
          INFERRED_HEMATOPOIETIC_SYSTEM,
          INFERRED_HEMATOPOIETIC_SYSTEM_CATEGORY_VALUE,
        ],
        [INFERRED_IMMUNE_SYSTEM, INFERRED_IMMUNE_SYSTEM_CATEGORY_VALUE],
        [INFERRED_RENAL_SYSTEM, INFERRED_RENAL_SYSTEM_CATEGORY_VALUE],
        [INFERRED_BLADDER_ORGAN, INFERRED_BLADDER_ORGAN_CATEGORY_VALUE],
        [INFERRED_BLOOD, INFERRED_BLOOD_CATEGORY_VALUE],
        [INFERRED_BONE_MARROW, INFERRED_BONE_MARROW_CATEGORY_VALUE],
        [INFERRED_KIDNEY, INFERRED_KIDNEY_CATEGORY_VALUE],
        [INFERRED_LYMPH_NODE, INFERRED_LYMPH_NODE_CATEGORY_VALUE],
        [INFERRED_SPLEEN, INFERRED_SPLEEN_CATEGORY_VALUE],
        [INFERRED_THYMUS, INFERRED_THYMUS_CATEGORY_VALUE],
        [EXPLICIT_BLADDER_LUMEN, EXPLICIT_BLADDER_LUMEN_CATEGORY_VALUE],
        [EXPLICIT_BLOOD, EXPLICIT_BLOOD_CATEGORY_VALUE],
        [EXPLICIT_BONE_MARROW, EXPLICIT_BONE_MARROW_CATEGORY_VALUE],
        [EXPLICIT_KIDNEY, EXPLICIT_KIDNEY_CATEGORY_VALUE],
        [EXPLICIT_RENAL_MEDULLA, EXPLICIT_RENAL_MEDULLA_CATEGORY_VALUE],
        [EXPLICIT_SPLEEN, EXPLICIT_SPLEEN_CATEGORY_VALUE],
        [EXPLICIT_THYMUS, EXPLICIT_THYMUS_CATEGORY_VALUE],
        [
          EXPLICIT_UMBILICAL_CORD_BLOOD,
          EXPLICIT_UMBILICAL_CORD_BLOOD_CATEGORY_VALUE,
        ],
        [EXPLICIT_URETER, EXPLICIT_URETER_CATEGORY_VALUE],
        [EXPLICIT_URETHRA, EXPLICIT_URETHRA_CATEGORY_VALUE],
        [EXPLICIT_VENOUS_BLOOD, EXPLICIT_VENOUS_BLOOD_CATEGORY_VALUE],
      ]
    );

    const PANEL_INDEX_SYSTEM = 0;
    const PANEL_INDEX_ORGAN = 1;
    const PANEL_INDEX_TISSUE = 2;

    const TISSUE_SYSTEMS = [
      INFERRED_HEMATOPOIETIC_SYSTEM,
      INFERRED_IMMUNE_SYSTEM,
      INFERRED_RENAL_SYSTEM,
    ];
    const TISSUE_ORGANS = [
      INFERRED_BLADDER_ORGAN,
      INFERRED_BLOOD,
      INFERRED_BONE_MARROW,
      INFERRED_KIDNEY,
      INFERRED_LYMPH_NODE,
      INFERRED_SPLEEN,
      INFERRED_THYMUS,
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

    const TISSUE_CATEGORY_VALUE_IDS_BY_PANEL = [
      TISSUE_SYSTEMS,
      TISSUE_ORGANS,
      TISSUES,
    ];

    const INFERRED_HEMATOPOIETIC_SYSTEM_UI_NODE = {
      categoryValueId: INFERRED_HEMATOPOIETIC_SYSTEM,
      panelIndex: PANEL_INDEX_SYSTEM,
      uiChildren: [
        INFERRED_BLOOD,
        INFERRED_BONE_MARROW,
        INFERRED_SPLEEN,
        INFERRED_THYMUS,
      ],
      uiParents: [],
    };

    const INFERRED_IMMUNE_SYSTEM_UI_NODE = {
      categoryValueId: INFERRED_IMMUNE_SYSTEM,
      panelIndex: PANEL_INDEX_SYSTEM,
      // Note, children list is not exact but has enough values for exercising different cases
      uiChildren: [
        INFERRED_BONE_MARROW,
        INFERRED_SPLEEN,
        INFERRED_THYMUS,
        EXPLICIT_THORACIC_LYMPH_NODE,
      ],
      uiParents: [],
    };

    const INFERRED_BLOOD_UI_NODE = {
      categoryValueId: INFERRED_BLOOD,
      panelIndex: PANEL_INDEX_ORGAN,
      uiChildren: [
        EXPLICIT_BLOOD,
        EXPLICIT_UMBILICAL_CORD_BLOOD,
        EXPLICIT_VENOUS_BLOOD,
      ],
      uiParents: [INFERRED_HEMATOPOIETIC_SYSTEM],
    };

    const INFERRED_BONE_MARROW_UI_NODE = {
      categoryValueId: INFERRED_BONE_MARROW,
      panelIndex: PANEL_INDEX_ORGAN,
      uiChildren: [EXPLICIT_BONE_MARROW],
      uiParents: [INFERRED_HEMATOPOIETIC_SYSTEM, INFERRED_IMMUNE_SYSTEM],
    };

    const INFERRED_LYMPH_NODE_UI_NODE = {
      categoryValueId: INFERRED_LYMPH_NODE,
      panelIndex: PANEL_INDEX_ORGAN,
      uiChildren: [], // Test data here does not match UBERON
      uiParents: [INFERRED_IMMUNE_SYSTEM],
    };

    const INFERRED_SPLEEN_UI_NODE = {
      categoryValueId: INFERRED_SPLEEN,
      panelIndex: PANEL_INDEX_ORGAN,
      uiChildren: [EXPLICIT_SPLEEN],
      uiParents: [INFERRED_HEMATOPOIETIC_SYSTEM, INFERRED_IMMUNE_SYSTEM],
    };

    const INFERRED_THYMUS_UI_NODE = {
      categoryValueId: INFERRED_THYMUS,
      panelIndex: PANEL_INDEX_ORGAN,
      uiChildren: [EXPLICIT_THYMUS],
      uiParents: [INFERRED_HEMATOPOIETIC_SYSTEM, INFERRED_IMMUNE_SYSTEM],
    };

    const EXPLICIT_BLOOD_UI_NODE = {
      categoryValueId: EXPLICIT_BLOOD,
      panelIndex: PANEL_INDEX_TISSUE,
      uiChildren: [],
      uiParents: [INFERRED_BLOOD],
    };

    const EXPLICIT_BONE_MARROW_UI_NODE = {
      categoryValueId: EXPLICIT_BONE_MARROW,
      panelIndex: PANEL_INDEX_TISSUE,
      uiChildren: [],
      uiParents: [INFERRED_BONE_MARROW],
    };

    const EXPLICIT_SPLEEN_UI_NODE = {
      categoryValueId: EXPLICIT_SPLEEN,
      panelIndex: PANEL_INDEX_TISSUE,
      uiChildren: [],
      uiParents: [INFERRED_SPLEEN],
    };

    const EXPLICIT_THYMUS_UI_NODE = {
      categoryValueId: EXPLICIT_THYMUS,
      panelIndex: PANEL_INDEX_TISSUE,
      uiChildren: [],
      uiParents: [INFERRED_THYMUS],
    };

    const EXPLICIT_UMBILICAL_CORD_BLOOD_UI_NODE = {
      categoryValueId: EXPLICIT_UMBILICAL_CORD_BLOOD,
      panelIndex: PANEL_INDEX_TISSUE,
      uiChildren: [],
      uiParents: [INFERRED_BLOOD],
    };

    const EXPLICIT_VENOUS_BLOOD_UI_NODE = {
      categoryValueId: EXPLICIT_VENOUS_BLOOD,
      panelIndex: PANEL_INDEX_TISSUE,
      uiChildren: [],
      uiParents: [INFERRED_BLOOD],
    };

    const EXPLICIT_THORACIC_LYMPH_NODE_UI_NODE = {
      categoryValueId: EXPLICIT_THORACIC_LYMPH_NODE,
      panelIndex: PANEL_INDEX_TISSUE,
      uiChildren: [],
      uiParents: [INFERRED_IMMUNE_SYSTEM],
    };

    const UI_NODES_BY_CATEGORY_VALUE_ID = new Map<
      CategoryValueId,
      MultiPanelUINode
    >([
      [INFERRED_HEMATOPOIETIC_SYSTEM, INFERRED_HEMATOPOIETIC_SYSTEM_UI_NODE],
      [INFERRED_IMMUNE_SYSTEM, INFERRED_IMMUNE_SYSTEM_UI_NODE],
      [INFERRED_BLOOD, INFERRED_BLOOD_UI_NODE],
      [INFERRED_BONE_MARROW, INFERRED_BONE_MARROW_UI_NODE],
      [INFERRED_LYMPH_NODE, INFERRED_LYMPH_NODE_UI_NODE],
      [INFERRED_SPLEEN, INFERRED_SPLEEN_UI_NODE],
      [INFERRED_THYMUS, INFERRED_THYMUS_UI_NODE],
      [EXPLICIT_BLOOD, EXPLICIT_BLOOD_UI_NODE],
      [EXPLICIT_BONE_MARROW, EXPLICIT_BONE_MARROW_UI_NODE],
      [EXPLICIT_SPLEEN, EXPLICIT_SPLEEN_UI_NODE],
      [EXPLICIT_THYMUS, EXPLICIT_THYMUS_UI_NODE],
      [EXPLICIT_UMBILICAL_CORD_BLOOD, EXPLICIT_UMBILICAL_CORD_BLOOD_UI_NODE],
      [EXPLICIT_VENOUS_BLOOD, EXPLICIT_VENOUS_BLOOD_UI_NODE],
      [EXPLICIT_THORACIC_LYMPH_NODE, EXPLICIT_THORACIC_LYMPH_NODE_UI_NODE],
    ]);

    const BASE_CATEGORY_FILTER_UI_STATE: MultiPanelCategoryFilterUIState = {
      selected: [],
      selectedPartial: [],
      uiNodesByCategoryValueId: UI_NODES_BY_CATEGORY_VALUE_ID,
    };

    const INFERRED_HEMATOPOIETIC_CATEGORY_VALUE_VIEW = {
      categoryValueId: INFERRED_HEMATOPOIETIC_SYSTEM,
      count: 0,
      label: "hematopoietic system",
      selected: false,
      selectedPartial: false,
      value: INFERRED_HEMATOPOIETIC_SYSTEM,
      visible: true,
    };

    const INFERRED_IMMUNE_CATEGORY_VALUE_VIEW = {
      categoryValueId: INFERRED_IMMUNE_SYSTEM,
      count: 0,
      label: "immune system",
      selected: false,
      selectedPartial: false,
      value: INFERRED_IMMUNE_SYSTEM,
      visible: true,
    };

    const INFERRED_BLOOD_CATEGORY_VALUE_VIEW = {
      categoryValueId: INFERRED_BLOOD,
      count: 0,
      label: "blood",
      selected: false,
      selectedPartial: false,
      value: INFERRED_BLOOD,
      visible: true,
    };

    const INFERRED_BONE_MARROW_CATEGORY_VALUE_VIEW = {
      categoryValueId: INFERRED_BONE_MARROW,
      count: 0,
      label: "bone marrow",
      selected: false,
      selectedPartial: false,
      value: INFERRED_BONE_MARROW,
      visible: true,
    };

    const INFERRED_LYMPH_NODE_CATEGORY_VALUE_VIEW = {
      categoryValueId: INFERRED_LYMPH_NODE,
      count: 0,
      label: "lymph node",
      selected: false,
      selectedPartial: false,
      value: INFERRED_LYMPH_NODE,
      visible: true,
    };

    const INFERRED_SPLEEN_CATEGORY_VALUE_VIEW = {
      categoryValueId: INFERRED_SPLEEN,
      count: 0,
      label: "spleen",
      selected: false,
      selectedPartial: false,
      value: INFERRED_SPLEEN,
      visible: true,
    };

    const INFERRED_THYMUS_CATEGORY_VALUE_VIEW = {
      categoryValueId: INFERRED_THYMUS,
      count: 0,
      label: "thymus",
      selected: false,
      selectedPartial: false,
      value: INFERRED_THYMUS,
      visible: true,
    };

    const EXPLICIT_BLOOD_CATEGORY_VALUE_VIEW = {
      categoryValueId: EXPLICIT_BLOOD,
      count: 0,
      label: "blood, non-specific",
      selected: false,
      selectedPartial: false,
      value: EXPLICIT_BLOOD,
      visible: true,
    };

    const EXPLICIT_BONE_MARROW_CATEGORY_VALUE_VIEW = {
      categoryValueId: EXPLICIT_BONE_MARROW,
      count: 0,
      label: "bone marrow, non-specific",
      selected: false,
      selectedPartial: false,
      value: EXPLICIT_BONE_MARROW,
      visible: true,
    };

    const EXPLICIT_SPLEEN_CATEGORY_VALUE_VIEW = {
      categoryValueId: EXPLICIT_SPLEEN,
      count: 0,
      label: "spleen, non-specific",
      selected: false,
      selectedPartial: false,
      value: EXPLICIT_SPLEEN,
      visible: true,
    };

    const EXPLICIT_THYMUS_CATEGORY_VALUE_VIEW = {
      categoryValueId: EXPLICIT_THYMUS,
      count: 0,
      label: "thymus, non-specific",
      selected: false,
      selectedPartial: false,
      value: EXPLICIT_THYMUS,
      visible: true,
    };

    const EXPLICIT_UMBILICAL_CORD_CATEGORY_VALUE_VIEW = {
      categoryValueId: EXPLICIT_UMBILICAL_CORD_BLOOD,
      count: 0,
      label: "umbilical cord blood",
      selected: false,
      selectedPartial: false,
      value: EXPLICIT_UMBILICAL_CORD_BLOOD,
      visible: true,
    };

    const EXPLICIT_VENOUS_BLOOD_CATEGORY_VALUE_VIEW = {
      categoryValueId: EXPLICIT_VENOUS_BLOOD,
      count: 0,
      label: "venous blood",
      selected: false,
      selectedPartial: false,
      value: EXPLICIT_VENOUS_BLOOD,
      visible: true,
    };

    const EXPLICIT_THORACIC_LYMPH_NODE_VALUE_VIEW = {
      categoryValueId: EXPLICIT_THORACIC_LYMPH_NODE,
      count: 0,
      label: "thoracic lymph node",
      selected: false,
      selectedPartial: false,
      value: EXPLICIT_THORACIC_LYMPH_NODE,
      visible: true,
    };

    const VALUE_VIEWS_BY_CATEGORY_VALUE_ID = new Map<
      CategoryValueId,
      SelectCategoryValueView
    >([
      [
        INFERRED_HEMATOPOIETIC_SYSTEM,
        INFERRED_HEMATOPOIETIC_CATEGORY_VALUE_VIEW,
      ],
      [INFERRED_IMMUNE_SYSTEM, INFERRED_IMMUNE_CATEGORY_VALUE_VIEW],
      [INFERRED_BLOOD, INFERRED_BLOOD_CATEGORY_VALUE_VIEW],
      [INFERRED_BONE_MARROW, INFERRED_BONE_MARROW_CATEGORY_VALUE_VIEW],
      [INFERRED_LYMPH_NODE, INFERRED_LYMPH_NODE_CATEGORY_VALUE_VIEW],
      [INFERRED_SPLEEN, INFERRED_SPLEEN_CATEGORY_VALUE_VIEW],
      [INFERRED_THYMUS, INFERRED_THYMUS_CATEGORY_VALUE_VIEW],
      [EXPLICIT_BLOOD, EXPLICIT_BLOOD_CATEGORY_VALUE_VIEW],
      [EXPLICIT_BONE_MARROW, EXPLICIT_BONE_MARROW_CATEGORY_VALUE_VIEW],
      [EXPLICIT_SPLEEN, EXPLICIT_SPLEEN_CATEGORY_VALUE_VIEW],
      [EXPLICIT_THYMUS, EXPLICIT_THYMUS_CATEGORY_VALUE_VIEW],
      [
        EXPLICIT_UMBILICAL_CORD_BLOOD,
        EXPLICIT_UMBILICAL_CORD_CATEGORY_VALUE_VIEW,
      ],
      [EXPLICIT_VENOUS_BLOOD, EXPLICIT_VENOUS_BLOOD_CATEGORY_VALUE_VIEW],
      [EXPLICIT_THORACIC_LYMPH_NODE, EXPLICIT_THORACIC_LYMPH_NODE_VALUE_VIEW],
    ]);

    const TERM_ID_HEMATOPOIETIC_CELL = "CL:0000988";
    const TERM_ID_LEUKOCYTE = "CL:0000738";
    const TERM_ID_T_CELL = "CL:0000084";
    const TERM_ID_ALPHA_BETA_T_CELL = "CL:0000789";

    const INFERRED_HEMATOPOIETIC_CELL = `${OrFilterPrefix.INFERRED}:${TERM_ID_HEMATOPOIETIC_CELL}`;
    const INFERRED_LEUKOCYTE = `${OrFilterPrefix.INFERRED}:${TERM_ID_LEUKOCYTE}`;
    const INFERRED_T_CELL = `${OrFilterPrefix.INFERRED}:${TERM_ID_T_CELL}`;

    const EXPLICIT_ALPHA_BETA_T_CELL = `${OrFilterPrefix.EXPLICIT}:${TERM_ID_ALPHA_BETA_T_CELL}`;

    const CELL_CLASSES = [INFERRED_HEMATOPOIETIC_CELL];
    const CELL_SUBCLASSES = [INFERRED_LEUKOCYTE, INFERRED_T_CELL];
    const CELL_TYPES = [EXPLICIT_ALPHA_BETA_T_CELL];

    const CELL_TYPE_CATEGORY_VALUE_IDS_BY_PANEL = [
      CELL_CLASSES,
      CELL_SUBCLASSES,
      CELL_TYPES,
    ];

    describe("keyCategoryValueIdsByPanel", () => {
      /**
       * Panels: curated, curated, explicit.
       */
      test("keys panels", () => {
        const config =
          CATEGORY_FILTER_CONFIGS_BY_ID[CATEGORY_FILTER_ID.TISSUE_CALCULATED];
        const categoryValueIdsByPanel = keyCategoryValueIdsByPanel(
          config as MultiPanelOntologyFilterConfig,
          []
        );
        expect(categoryValueIdsByPanel.length).toEqual(3);
      });

      /**
       * Curated panel: system
       */
      test("keys curated panel - Curated panel: system", () => {
        const config =
          CATEGORY_FILTER_CONFIGS_BY_ID[CATEGORY_FILTER_ID.TISSUE_CALCULATED];
        const categoryValueIdsByPanel = keyCategoryValueIdsByPanel(
          config as MultiPanelOntologyFilterConfig,
          []
        );

        const systemPanel = categoryValueIdsByPanel[0];
        expect(systemPanel).toBeTruthy();

        const systems = TISSUE_SYSTEM_ONTOLOGY_TERM_SET["UBERON"];
        // eslint-disable-next-line @typescript-eslint/no-non-null-assertion -- truthy check above
        expect(systemPanel.length).toEqual(systems!.length);

        // eslint-disable-next-line @typescript-eslint/no-non-null-assertion -- truthy check above
        const isEverySystemIncluded = systems!.every((system) =>
          systemPanel.includes(
            buildInferredOntologyTermId(system.ontology_term_id)
          )
        );
        expect(isEverySystemIncluded).toBeTruthy();
      });

      /**
       * Explicit panel: tissue
       */
      test("keys curated panel - Explicit panel: tissue", () => {
        const config =
          CATEGORY_FILTER_CONFIGS_BY_ID[CATEGORY_FILTER_ID.TISSUE_CALCULATED];
        const categoryValueIdsByPanel = keyCategoryValueIdsByPanel(
          config as MultiPanelOntologyFilterConfig,
          [
            {
              original: { tissueCalculated: [EXPLICIT_BLOOD] },
            } as unknown as Row<DatasetRow>,
            {
              original: { tissueCalculated: [EXPLICIT_BLOOD, EXPLICIT_SPLEEN] },
            } as unknown as Row<DatasetRow>,
          ]
        );

        const tissuePanel = categoryValueIdsByPanel[2];
        expect(tissuePanel).toBeTruthy();

        expect(tissuePanel.length).toEqual(2);
        expect(tissuePanel.includes(EXPLICIT_BLOOD)).toBeTruthy();
        expect(tissuePanel.includes(EXPLICIT_SPLEEN)).toBeTruthy();
      });
    });

    describe("overrideSelectedParents", () => {
      let uiNodesByCategoryValueId: Map<CategoryValueId, MultiPanelUINode>;
      test.beforeAll(() => {
        uiNodesByCategoryValueId = buildUINodesByCategoryValueId(
          TISSUE_CATEGORY_VALUE_IDS_BY_PANEL,
          TISSUE_DESCENDANTS
        );
      });

      /**
       * Selected: blood non-specific
       * Post-overrides: blood non-specific
       */
      test("doesn't override tissue", () => {
        const overriddenSelectedValues = overrideSelectedParents(
          [EXPLICIT_BLOOD],
          TISSUE_DESCENDANTS,
          uiNodesByCategoryValueId
        );
        expect(overriddenSelectedValues.length).toEqual(1);
      });

      /**
       * Selected: blood, blood non-specific
       * Post-overrides: blood non-specific
       */
      test("overrides inferred organ with explicit tissue", () => {
        const overriddenSelectedValues = overrideSelectedParents(
          [INFERRED_BLOOD, EXPLICIT_BLOOD],
          TISSUE_DESCENDANTS,
          uiNodesByCategoryValueId
        );
        expect(overriddenSelectedValues.length).toEqual(1);
        expect(overriddenSelectedValues[0]).toEqual(EXPLICIT_BLOOD);
      });

      /**
       * Selected: blood, umbilical cord blood
       * Post-overrides: umbilical cord blood
       */
      test("overrides organ with tissue", () => {
        const overriddenSelectedValues = overrideSelectedParents(
          [INFERRED_BLOOD, EXPLICIT_UMBILICAL_CORD_BLOOD],
          TISSUE_DESCENDANTS,
          uiNodesByCategoryValueId
        );
        expect(overriddenSelectedValues.length).toEqual(1);
        expect(overriddenSelectedValues[0]).toEqual(
          EXPLICIT_UMBILICAL_CORD_BLOOD
        );
      });

      /**
       * Selected: hematopoietic system, blood
       * Post-overrides: blood
       */
      test("overrides system with organ", () => {
        const overriddenSelectedValues = overrideSelectedParents(
          [INFERRED_HEMATOPOIETIC_SYSTEM, INFERRED_BLOOD],
          TISSUE_DESCENDANTS,
          uiNodesByCategoryValueId
        );
        expect(overriddenSelectedValues.length).toEqual(1);
        expect(overriddenSelectedValues[0]).toEqual(INFERRED_BLOOD);
      });

      /**
       * Selected: hematopoietic system, blood non-specific
       * Post-overrides: blood non-specific
       */
      test("overrides system with tissue", () => {
        const overriddenSelectedValues = overrideSelectedParents(
          [INFERRED_HEMATOPOIETIC_SYSTEM, EXPLICIT_BLOOD],
          TISSUE_DESCENDANTS,
          uiNodesByCategoryValueId
        );
        expect(overriddenSelectedValues.length).toEqual(1);
        expect(overriddenSelectedValues[0]).toEqual(EXPLICIT_BLOOD);
      });

      /**
       * Selected: hematopoietic system, immune system, spleen
       * Post-overrides: spleen
       */
      test("overrides multiple systems with organ", () => {
        const overriddenSelectedValues = overrideSelectedParents(
          [
            INFERRED_HEMATOPOIETIC_SYSTEM,
            INFERRED_IMMUNE_SYSTEM,
            INFERRED_SPLEEN,
          ],
          TISSUE_DESCENDANTS,
          uiNodesByCategoryValueId
        );
        expect(overriddenSelectedValues.length).toEqual(1);
        expect(overriddenSelectedValues[0]).toEqual(INFERRED_SPLEEN);
      });

      describe("cell type", () => {
        let cellTypeUINodesByCategoryValueId: Map<
          CategoryValueId,
          MultiPanelUINode
        >;
        test.beforeAll(() => {
          cellTypeUINodesByCategoryValueId = buildUINodesByCategoryValueId(
            CELL_TYPE_CATEGORY_VALUE_IDS_BY_PANEL,
            CELL_TYPE_DESCENDANTS
          );
        });

        /**
         * Selected: leukocyte, T cell
         * Post-overrides: leukocyte, T cell
         */
        test("doesn't override within cell subclasses", () => {
          const overriddenSelectedValues = overrideSelectedParents(
            [INFERRED_LEUKOCYTE, INFERRED_T_CELL],
            CELL_TYPE_DESCENDANTS,
            cellTypeUINodesByCategoryValueId
          );
          expect(overriddenSelectedValues.length).toEqual(2);
          expect(
            overriddenSelectedValues.includes(INFERRED_LEUKOCYTE)
          ).toBeTruthy();
          expect(
            overriddenSelectedValues.includes(INFERRED_T_CELL)
          ).toBeTruthy();
        });

        /**
         * Selected: leukocyte, T cell, alpha-beta T cell
         * Post-overrides: alpha-beta T cell
         */
        test("overrides cell subclasses with selected cell type", () => {
          const overriddenSelectedValues = overrideSelectedParents(
            [INFERRED_LEUKOCYTE, INFERRED_T_CELL, EXPLICIT_ALPHA_BETA_T_CELL],
            CELL_TYPE_DESCENDANTS,
            cellTypeUINodesByCategoryValueId
          );
          expect(overriddenSelectedValues.length).toEqual(1);
          expect(
            overriddenSelectedValues.includes(EXPLICIT_ALPHA_BETA_T_CELL)
          ).toBeTruthy();
        });

        /**
         * Selected: hemtapoietic cell, leukocyte, T cell, alpha-beta T cell
         * Post-overrides: alpha-beta T cell
         */
        test("overrides cell classes and subclasses with selected cell type", () => {
          const overriddenSelectedValues = overrideSelectedParents(
            [
              INFERRED_HEMATOPOIETIC_CELL,
              INFERRED_LEUKOCYTE,
              INFERRED_T_CELL,
              EXPLICIT_ALPHA_BETA_T_CELL,
            ],
            CELL_TYPE_DESCENDANTS,
            cellTypeUINodesByCategoryValueId
          );
          expect(overriddenSelectedValues.length).toEqual(1);
          expect(
            overriddenSelectedValues.includes(EXPLICIT_ALPHA_BETA_T_CELL)
          ).toBeTruthy();
        });

        /**
         * Selected: hemtapoietic cell, alpha-beta T cell
         * Post-overrides: alpha-beta T cell
         */
        test("overrides cell classes with selected cell type", () => {
          const overriddenSelectedValues = overrideSelectedParents(
            [INFERRED_HEMATOPOIETIC_CELL, EXPLICIT_ALPHA_BETA_T_CELL],
            CELL_TYPE_DESCENDANTS,
            cellTypeUINodesByCategoryValueId
          );
          expect(overriddenSelectedValues.length).toEqual(1);
          expect(
            overriddenSelectedValues.includes(EXPLICIT_ALPHA_BETA_T_CELL)
          ).toBeTruthy();
        });
      });
    });

    describe("onRemoveMultiPanelCategoryValueTag", () => {
      /**
       * Currently selected - blood non-specific
       * Remove - blood non-specific
       * Still selected - none
       */
      test("removes tissue", () => {
        const categoryFilterUIState = {
          ...BASE_CATEGORY_FILTER_UI_STATE,
          selected: [EXPLICIT_BLOOD],
          selectedPartial: [],
        };

        const selectedCategoryValueIds = onRemoveMultiPanelCategoryValueTag(
          EXPLICIT_BLOOD,
          [],
          categoryFilterUIState,
          TISSUE_DESCENDANTS
        );
        expect(selectedCategoryValueIds.length).toEqual(0);
      });

      /**
       * Currently selected - blood
       * Remove - blood
       * Still selected - none
       */
      test("removes organ", () => {
        const categoryFilterUIState = {
          ...BASE_CATEGORY_FILTER_UI_STATE,
          selected: [INFERRED_BLOOD],
          selectedPartial: [],
        };

        const selectedCategoryValueIds = onRemoveMultiPanelCategoryValueTag(
          INFERRED_BLOOD,
          [],
          categoryFilterUIState,
          TISSUE_DESCENDANTS
        );
        expect(selectedCategoryValueIds.length).toEqual(0);
      });

      /**
       * Currently selected - hematopoietic system
       * Remove - hematopoietic system
       * Still selected - none
       */
      test("removes system", () => {
        const categoryFilterUIState = {
          ...BASE_CATEGORY_FILTER_UI_STATE,
          selected: [INFERRED_HEMATOPOIETIC_SYSTEM],
          selectedPartial: [],
        };

        const selectedCategoryValueIds = onRemoveMultiPanelCategoryValueTag(
          INFERRED_HEMATOPOIETIC_SYSTEM,
          [],
          categoryFilterUIState,
          TISSUE_DESCENDANTS
        );

        expect(selectedCategoryValueIds.length).toEqual(0);
      });

      /**
       * Currently selected - blood, blood no-specific, umbilical cord blood, venous blood
       * Remove - blood
       * Still selected - none
       */
      test("removes organ and all descendants of organ", () => {
        const categoryFilterUIState = {
          ...BASE_CATEGORY_FILTER_UI_STATE,
          selected: [
            INFERRED_BLOOD,
            EXPLICIT_BLOOD,
            EXPLICIT_UMBILICAL_CORD_BLOOD,
            EXPLICIT_VENOUS_BLOOD,
          ],
          selectedPartial: [],
        };

        const selectedCategoryValueIds = onRemoveMultiPanelCategoryValueTag(
          INFERRED_BLOOD,
          [
            EXPLICIT_BLOOD,
            EXPLICIT_UMBILICAL_CORD_BLOOD,
            EXPLICIT_VENOUS_BLOOD,
          ],
          categoryFilterUIState,
          TISSUE_DESCENDANTS
        );

        expect(selectedCategoryValueIds.length).toEqual(0);
      });

      /**
       * Currently selected - hematopoietic system, blood, bone marrow spleen, thymus
       * Remove - hematopoietic system
       * Still selected - none
       */
      test("removes system and all descendant organs of system", () => {
        const categoryFilterUIState = {
          ...BASE_CATEGORY_FILTER_UI_STATE,
          selected: [
            INFERRED_HEMATOPOIETIC_SYSTEM,
            INFERRED_BLOOD,
            INFERRED_BONE_MARROW,
            INFERRED_SPLEEN,
            INFERRED_THYMUS,
          ],
          selectedPartial: [],
        };

        const selectedCategoryValueIds = onRemoveMultiPanelCategoryValueTag(
          INFERRED_HEMATOPOIETIC_SYSTEM,
          [
            INFERRED_BLOOD,
            INFERRED_BONE_MARROW,
            INFERRED_SPLEEN,
            INFERRED_THYMUS,
          ],
          categoryFilterUIState,
          TISSUE_DESCENDANTS
        );

        expect(selectedCategoryValueIds.length).toEqual(0);
      });

      /**
       * Currently selected - hematopoietic system, bone marrow spleen, thymus, blood non-specific, umbilical cord
       * blood, venous blood
       * Remove - hematopoietic system
       * Still selected - none
       */
      test("removes system and all descendants of system", () => {
        const categoryFilterUIState = {
          ...BASE_CATEGORY_FILTER_UI_STATE,
          selected: [
            INFERRED_HEMATOPOIETIC_SYSTEM,
            INFERRED_BONE_MARROW,
            INFERRED_SPLEEN,
            INFERRED_THYMUS,
            EXPLICIT_BLOOD,
            EXPLICIT_UMBILICAL_CORD_BLOOD,
            EXPLICIT_VENOUS_BLOOD,
          ],
          selectedPartial: [],
        };

        const selectedCategoryValueIds = onRemoveMultiPanelCategoryValueTag(
          INFERRED_HEMATOPOIETIC_SYSTEM,
          [
            INFERRED_BONE_MARROW,
            INFERRED_SPLEEN,
            INFERRED_THYMUS,
            EXPLICIT_BLOOD,
            EXPLICIT_UMBILICAL_CORD_BLOOD,
            EXPLICIT_VENOUS_BLOOD,
          ],
          categoryFilterUIState,
          TISSUE_DESCENDANTS
        );

        expect(selectedCategoryValueIds.length).toEqual(0);
      });

      /**
       * Currently selected - immune system, hematopoietic system, blood, bone marrow, spleen, thymus
       * Remove - hematopoietic system
       * Still selected - immune system, bone marrow, spleen, thymus
       */
      test("removes system and all selected descendants of system that are not children of another selected system", () => {
        const categoryFilterUIState = {
          ...BASE_CATEGORY_FILTER_UI_STATE,
          selected: [
            INFERRED_IMMUNE_SYSTEM,
            INFERRED_HEMATOPOIETIC_SYSTEM,
            INFERRED_BLOOD,
            INFERRED_BONE_MARROW,
            INFERRED_SPLEEN,
            INFERRED_THYMUS,
          ],
          selectedPartial: [INFERRED_IMMUNE_SYSTEM],
        };

        const selectedCategoryValueIds = onRemoveMultiPanelCategoryValueTag(
          INFERRED_HEMATOPOIETIC_SYSTEM,
          [
            INFERRED_IMMUNE_SYSTEM,
            INFERRED_BLOOD,
            INFERRED_BONE_MARROW,
            INFERRED_SPLEEN,
            INFERRED_THYMUS,
          ],
          categoryFilterUIState,
          TISSUE_DESCENDANTS
        );

        expect(selectedCategoryValueIds.length).toEqual(4);
        expect(
          selectedCategoryValueIds.includes(INFERRED_IMMUNE_SYSTEM)
        ).toBeTruthy();
        expect(
          selectedCategoryValueIds.includes(INFERRED_BONE_MARROW)
        ).toBeTruthy();
        expect(selectedCategoryValueIds.includes(INFERRED_SPLEEN)).toBeTruthy();
        expect(selectedCategoryValueIds.includes(INFERRED_THYMUS)).toBeTruthy();
      });

      /**
       * Currently selected - blood, blood non-specific
       * Remove - blood non-specific
       * Still selected - blood
       */
      test("unravels to partially selected organ on remove of tissue", () => {
        const categoryFilterUIState = {
          ...BASE_CATEGORY_FILTER_UI_STATE,
          selected: [INFERRED_BLOOD, EXPLICIT_BLOOD],
          selectedPartial: [],
        };

        const selectedCategoryValueIds = onRemoveMultiPanelCategoryValueTag(
          EXPLICIT_BLOOD,
          [INFERRED_BLOOD],
          categoryFilterUIState,
          TISSUE_DESCENDANTS
        );

        expect(selectedCategoryValueIds.length).toEqual(1);
        expect(selectedCategoryValueIds.includes(INFERRED_BLOOD)).toBeTruthy();
      });

      /**
       * Currently selected - hematopoietic system, blood
       * Remove - blood
       * Still selected - hematopoietic system
       */
      test("unravels to partially selected system on remove of organ", () => {
        const categoryFilterUIState = {
          ...BASE_CATEGORY_FILTER_UI_STATE,
          selected: [INFERRED_HEMATOPOIETIC_SYSTEM, INFERRED_BLOOD],
          selectedPartial: [],
        };

        const selectedCategoryValueIds = onRemoveMultiPanelCategoryValueTag(
          INFERRED_BLOOD,
          [INFERRED_HEMATOPOIETIC_SYSTEM],
          categoryFilterUIState,
          TISSUE_DESCENDANTS
        );

        expect(selectedCategoryValueIds.length).toEqual(1);
        expect(
          selectedCategoryValueIds.includes(INFERRED_HEMATOPOIETIC_SYSTEM)
        ).toBeTruthy();
      });

      /**
       * Currently selected - hematopoietic system, blood non-specifc
       * Remove - blood non-specifc
       * Still selected - hematopoietic system
       */
      test("unravels to partially selected system on remove of tissue", () => {
        const categoryFilterUIState = {
          ...BASE_CATEGORY_FILTER_UI_STATE,
          selected: [INFERRED_HEMATOPOIETIC_SYSTEM, EXPLICIT_BLOOD],
          selectedPartial: [],
        };

        const selectedCategoryValueIds = onRemoveMultiPanelCategoryValueTag(
          EXPLICIT_BLOOD,
          [INFERRED_HEMATOPOIETIC_SYSTEM],
          categoryFilterUIState,
          TISSUE_DESCENDANTS
        );

        expect(selectedCategoryValueIds.length).toEqual(1);
        expect(
          selectedCategoryValueIds.includes(INFERRED_HEMATOPOIETIC_SYSTEM)
        ).toBeTruthy();
      });
    });

    describe("listPartiallySelectedCategoryValues", () => {
      /**
       * Selected - blood non-specific
       * Selected partial - none
       */
      test("doesn't lists tissue as partially selected if only tissue is selected", () => {
        const selectedPartial = listPartiallySelectedCategoryValues(
          [EXPLICIT_BLOOD],
          UI_NODES_BY_CATEGORY_VALUE_ID
        );

        expect(selectedPartial.length).toEqual(0);
      });

      /**
       * Selected - blood, blood non-specific
       * Selected partial - blood
       */
      test("lists organ as partially selected when tissue is selected", () => {
        const selectedPartial = listPartiallySelectedCategoryValues(
          [INFERRED_BLOOD, EXPLICIT_BLOOD],
          UI_NODES_BY_CATEGORY_VALUE_ID
        );

        expect(selectedPartial.length).toEqual(1);
        expect(selectedPartial?.[0]).toEqual(INFERRED_BLOOD);
      });

      /**
       * Selected - blood, blood non-specific, umbilical cord blood, venous blood
       * Selected partial - none
       */
      test("doesn't lists organ as partially selected if all children tissue are also selected", () => {
        const selectedPartial = listPartiallySelectedCategoryValues(
          [
            INFERRED_BLOOD,
            EXPLICIT_BLOOD,
            EXPLICIT_UMBILICAL_CORD_BLOOD,
            EXPLICIT_VENOUS_BLOOD,
          ],
          UI_NODES_BY_CATEGORY_VALUE_ID
        );

        expect(selectedPartial.length).toEqual(0);
      });

      /**
       * Selected - spleen, spleen non-specific
       * Selected partial - none
       */
      test("doesn't list organ as partially selected when single child tissue is selected", () => {
        const selectedPartial = listPartiallySelectedCategoryValues(
          [INFERRED_SPLEEN, EXPLICIT_SPLEEN],
          UI_NODES_BY_CATEGORY_VALUE_ID
        );

        expect(selectedPartial.length).toEqual(0);
      });

      /**
       * Selected - hematopoietic system, blood
       * Selected partial - hematopoietic system
       */
      test("lists system as partially selected when not all children organs are selected", () => {
        const selectedPartial = listPartiallySelectedCategoryValues(
          [INFERRED_HEMATOPOIETIC_SYSTEM, INFERRED_BLOOD],
          UI_NODES_BY_CATEGORY_VALUE_ID
        );

        expect(selectedPartial.length).toEqual(1);
        expect(selectedPartial?.[0]).toEqual(INFERRED_HEMATOPOIETIC_SYSTEM);
      });

      /**
       * Selected - hematopoietic system, blood, bone marrow, spleen, thymus
       * Selected partial - none
       */
      test("doesn't list system as partially selected if all children organs are also selected", () => {
        const selectedPartial = listPartiallySelectedCategoryValues(
          [
            INFERRED_HEMATOPOIETIC_SYSTEM,
            INFERRED_BLOOD,
            INFERRED_BONE_MARROW,
            INFERRED_SPLEEN,
            INFERRED_THYMUS,
          ],
          UI_NODES_BY_CATEGORY_VALUE_ID
        );

        expect(selectedPartial.length).toEqual(0);
      });

      /**
       * Selected - hematopoietic system, bone marrow, spleen, thymus, blood non-specific, umbilical cord blood,
       * venous blood.
       * Selected partial - none
       */
      test("doesn't lists system as partially selected if all children organs or children tissues are selected", () => {
        const selectedPartial = listPartiallySelectedCategoryValues(
          [
            INFERRED_HEMATOPOIETIC_SYSTEM,
            INFERRED_BONE_MARROW,
            INFERRED_SPLEEN,
            INFERRED_THYMUS,
            EXPLICIT_BLOOD,
            EXPLICIT_UMBILICAL_CORD_BLOOD,
            EXPLICIT_VENOUS_BLOOD,
          ],
          UI_NODES_BY_CATEGORY_VALUE_ID
        );

        expect(selectedPartial.length).toEqual(0);
      });

      /**
       * Selected - hematopoietic system, blood, bone marrow, spleen, thymus, blood non-specific
       * Selected partial - hematopoietic system, blood
       */
      test("lists system as partially selected if not all children organs or children tissues are selected", () => {
        const selectedPartial = listPartiallySelectedCategoryValues(
          [
            INFERRED_HEMATOPOIETIC_SYSTEM,
            INFERRED_BLOOD,
            INFERRED_BONE_MARROW,
            INFERRED_SPLEEN,
            INFERRED_THYMUS,
            EXPLICIT_BLOOD,
          ],
          UI_NODES_BY_CATEGORY_VALUE_ID
        );

        expect(selectedPartial.length).toEqual(2);
        expect(
          selectedPartial.includes(INFERRED_HEMATOPOIETIC_SYSTEM)
        ).toBeTruthy();
        expect(selectedPartial.includes(INFERRED_BLOOD)).toBeTruthy();
      });
    });

    describe("listMultiPanelSelectedViews", () => {
      /**
       * Selected - blood non-specific
       * Selected partial - none
       * Selected values - blood non-specific
       */
      test("builds selected views for tissue", () => {
        const categoryFilterUIState = {
          ...BASE_CATEGORY_FILTER_UI_STATE,
          selected: [EXPLICIT_BLOOD],
          selectedPartial: [],
        };

        const valueViews = new Map(VALUE_VIEWS_BY_CATEGORY_VALUE_ID);
        updateCategoryValueViewSelected(
          valueViews,
          EXPLICIT_BLOOD,
          true,
          false
        );

        const selectedViews = listMultiPanelSelectedViews(
          [...valueViews.values()],
          categoryFilterUIState
        );

        expect(selectedViews.length).toEqual(1);
        expect(selectedViews?.[0]?.categoryValueId).toEqual(EXPLICIT_BLOOD);
      });

      /**
       * Selected - spleen, spleen non-specific
       * Selected partial - none
       * Selected values - spleen
       */
      test("builds selected views for selected organ with single selected tissue", () => {
        const categoryFilterUIState = {
          ...BASE_CATEGORY_FILTER_UI_STATE,
          selected: [INFERRED_SPLEEN, EXPLICIT_SPLEEN],
          selectedPartial: [],
        };

        const valueViews = new Map(VALUE_VIEWS_BY_CATEGORY_VALUE_ID);
        updateCategoryValueViewSelected(
          valueViews,
          INFERRED_SPLEEN,
          true,
          false
        );
        updateCategoryValueViewSelected(
          valueViews,
          EXPLICIT_SPLEEN,
          true,
          false
        );

        const selectedViews = listMultiPanelSelectedViews(
          [...valueViews.values()],
          categoryFilterUIState
        );

        expect(selectedViews.length).toEqual(1);
        expect(selectedViews?.[0]?.categoryValueId).toEqual(INFERRED_SPLEEN);
      });

      /**
       * Selected - blood non-specific, umbilical cord blood, venous blood
       * Selected partial - none
       * Selected values - blood non-specific, umbilical cord blood, venous blood
       */
      test("builds selected views for all children of blood selected", () => {
        const categoryFilterUIState = {
          ...BASE_CATEGORY_FILTER_UI_STATE,
          selected: [
            EXPLICIT_BLOOD,
            EXPLICIT_UMBILICAL_CORD_BLOOD,
            EXPLICIT_VENOUS_BLOOD,
          ],
          selectedPartial: [],
        };

        const valueViews = new Map(VALUE_VIEWS_BY_CATEGORY_VALUE_ID);
        updateCategoryValueViewSelected(
          valueViews,
          EXPLICIT_BLOOD,
          true,
          false
        );
        updateCategoryValueViewSelected(
          valueViews,
          EXPLICIT_UMBILICAL_CORD_BLOOD,
          true,
          false
        );
        updateCategoryValueViewSelected(
          valueViews,
          EXPLICIT_VENOUS_BLOOD,
          true,
          false
        );

        const selectedViews = listMultiPanelSelectedViews(
          [...valueViews.values()],
          categoryFilterUIState
        );

        expect(selectedViews.length).toEqual(3);
        const selectedKeys = listCategoryValueIds(selectedViews);
        expect(selectedKeys.includes(EXPLICIT_BLOOD)).toBeTruthy();
        expect(
          selectedKeys.includes(EXPLICIT_UMBILICAL_CORD_BLOOD)
        ).toBeTruthy();
        expect(selectedKeys.includes(EXPLICIT_VENOUS_BLOOD)).toBeTruthy();
      });

      /**
       * Selected - blood, blood non-specific
       * Selected partial - blood
       * Selected values - blood non-specific
       */
      test("builds selected views for organ and corresponding specific tissue selected", () => {
        const categoryFilterUIState = {
          ...BASE_CATEGORY_FILTER_UI_STATE,
          selected: [INFERRED_BLOOD, EXPLICIT_BLOOD],
          selectedPartial: [INFERRED_BLOOD],
        };

        const valueViews = new Map(VALUE_VIEWS_BY_CATEGORY_VALUE_ID);
        updateCategoryValueViewSelected(
          valueViews,
          EXPLICIT_BLOOD,
          true,
          false
        );
        updateCategoryValueViewSelected(
          valueViews,
          INFERRED_BLOOD,
          false,
          true
        );

        const selectedViews = listMultiPanelSelectedViews(
          [...valueViews.values()],
          categoryFilterUIState
        );

        expect(selectedViews.length).toEqual(1);
        expect(selectedViews?.[0]?.categoryValueId).toEqual(EXPLICIT_BLOOD);
      });

      /**
       * Selected - blood, blood non-specific, umbilical cord blood, venous blood
       * Selected partial - none
       * Selected values - blood
       */
      test("rolls up selected organ and all selected organ's tissues to just organ", () => {
        const categoryFilterUIState = {
          ...BASE_CATEGORY_FILTER_UI_STATE,
          selected: [
            INFERRED_BLOOD,
            EXPLICIT_BLOOD,
            EXPLICIT_UMBILICAL_CORD_BLOOD,
            EXPLICIT_VENOUS_BLOOD,
          ],
          selectedPartial: [],
        };

        const valueViews = new Map(VALUE_VIEWS_BY_CATEGORY_VALUE_ID);
        updateCategoryValueViewSelected(
          valueViews,
          INFERRED_BLOOD,
          true,
          false
        );
        updateCategoryValueViewSelected(
          valueViews,
          EXPLICIT_BLOOD,
          true,
          false
        );
        updateCategoryValueViewSelected(
          valueViews,
          EXPLICIT_UMBILICAL_CORD_BLOOD,
          true,
          false
        );
        updateCategoryValueViewSelected(
          valueViews,
          EXPLICIT_VENOUS_BLOOD,
          true,
          false
        );

        const selectedViews = listMultiPanelSelectedViews(
          [...valueViews.values()],
          categoryFilterUIState
        );

        expect(selectedViews.length).toEqual(1);
        expect(selectedViews?.[0]?.categoryValueId).toEqual(INFERRED_BLOOD);
      });

      /**
       * Selected - hematopoietic system, blood, bone marrow, spleen, thymus
       * Selected partial - none
       * Selected values - hematopoietic system
       */
      test("rolls up selected system and all selected system's organs to just system", () => {
        const categoryFilterUIState = {
          ...BASE_CATEGORY_FILTER_UI_STATE,
          selected: [
            INFERRED_HEMATOPOIETIC_SYSTEM,
            INFERRED_BLOOD,
            INFERRED_BONE_MARROW,
            INFERRED_SPLEEN,
            INFERRED_THYMUS,
          ],
          selectedPartial: [],
        };

        const valueViews = new Map(VALUE_VIEWS_BY_CATEGORY_VALUE_ID);
        updateCategoryValueViewSelected(
          valueViews,
          INFERRED_HEMATOPOIETIC_SYSTEM,
          true,
          false
        );
        updateCategoryValueViewSelected(
          valueViews,
          INFERRED_BLOOD,
          true,
          false
        );
        updateCategoryValueViewSelected(
          valueViews,
          INFERRED_BONE_MARROW,
          true,
          false
        );
        updateCategoryValueViewSelected(
          valueViews,
          INFERRED_SPLEEN,
          true,
          false
        );
        updateCategoryValueViewSelected(
          valueViews,
          INFERRED_THYMUS,
          true,
          false
        );

        const selectedViews = listMultiPanelSelectedViews(
          [...valueViews.values()],
          categoryFilterUIState
        );

        expect(selectedViews.length).toEqual(1);
        expect(selectedViews?.[0]?.categoryValueId).toEqual(
          INFERRED_HEMATOPOIETIC_SYSTEM
        );
      });

      /**
       * Selected - hematopoietic system, blood, bone marrow, spleen, thymus
       * Selected partial - immune system
       * Selected values - hematopoietic system, bone marrow, spleen, thymus
       */
      test("builds selected views for one selected system and one partially selected system, with shared selected organs", () => {
        const categoryFilterUIState = {
          ...BASE_CATEGORY_FILTER_UI_STATE,
          selected: [
            INFERRED_HEMATOPOIETIC_SYSTEM,
            INFERRED_IMMUNE_SYSTEM,
            INFERRED_BLOOD,
            INFERRED_BONE_MARROW,
            INFERRED_SPLEEN,
            INFERRED_THYMUS,
          ],
          selectedPartial: [INFERRED_IMMUNE_SYSTEM],
        };

        const valueViews = new Map(VALUE_VIEWS_BY_CATEGORY_VALUE_ID);
        updateCategoryValueViewSelected(
          valueViews,
          INFERRED_HEMATOPOIETIC_SYSTEM,
          true,
          false
        );
        updateCategoryValueViewSelected(
          valueViews,
          INFERRED_IMMUNE_SYSTEM,
          false,
          true
        );
        updateCategoryValueViewSelected(
          valueViews,
          INFERRED_BLOOD,
          true,
          false
        );
        updateCategoryValueViewSelected(
          valueViews,
          INFERRED_BONE_MARROW,
          true,
          false
        );
        updateCategoryValueViewSelected(
          valueViews,
          INFERRED_SPLEEN,
          true,
          false
        );
        updateCategoryValueViewSelected(
          valueViews,
          INFERRED_THYMUS,
          true,
          false
        );

        const selectedViews = listMultiPanelSelectedViews(
          [...valueViews.values()],
          categoryFilterUIState
        );

        expect(selectedViews.length).toEqual(4);

        const selectedKeys = listCategoryValueIds(selectedViews);
        expect(
          selectedKeys.includes(INFERRED_HEMATOPOIETIC_SYSTEM)
        ).toBeTruthy();
        expect(selectedKeys.includes(INFERRED_BONE_MARROW)).toBeTruthy();
        expect(selectedKeys.includes(INFERRED_SPLEEN)).toBeTruthy();
        expect(selectedKeys.includes(INFERRED_THYMUS)).toBeTruthy();
      });

      /**
       * Selected - hematopoietic system, bone marrow, spleen, thymus, blood non-specific, umbilical cord blood, venous blood
       * Selected partial - none
       * Selected values - hematopoietic system
       */
      test("builds selected views for selected system and selected organs except those for a selected tissue", () => {
        const categoryFilterUIState = {
          ...BASE_CATEGORY_FILTER_UI_STATE,
          selected: [
            INFERRED_HEMATOPOIETIC_SYSTEM,
            INFERRED_BONE_MARROW,
            INFERRED_SPLEEN,
            INFERRED_THYMUS,
            EXPLICIT_BLOOD,
            EXPLICIT_UMBILICAL_CORD_BLOOD,
            EXPLICIT_VENOUS_BLOOD,
          ],
          selectedPartial: [],
        };

        const valueViews = new Map(VALUE_VIEWS_BY_CATEGORY_VALUE_ID);
        updateCategoryValueViewSelected(
          valueViews,
          INFERRED_HEMATOPOIETIC_SYSTEM,
          true,
          false
        );
        updateCategoryValueViewSelected(
          valueViews,
          INFERRED_BONE_MARROW,
          true,
          false
        );
        updateCategoryValueViewSelected(
          valueViews,
          INFERRED_SPLEEN,
          true,
          false
        );
        updateCategoryValueViewSelected(
          valueViews,
          INFERRED_THYMUS,
          true,
          false
        );
        updateCategoryValueViewSelected(
          valueViews,
          EXPLICIT_BLOOD,
          true,
          false
        );
        updateCategoryValueViewSelected(
          valueViews,
          EXPLICIT_UMBILICAL_CORD_BLOOD,
          true,
          false
        );
        updateCategoryValueViewSelected(
          valueViews,
          EXPLICIT_VENOUS_BLOOD,
          true,
          false
        );

        const selectedViews = listMultiPanelSelectedViews(
          [...valueViews.values()],
          categoryFilterUIState
        );

        expect(selectedViews.length).toEqual(1);
        expect(selectedViews?.[0]?.categoryValueId).toEqual(
          INFERRED_HEMATOPOIETIC_SYSTEM
        );
      });

      /**
       * Selected - bone marrow, spleen, thymus, blood non-specific
       * Selected partial - hematopoietic system, blood
       * Selected values - bone marrow, spleen, thymus, blood non-specific
       */
      test("builds selected views for selected system, partially selected organs, selected tissue", () => {
        const categoryFilterUIState = {
          ...BASE_CATEGORY_FILTER_UI_STATE,
          selected: [
            INFERRED_BONE_MARROW,
            INFERRED_SPLEEN,
            INFERRED_THYMUS,
            EXPLICIT_BLOOD,
          ],
          selectedPartial: [INFERRED_HEMATOPOIETIC_SYSTEM, INFERRED_BLOOD],
        };

        const valueViews = new Map(VALUE_VIEWS_BY_CATEGORY_VALUE_ID);
        updateCategoryValueViewSelected(
          valueViews,
          INFERRED_HEMATOPOIETIC_SYSTEM,
          false,
          true
        );
        updateCategoryValueViewSelected(
          valueViews,
          INFERRED_BLOOD,
          false,
          true
        );
        updateCategoryValueViewSelected(
          valueViews,
          INFERRED_BONE_MARROW,
          true,
          false
        );
        updateCategoryValueViewSelected(
          valueViews,
          INFERRED_SPLEEN,
          true,
          false
        );
        updateCategoryValueViewSelected(
          valueViews,
          INFERRED_THYMUS,
          true,
          false
        );
        updateCategoryValueViewSelected(
          valueViews,
          EXPLICIT_BLOOD,
          true,
          false
        );

        const selectedViews = listMultiPanelSelectedViews(
          [...valueViews.values()],
          categoryFilterUIState
        );

        expect(selectedViews.length).toEqual(4);

        const selectedKeys = listCategoryValueIds(selectedViews);
        expect(selectedKeys.includes(INFERRED_BONE_MARROW)).toBeTruthy();
        expect(selectedKeys.includes(INFERRED_SPLEEN)).toBeTruthy();
        expect(selectedKeys.includes(INFERRED_THYMUS)).toBeTruthy();
        expect(selectedKeys.includes(EXPLICIT_BLOOD)).toBeTruthy();
      });
    });

    describe("buildUINodesByCategoryValueId", () => {
      let uiNodesByCategoryValueId: Map<CategoryValueId, MultiPanelUINode>;
      test.beforeAll(() => {
        uiNodesByCategoryValueId = buildUINodesByCategoryValueId(
          TISSUE_CATEGORY_VALUE_IDS_BY_PANEL,
          TISSUE_DESCENDANTS
        );
      });

      /**
       * Build parent child relationships for hematopoietic system. hematopoietic system has "UI" children in both organ
       * and tissue panels.
       * - Parents: -
       * - Children: blood non-specific, umbilical cord blood, venous blood
       */
      test("builds parent child relationships for hematopoietic system", () => {
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
        expect(uiChildren?.includes(INFERRED_THYMUS)).toBeTruthy();
      });

      /**
       * Build parent child relationships for renal system. renal system only has organs as direct "UI" children.
       * - Parents: -
       * - Children: kidney, urethra, bladder lumen, ureter
       */
      test("builds parent child relationships for renal system", () => {
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
      test("builds parent child relationships for blood", () => {
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
      test("builds parent child relationships for bladder lumen", () => {
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
      test("builds parent child relationships for ureter", () => {
        const uiNode = uiNodesByCategoryValueId?.get(EXPLICIT_URETER);

        const uiParents = uiNode?.uiParents;
        expect(uiParents?.length).toEqual(1);
        expect(uiParents?.includes(INFERRED_RENAL_SYSTEM)).toBeTruthy();

        const uiChildren = uiNode?.uiChildren;
        expect(uiChildren?.length).toEqual(0);
      });
    });

    describe("buildMultiPanelCategoryView", () => {
      let multiPanelUIState: MultiPanelUIState;
      test.beforeAll(() => {
        const uiNodesByCategoryValueId = buildUINodesByCategoryValueId(
          TISSUE_CATEGORY_VALUE_IDS_BY_PANEL,
          TISSUE_DESCENDANTS
        );

        const categoryFilterUISTate = {
          selected: [],
          selectedPartial: [],
          uiNodesByCategoryValueId,
        };

        multiPanelUIState = new Map<
          CATEGORY_FILTER_ID,
          MultiPanelCategoryFilterUIState
        >([[CATEGORY_FILTER_ID.TISSUE_CALCULATED, categoryFilterUISTate]]);
      });

      /**
       * Selected: none
       */
      describe("nothing selected", () => {
        let categoryView: MultiPanelOntologyCategoryView;
        test.beforeAll(() => {
          categoryView = buildMultiPanelCategoryView(
            CATEGORY_FILTER_CONFIGS_BY_ID[
              CATEGORY_FILTER_ID.TISSUE_CALCULATED
            ] as MultiPanelOntologyFilterConfig,
            KEYED_CATEGORY_VALUES,
            multiPanelUIState,
            new Map<string, string>()
          );
        });

        /**
         * Views: all
         */
        test("includes all views", () => {
          expect(categoryView.panels.length).toEqual(3);

          const systemPanelViews =
            // eslint-disable-next-line @typescript-eslint/no-non-null-assertion -- length check above
            categoryView.panels[PANEL_INDEX_SYSTEM]!.views;
          expect(systemPanelViews.length).toEqual(3);
          const systemCategoryValueIds = listCategoryValueIds(systemPanelViews);
          TISSUE_SYSTEMS.forEach((categoryValueId) =>
            systemCategoryValueIds.includes(categoryValueId)
          );

          // eslint-disable-next-line @typescript-eslint/no-non-null-assertion -- length check above
          const organPanelViews = categoryView.panels[PANEL_INDEX_ORGAN]!.views;
          expect(organPanelViews.length).toEqual(7);
          const organCategoryValueIds = listCategoryValueIds(organPanelViews);
          TISSUE_ORGANS.forEach((categoryValueId) =>
            organCategoryValueIds.includes(categoryValueId)
          );

          const tissuePanelViews =
            // eslint-disable-next-line @typescript-eslint/no-non-null-assertion -- length check above
            categoryView.panels[PANEL_INDEX_TISSUE]!.views;
          expect(tissuePanelViews.length).toEqual(11);
          const tissueCategoryValueIds = listCategoryValueIds(tissuePanelViews);
          TISSUES.forEach((categoryValueId) =>
            tissueCategoryValueIds.includes(categoryValueId)
          );
        });

        /**
         * Selected - none
         */
        test("includes no selected values", () => {
          categoryView.panels.forEach((panel) => {
            panel.views.forEach((view) => {
              expect(view.selected).toBeFalsy();
              expect(view.selectedPartial).toBeFalsy();
            });
          });
        });
      });

      /**
       * Selected: blood non-specific
       */
      describe("tissue selected", () => {
        let categoryView: MultiPanelOntologyCategoryView;
        test.beforeAll(() => {
          const updatedKeyedCategoryValues = new Map(KEYED_CATEGORY_VALUES);
          updateSelectCategoryValueSelected(
            updatedKeyedCategoryValues,
            EXPLICIT_BLOOD,
            true,
            false
          );

          const updatedMultiPanelUIState = new Map(multiPanelUIState);
          updateCategoryFilterUIStateSelected(
            updatedMultiPanelUIState,
            CATEGORY_FILTER_ID.TISSUE_CALCULATED,
            [EXPLICIT_BLOOD],
            []
          );

          categoryView = buildMultiPanelCategoryView(
            CATEGORY_FILTER_CONFIGS_BY_ID[
              CATEGORY_FILTER_ID.TISSUE_CALCULATED
            ] as MultiPanelOntologyFilterConfig,
            updatedKeyedCategoryValues,
            updatedMultiPanelUIState,
            new Map<string, string>()
          );
        });

        /**
         * Views: all systems, all organs, all tissues
         */
        test("includes all systems, organs and tissue views", () => {
          expect(categoryView.panels.length).toEqual(3);

          const systemPanelViews =
            // eslint-disable-next-line @typescript-eslint/no-non-null-assertion -- length check above
            categoryView.panels[PANEL_INDEX_SYSTEM]!.views;
          expect(systemPanelViews.length).toEqual(3);
          const systemCategoryValueIds = listCategoryValueIds(systemPanelViews);
          TISSUE_SYSTEMS.forEach((categoryValueId) =>
            systemCategoryValueIds.includes(categoryValueId)
          );

          // eslint-disable-next-line @typescript-eslint/no-non-null-assertion -- length check above
          const organPanelViews = categoryView.panels[PANEL_INDEX_ORGAN]!.views;
          expect(organPanelViews.length).toEqual(7);
          const organCategoryValueIds = listCategoryValueIds(organPanelViews);
          TISSUE_ORGANS.forEach((categoryValueId) =>
            organCategoryValueIds.includes(categoryValueId)
          );

          const tissuePanelViews =
            // eslint-disable-next-line @typescript-eslint/no-non-null-assertion -- length check above
            categoryView.panels[PANEL_INDEX_TISSUE]!.views;
          expect(tissuePanelViews.length).toEqual(11);
          const tissueCategoryValueIds = listCategoryValueIds(tissuePanelViews);
          TISSUES.forEach((categoryValueId) =>
            tissueCategoryValueIds.includes(categoryValueId)
          );
        });

        /**
         * Selected - blood non-specific
         */
        test("includes blood non-specific selected only", () => {
          categoryView.panels.forEach((panel) => {
            panel.views.forEach((view) => {
              const selected = view.categoryValueId === EXPLICIT_BLOOD;
              expect(view.selected).toEqual(selected);
              expect(view.selectedPartial).toBeFalsy();
            });
          });
        });
      });

      /**
       * Selected: blood
       */
      describe("organ selected", () => {
        let categoryView: MultiPanelOntologyCategoryView;
        test.beforeAll(() => {
          const updatedKeyedCategoryValues = new Map(KEYED_CATEGORY_VALUES);
          updateSelectCategoryValueSelected(
            updatedKeyedCategoryValues,
            INFERRED_BLOOD,
            true,
            false
          );

          const updatedMultiPanelUIState = new Map(multiPanelUIState);
          updateCategoryFilterUIStateSelected(
            updatedMultiPanelUIState,
            CATEGORY_FILTER_ID.TISSUE_CALCULATED,
            [INFERRED_BLOOD],
            []
          );

          categoryView = buildMultiPanelCategoryView(
            CATEGORY_FILTER_CONFIGS_BY_ID[
              CATEGORY_FILTER_ID.TISSUE_CALCULATED
            ] as MultiPanelOntologyFilterConfig,
            updatedKeyedCategoryValues,
            updatedMultiPanelUIState,
            new Map<string, string>()
          );
        });

        /**
         * Views: all systems, all organs, blood non-specific, umbilical cord blood, venous blood
         */
        test("includes all systems, all organs and filtered tissues", () => {
          expect(categoryView.panels.length).toEqual(3);

          const systemPanelViews =
            // eslint-disable-next-line @typescript-eslint/no-non-null-assertion -- length check above
            categoryView.panels[PANEL_INDEX_SYSTEM]!.views;
          expect(systemPanelViews.length).toEqual(3);
          const systemCategoryValueIds = listCategoryValueIds(systemPanelViews);
          TISSUE_SYSTEMS.forEach((categoryValueId) =>
            systemCategoryValueIds.includes(categoryValueId)
          );

          // eslint-disable-next-line @typescript-eslint/no-non-null-assertion -- length check above
          const organPanelViews = categoryView.panels[PANEL_INDEX_ORGAN]!.views;
          expect(organPanelViews.length).toEqual(7);
          const organCategoryValueIds = listCategoryValueIds(organPanelViews);
          TISSUE_ORGANS.forEach((categoryValueId) =>
            organCategoryValueIds.includes(categoryValueId)
          );

          const tissuePanelViews =
            // eslint-disable-next-line @typescript-eslint/no-non-null-assertion -- length check above
            categoryView.panels[PANEL_INDEX_TISSUE]!.views;
          expect(
            tissuePanelViews.filter((view) => view.visible).length
          ).toEqual(3);
          const tissueCategoryValueIds = listCategoryValueIds(tissuePanelViews);
          [
            EXPLICIT_BLOOD,
            EXPLICIT_UMBILICAL_CORD_BLOOD,
            EXPLICIT_VENOUS_BLOOD,
          ].forEach((categoryValueId) =>
            tissueCategoryValueIds.includes(categoryValueId)
          );
        });

        /**
         * Selected - blood
         */
        test("includes blood selected only", () => {
          categoryView.panels.forEach((panel) => {
            panel.views.forEach((view) => {
              const selected = view.categoryValueId === INFERRED_BLOOD;
              expect(view.selected).toEqual(selected);
              expect(view.selectedPartial).toBeFalsy();
            });
          });
        });
      });

      describe("system selected", () => {
        let categoryView: MultiPanelOntologyCategoryView;
        test.beforeAll(() => {
          const updatedKeyedCategoryValues = new Map(KEYED_CATEGORY_VALUES);
          updateSelectCategoryValueSelected(
            updatedKeyedCategoryValues,
            INFERRED_IMMUNE_SYSTEM,
            true,
            false
          );

          const updatedMultiPanelUIState = new Map(multiPanelUIState);
          updateCategoryFilterUIStateSelected(
            updatedMultiPanelUIState,
            CATEGORY_FILTER_ID.TISSUE_CALCULATED,
            [INFERRED_IMMUNE_SYSTEM],
            []
          );

          categoryView = buildMultiPanelCategoryView(
            CATEGORY_FILTER_CONFIGS_BY_ID[
              CATEGORY_FILTER_ID.TISSUE_CALCULATED
            ] as MultiPanelOntologyFilterConfig,
            updatedKeyedCategoryValues,
            updatedMultiPanelUIState,
            new Map<string, string>()
          );
        });

        /**
         * Selected: immune system
         * Views: all systems, organs under immune system, tissues under immune system
         */
        test("includes all systems and filtered organs and tissues for immune system", () => {
          expect(categoryView.panels.length).toEqual(3);

          const systemPanelViews =
            // eslint-disable-next-line @typescript-eslint/no-non-null-assertion -- length check above
            categoryView.panels[PANEL_INDEX_SYSTEM]!.views;
          expect(systemPanelViews.length).toEqual(3);
          const systemCategoryValueIds = listCategoryValueIds(systemPanelViews);
          TISSUE_SYSTEMS.forEach((categoryValueId) =>
            systemCategoryValueIds.includes(categoryValueId)
          );

          // eslint-disable-next-line @typescript-eslint/no-non-null-assertion -- length check above
          const organPanelViews = categoryView.panels[PANEL_INDEX_ORGAN]!.views;
          expect(organPanelViews.filter((view) => view.visible).length).toEqual(
            4
          );
          const organCategoryValueIds = listCategoryValueIds(organPanelViews);
          [
            INFERRED_BONE_MARROW,
            INFERRED_LYMPH_NODE,
            INFERRED_SPLEEN,
            INFERRED_THYMUS,
          ].forEach((categoryValueId) =>
            organCategoryValueIds.includes(categoryValueId)
          );

          const tissuePanelViews =
            // eslint-disable-next-line @typescript-eslint/no-non-null-assertion -- length check above
            categoryView.panels[PANEL_INDEX_TISSUE]!.views;
          expect(
            tissuePanelViews.filter((view) => view.visible).length
          ).toEqual(3);
          const tissueCategoryValueIds = listCategoryValueIds(tissuePanelViews);
          [EXPLICIT_BONE_MARROW, EXPLICIT_SPLEEN, EXPLICIT_THYMUS].forEach(
            (categoryValueId) =>
              tissueCategoryValueIds.includes(categoryValueId)
          );
        });

        /**
         * Selected - immune system
         */
        test("includes immune system selected only", () => {
          categoryView.panels.forEach((panel) => {
            panel.views.forEach((view) => {
              const selected = view.categoryValueId === INFERRED_IMMUNE_SYSTEM;
              expect(view.selected).toEqual(selected);
              expect(view.selectedPartial).toBeFalsy();
            });
          });
        });
      });

      /**
       * Selected: blood, blood non-specific
       */
      describe("organ and tissue selected", () => {
        let categoryView: MultiPanelOntologyCategoryView;
        test.beforeAll(() => {
          const updatedKeyedCategoryValues = new Map(KEYED_CATEGORY_VALUES);
          updateSelectCategoryValueSelected(
            updatedKeyedCategoryValues,
            INFERRED_BLOOD,
            false,
            true
          );
          updateSelectCategoryValueSelected(
            updatedKeyedCategoryValues,
            EXPLICIT_BLOOD,
            true,
            false
          );

          const updatedMultiPanelUIState = new Map(multiPanelUIState);
          updateCategoryFilterUIStateSelected(
            updatedMultiPanelUIState,
            CATEGORY_FILTER_ID.TISSUE_CALCULATED,
            [INFERRED_BLOOD, EXPLICIT_BLOOD],
            [INFERRED_BLOOD]
          );

          categoryView = buildMultiPanelCategoryView(
            CATEGORY_FILTER_CONFIGS_BY_ID[
              CATEGORY_FILTER_ID.TISSUE_CALCULATED
            ] as MultiPanelOntologyFilterConfig,
            updatedKeyedCategoryValues,
            updatedMultiPanelUIState,
            new Map<string, string>()
          );
        });

        /**
         * Views: all systems, all organs, tissues under blood
         */
        test("includes all systems and organs and filtered tissues", () => {
          expect(categoryView.panels.length).toEqual(3);

          const systemPanelViews =
            // eslint-disable-next-line @typescript-eslint/no-non-null-assertion -- length check above
            categoryView.panels[PANEL_INDEX_SYSTEM]!.views;
          expect(systemPanelViews.length).toEqual(3);
          const systemCategoryValueIds = listCategoryValueIds(systemPanelViews);
          TISSUE_SYSTEMS.forEach((categoryValueId) =>
            systemCategoryValueIds.includes(categoryValueId)
          );

          // eslint-disable-next-line @typescript-eslint/no-non-null-assertion -- length check above
          const organPanelViews = categoryView.panels[PANEL_INDEX_ORGAN]!.views;
          expect(organPanelViews.length).toEqual(7);
          const organCategoryValueIds = listCategoryValueIds(organPanelViews);
          TISSUE_ORGANS.forEach((categoryValueId) =>
            organCategoryValueIds.includes(categoryValueId)
          );

          const tissuePanelViews =
            // eslint-disable-next-line @typescript-eslint/no-non-null-assertion -- length check above
            categoryView.panels[PANEL_INDEX_TISSUE]!.views;
          expect(
            tissuePanelViews.filter((view) => view.visible).length
          ).toEqual(3);
          const tissueCategoryValueIds = listCategoryValueIds(tissuePanelViews);
          [
            EXPLICIT_BLOOD,
            EXPLICIT_UMBILICAL_CORD_BLOOD,
            EXPLICIT_VENOUS_BLOOD,
          ].forEach((categoryValueId) =>
            tissueCategoryValueIds.includes(categoryValueId)
          );
        });

        /**
         * Selected - blood
         * Selected partial - hematopoietic system
         */
        test("includes blood partially selected, blood non-specific selected", () => {
          categoryView.panels.forEach((panel) => {
            panel.views.forEach((view) => {
              const selected = view.categoryValueId === EXPLICIT_BLOOD;
              const selectedPartial = view.categoryValueId === INFERRED_BLOOD;
              expect(view.selected).toEqual(selected);
              expect(view.selectedPartial).toEqual(selectedPartial);
            });
          });
        });
      });

      /**
       * Selected: hematopoietic system, blood
       */
      describe("system and organ selected", () => {
        let categoryView: MultiPanelOntologyCategoryView;
        test.beforeAll(() => {
          const updatedKeyedCategoryValues = new Map(KEYED_CATEGORY_VALUES);
          updateSelectCategoryValueSelected(
            updatedKeyedCategoryValues,
            INFERRED_HEMATOPOIETIC_SYSTEM,
            false,
            true
          );
          updateSelectCategoryValueSelected(
            updatedKeyedCategoryValues,
            INFERRED_BLOOD,
            true,
            false
          );

          const updatedMultiPanelUIState = new Map(multiPanelUIState);
          updateCategoryFilterUIStateSelected(
            updatedMultiPanelUIState,
            CATEGORY_FILTER_ID.TISSUE_CALCULATED,
            [INFERRED_HEMATOPOIETIC_SYSTEM, INFERRED_BLOOD],
            [INFERRED_HEMATOPOIETIC_SYSTEM]
          );

          categoryView = buildMultiPanelCategoryView(
            CATEGORY_FILTER_CONFIGS_BY_ID[
              CATEGORY_FILTER_ID.TISSUE_CALCULATED
            ] as MultiPanelOntologyFilterConfig,
            updatedKeyedCategoryValues,
            updatedMultiPanelUIState,
            new Map<string, string>()
          );
        });

        /**
         * Views: all systems, organs under hematopoietic system, tissues under blood
         */
        test("includes all systems and filtered organs and tissues for hematopoietic system", () => {
          expect(categoryView.panels.length).toEqual(3);

          const systemPanelViews =
            // eslint-disable-next-line @typescript-eslint/no-non-null-assertion -- length check above
            categoryView.panels[PANEL_INDEX_SYSTEM]!.views;
          expect(systemPanelViews.length).toEqual(3);
          const systemCategoryValueIds = listCategoryValueIds(systemPanelViews);
          TISSUE_SYSTEMS.forEach((categoryValueId) =>
            systemCategoryValueIds.includes(categoryValueId)
          );

          // eslint-disable-next-line @typescript-eslint/no-non-null-assertion -- length check above
          const organPanelViews = categoryView.panels[PANEL_INDEX_ORGAN]!.views;
          expect(organPanelViews.filter((view) => view.visible).length).toEqual(
            4
          );
          const organCategoryValueIds = listCategoryValueIds(organPanelViews);
          [
            INFERRED_BLOOD,
            INFERRED_BONE_MARROW,
            INFERRED_SPLEEN,
            INFERRED_THYMUS,
          ].forEach((categoryValueId) =>
            organCategoryValueIds.includes(categoryValueId)
          );

          const tissuePanelViews =
            // eslint-disable-next-line @typescript-eslint/no-non-null-assertion -- length check above
            categoryView.panels[PANEL_INDEX_TISSUE]!.views;
          expect(
            tissuePanelViews.filter((view) => view.visible).length
          ).toEqual(3);
          const tissueCategoryValueIds = listCategoryValueIds(tissuePanelViews);
          [
            EXPLICIT_BLOOD,
            EXPLICIT_UMBILICAL_CORD_BLOOD,
            EXPLICIT_VENOUS_BLOOD,
          ].forEach((categoryValueId) =>
            tissueCategoryValueIds.includes(categoryValueId)
          );
        });

        /**
         * Selected - blood
         * Selected partial - hematopoietic system
         */
        test("includes hematopoietic system partially selected, blood selected", () => {
          categoryView.panels.forEach((panel) => {
            panel.views.forEach((view) => {
              const selected = view.categoryValueId === INFERRED_BLOOD;
              const selectedPartial =
                view.categoryValueId === INFERRED_HEMATOPOIETIC_SYSTEM;
              expect(view.selected).toEqual(selected);
              expect(view.selectedPartial).toEqual(selectedPartial);
            });
          });
        });
      });

      /**
       * Selected: bone marrow, blood non-specific
       */
      describe("unrelated organ and tissue selected", () => {
        let categoryView: MultiPanelOntologyCategoryView;
        test.beforeAll(() => {
          const updatedKeyedCategoryValues = new Map(KEYED_CATEGORY_VALUES);
          updateSelectCategoryValueSelected(
            updatedKeyedCategoryValues,
            INFERRED_BONE_MARROW,
            true,
            false
          );
          updateSelectCategoryValueSelected(
            updatedKeyedCategoryValues,
            EXPLICIT_BLOOD,
            true,
            false
          );

          const updatedMultiPanelUIState = new Map(multiPanelUIState);
          updateCategoryFilterUIStateSelected(
            updatedMultiPanelUIState,
            CATEGORY_FILTER_ID.TISSUE_CALCULATED,
            [INFERRED_BONE_MARROW, EXPLICIT_BLOOD],
            []
          );

          categoryView = buildMultiPanelCategoryView(
            CATEGORY_FILTER_CONFIGS_BY_ID[
              CATEGORY_FILTER_ID.TISSUE_CALCULATED
            ] as MultiPanelOntologyFilterConfig,
            updatedKeyedCategoryValues,
            updatedMultiPanelUIState,
            new Map<string, string>()
          );
        });

        /**
         * Views: all systems, all organs, tissues under kidney as well as blood specific
         */
        test("includes all systems and organs, filtered tissues and blood non specific", () => {
          expect(categoryView.panels.length).toEqual(3);

          const systemPanelViews =
            // eslint-disable-next-line @typescript-eslint/no-non-null-assertion -- length check above
            categoryView.panels[PANEL_INDEX_SYSTEM]!.views;
          expect(systemPanelViews.length).toEqual(3);
          const systemCategoryValueIds = listCategoryValueIds(systemPanelViews);
          TISSUE_SYSTEMS.forEach((categoryValueId) =>
            systemCategoryValueIds.includes(categoryValueId)
          );

          // eslint-disable-next-line @typescript-eslint/no-non-null-assertion -- length check above
          const organPanelViews = categoryView.panels[PANEL_INDEX_ORGAN]!.views;
          expect(organPanelViews.length).toEqual(7);
          const organCategoryValueIds = listCategoryValueIds(organPanelViews);
          TISSUE_ORGANS.forEach((categoryValueId) =>
            organCategoryValueIds.includes(categoryValueId)
          );

          const tissuePanelViews =
            // eslint-disable-next-line @typescript-eslint/no-non-null-assertion -- length check above
            categoryView.panels[PANEL_INDEX_TISSUE]!.views;
          expect(
            tissuePanelViews.filter((view) => view.visible).length
          ).toEqual(2);
          const tissueCategoryValueIds = listCategoryValueIds(tissuePanelViews);
          [EXPLICIT_BLOOD, EXPLICIT_BONE_MARROW].forEach((categoryValueId) =>
            tissueCategoryValueIds.includes(categoryValueId)
          );
        });
      });

      /**
       * Selected: renal system
       */
      describe("system (containing no organs) selected", () => {
        let categoryView: MultiPanelOntologyCategoryView;
        test.beforeAll(() => {
          const updatedKeyedCategoryValues = new Map(KEYED_CATEGORY_VALUES);
          updateSelectCategoryValueSelected(
            updatedKeyedCategoryValues,
            INFERRED_RENAL_SYSTEM,
            true,
            false
          );

          const updatedMultiPanelUIState = new Map(multiPanelUIState);
          updateCategoryFilterUIStateSelected(
            updatedMultiPanelUIState,
            CATEGORY_FILTER_ID.TISSUE_CALCULATED,
            [INFERRED_RENAL_SYSTEM],
            []
          );

          categoryView = buildMultiPanelCategoryView(
            CATEGORY_FILTER_CONFIGS_BY_ID[
              CATEGORY_FILTER_ID.TISSUE_CALCULATED
            ] as MultiPanelOntologyFilterConfig,
            updatedKeyedCategoryValues,
            updatedMultiPanelUIState,
            new Map<string, string>()
          );
        });

        /**
         * Views: all systems, all organs under renal system, tissues under renal system
         */
        test("includes all systems, no organs and filtered tissues", () => {
          expect(categoryView.panels.length).toEqual(3);

          const systemPanelViews =
            // eslint-disable-next-line @typescript-eslint/no-non-null-assertion -- length check above
            categoryView.panels[PANEL_INDEX_SYSTEM]!.views;
          expect(systemPanelViews.length).toEqual(3);
          const systemCategoryValueIds = listCategoryValueIds(systemPanelViews);
          TISSUE_SYSTEMS.forEach((categoryValueId) =>
            systemCategoryValueIds.includes(categoryValueId)
          );

          // eslint-disable-next-line @typescript-eslint/no-non-null-assertion -- length check above
          const organPanelViews = categoryView.panels[PANEL_INDEX_ORGAN]!.views;
          expect(organPanelViews.filter((view) => view.visible).length).toEqual(
            1
          );
          const organCategoryValueIds = listCategoryValueIds(organPanelViews);
          expect(organCategoryValueIds.includes(INFERRED_KIDNEY)).toBeTruthy();

          const tissuePanelViews =
            // eslint-disable-next-line @typescript-eslint/no-non-null-assertion -- length check above
            categoryView.panels[PANEL_INDEX_TISSUE]!.views;
          expect(
            tissuePanelViews.filter((view) => view.visible).length
          ).toEqual(5);
          const tissueCategoryValueIds = listCategoryValueIds(tissuePanelViews);
          [
            EXPLICIT_BLADDER_LUMEN,
            EXPLICIT_KIDNEY,
            EXPLICIT_RENAL_MEDULLA,
            EXPLICIT_URETER,
            EXPLICIT_URETHRA,
          ].forEach((categoryValueId) =>
            tissueCategoryValueIds.includes(categoryValueId)
          );
        });

        /**
         * Selected - renal system
         */
        test("includes renal system partially selected, blood selected", () => {
          categoryView.panels.forEach((panel) => {
            panel.views.forEach((view) => {
              const selected = view.categoryValueId === INFERRED_RENAL_SYSTEM;
              expect(view.selected).toEqual(selected);
              expect(view.selectedPartial).toBeFalsy();
            });
          });
        });
      });

      /**
       * Selected: blood, blood non-specific, umbilical cord blood, venous blood
       */
      describe("organ and all organ's tissues selected", () => {
        let selectedCategoryValueIds: CategoryValueId[];
        let categoryView: MultiPanelOntologyCategoryView;
        test.beforeAll(() => {
          selectedCategoryValueIds = [
            INFERRED_BLOOD,
            EXPLICIT_BLOOD,
            EXPLICIT_UMBILICAL_CORD_BLOOD,
            EXPLICIT_VENOUS_BLOOD,
          ];

          const updatedKeyedCategoryValues = new Map(KEYED_CATEGORY_VALUES);
          selectedCategoryValueIds.forEach((categoryValueId) =>
            updateSelectCategoryValueSelected(
              updatedKeyedCategoryValues,
              categoryValueId,
              true,
              false
            )
          );

          const updatedMultiPanelUIState = new Map(multiPanelUIState);
          updateCategoryFilterUIStateSelected(
            updatedMultiPanelUIState,
            CATEGORY_FILTER_ID.TISSUE_CALCULATED,
            selectedCategoryValueIds,
            []
          );

          categoryView = buildMultiPanelCategoryView(
            CATEGORY_FILTER_CONFIGS_BY_ID[
              CATEGORY_FILTER_ID.TISSUE_CALCULATED
            ] as MultiPanelOntologyFilterConfig,
            updatedKeyedCategoryValues,
            updatedMultiPanelUIState,
            new Map<string, string>()
          );
        });

        /**
         * Selected - blood, blood non-specific, umbilical cord blood, venous blood
         */
        test("includes blood, blood non-specific, umbilical cord blood, venous blood selected", () => {
          categoryView.panels.forEach((panel) => {
            panel.views.forEach((view) => {
              const selected = selectedCategoryValueIds.includes(
                view.categoryValueId
              );
              expect(view.selected).toEqual(selected);
              expect(view.selectedPartial).toBeFalsy();
            });
          });
        });
      });

      /**
       * Selected: hematopoietic system, bone marrow, spleen, thymus, blood non-specific, umbilical cord blood,
       * venous blood
       */
      describe("system and a mix of system's organs and tissues all selected", () => {
        let selectedCategoryValueIds: CategoryValueId[];
        let categoryView: MultiPanelOntologyCategoryView;
        test.beforeAll(() => {
          selectedCategoryValueIds = [
            INFERRED_HEMATOPOIETIC_SYSTEM,
            INFERRED_BONE_MARROW,
            INFERRED_SPLEEN,
            INFERRED_THYMUS,
            EXPLICIT_BLOOD,
            EXPLICIT_UMBILICAL_CORD_BLOOD,
            EXPLICIT_VENOUS_BLOOD,
          ];

          const updatedKeyedCategoryValues = new Map(KEYED_CATEGORY_VALUES);
          selectedCategoryValueIds.forEach((categoryValueId) =>
            updateSelectCategoryValueSelected(
              updatedKeyedCategoryValues,
              categoryValueId,
              true,
              false
            )
          );

          const updatedMultiPanelUIState = new Map(multiPanelUIState);
          updateCategoryFilterUIStateSelected(
            updatedMultiPanelUIState,
            CATEGORY_FILTER_ID.TISSUE_CALCULATED,
            selectedCategoryValueIds,
            []
          );

          categoryView = buildMultiPanelCategoryView(
            CATEGORY_FILTER_CONFIGS_BY_ID[
              CATEGORY_FILTER_ID.TISSUE_CALCULATED
            ] as MultiPanelOntologyFilterConfig,
            updatedKeyedCategoryValues,
            updatedMultiPanelUIState,
            new Map<string, string>()
          );
        });

        /**
         * Selected - hematopoietic system, blood non-specific, umbilical cord blood, venous blood
         */
        test(
          "includes hematopoietic system, bone marrow, spleen, thymus, blood non-specific, umbilical cord " +
            "blood, venous blood selected",
          () => {
            categoryView.panels.forEach((panel) => {
              panel.views.forEach((view) => {
                const selected = selectedCategoryValueIds.includes(
                  view.categoryValueId
                );
                expect(view.selected).toEqual(selected);
                expect(view.selectedPartial).toBeFalsy();
              });
            });
          }
        );
      });

      /**
       * Selected: hematopoietic system, bone marrow, spleen, thymus, blood non-specific, umbilical cord blood
       */
      describe("system and a mix of system's organs and tissues all selected with one tissue not selected", () => {
        let selectedCategoryValueIds: CategoryValueId[];
        let categoryView: MultiPanelOntologyCategoryView;
        test.beforeAll(() => {
          selectedCategoryValueIds = [
            INFERRED_BONE_MARROW,
            INFERRED_SPLEEN,
            INFERRED_THYMUS,
            EXPLICIT_BLOOD,
            EXPLICIT_UMBILICAL_CORD_BLOOD,
          ];

          const updatedKeyedCategoryValues = new Map(KEYED_CATEGORY_VALUES);
          selectedCategoryValueIds.forEach((categoryValueId) =>
            updateSelectCategoryValueSelected(
              updatedKeyedCategoryValues,
              categoryValueId,
              true,
              false
            )
          );

          updateSelectCategoryValueSelected(
            updatedKeyedCategoryValues,
            INFERRED_HEMATOPOIETIC_SYSTEM,
            false,
            true
          );

          const updatedMultiPanelUIState = new Map(multiPanelUIState);
          updateCategoryFilterUIStateSelected(
            updatedMultiPanelUIState,
            CATEGORY_FILTER_ID.TISSUE_CALCULATED,
            [INFERRED_HEMATOPOIETIC_SYSTEM, ...selectedCategoryValueIds],
            [INFERRED_HEMATOPOIETIC_SYSTEM]
          );

          categoryView = buildMultiPanelCategoryView(
            CATEGORY_FILTER_CONFIGS_BY_ID[
              CATEGORY_FILTER_ID.TISSUE_CALCULATED
            ] as MultiPanelOntologyFilterConfig,
            updatedKeyedCategoryValues,
            updatedMultiPanelUIState,
            new Map<string, string>()
          );
        });

        /**
         * Selected - blood, blood non-specific, umbilical cord blood, venous blood
         * Selected partial - hematopoietic system
         */
        test(
          "includes hematopoietic system partially selected and bone marrow, spleen, thymus, blood " +
            "non-specific, umbilical cord blood selected",
          () => {
            categoryView.panels.forEach((panel) => {
              panel.views.forEach((view) => {
                const selected = selectedCategoryValueIds.includes(
                  view.categoryValueId
                );
                const selectedPartial =
                  view.categoryValueId === INFERRED_HEMATOPOIETIC_SYSTEM;
                expect(view.selected).toEqual(selected);
                expect(view.selectedPartial).toEqual(selectedPartial);
              });
            });
          }
        );
      });

      /**
       * Selected: hematopoietic system, blood, bone marrow, spleen, thymus, blood non-specific, umbilical cord blood
       */
      describe("system, all system's organs and tissues selected with one tissue not selected", () => {
        let selectedCategoryValueIds: CategoryValueId[];
        let selectedPartialCategoryValueIds: CategoryValueId[];
        let categoryView: MultiPanelOntologyCategoryView;
        test.beforeAll(() => {
          selectedCategoryValueIds = [
            INFERRED_BONE_MARROW,
            INFERRED_SPLEEN,
            INFERRED_THYMUS,
            EXPLICIT_BLOOD,
            EXPLICIT_UMBILICAL_CORD_BLOOD,
          ];

          selectedPartialCategoryValueIds = [
            INFERRED_HEMATOPOIETIC_SYSTEM,
            INFERRED_BLOOD,
          ];

          const updatedKeyedCategoryValues = new Map(KEYED_CATEGORY_VALUES);
          selectedCategoryValueIds.forEach((categoryValueId) =>
            updateSelectCategoryValueSelected(
              updatedKeyedCategoryValues,
              categoryValueId,
              true,
              false
            )
          );

          selectedPartialCategoryValueIds.forEach((categoryValueId) =>
            updateSelectCategoryValueSelected(
              updatedKeyedCategoryValues,
              categoryValueId,
              false,
              true
            )
          );

          const updatedMultiPanelUIState = new Map(multiPanelUIState);
          updateCategoryFilterUIStateSelected(
            updatedMultiPanelUIState,
            CATEGORY_FILTER_ID.TISSUE_CALCULATED,
            [...selectedPartialCategoryValueIds, ...selectedCategoryValueIds],
            [...selectedPartialCategoryValueIds]
          );

          categoryView = buildMultiPanelCategoryView(
            CATEGORY_FILTER_CONFIGS_BY_ID[
              CATEGORY_FILTER_ID.TISSUE_CALCULATED
            ] as MultiPanelOntologyFilterConfig,
            updatedKeyedCategoryValues,
            updatedMultiPanelUIState,
            new Map<string, string>()
          );
        });

        /**
         * Selected - blood, blood non-specific, umbilical cord blood, venous blood
         * Selected partial - hematopoietic system
         */
        test(
          "includes hematopoietic system and blood partially selected and bone marrow, spleen, thymus, blood " +
            "non-specific, umbilical cord blood selected",
          () => {
            categoryView.panels.forEach((panel) => {
              panel.views.forEach((view) => {
                const selected = selectedCategoryValueIds.includes(
                  view.categoryValueId
                );
                const selectedPartial =
                  selectedPartialCategoryValueIds.includes(
                    view.categoryValueId
                  );
                expect(view.selected).toEqual(selected);
                expect(view.selectedPartial).toEqual(selectedPartial);
              });
            });
          }
        );
      });
    });
  });

  describe("Curated Ontology Category", () => {
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
      test(`de-selects leaf ${ONTOLOGY_ID_HUMAN_CARNEGIE_CS1}, no siblings`, () => {
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
          ).source
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
      test(`de-selects leaf ${ONTOLOGY_ID_HUMAN_CARNEGIE_CS1}, with siblings`, () => {
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
          ).source
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
      test(`de-selects node ${ONTOLOGY_ID_HUMAN_EMBRYONIC_HUMAN}, sibling unselected`, () => {
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
          ).source
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
      test(`de-selects node ${ONTOLOGY_ID_HUMAN_EMBRYONIC_HUMAN}, sibling selected`, () => {
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
          ).source
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
      test(`selects leaf ${ONTOLOGY_ID_HUMAN_CARNEGIE_CS1}, no siblings currently selected`, () => {
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
          ).source
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
      test(`selects leaf ${ONTOLOGY_ID_HUMAN_INFANT}, some siblings currently selected`, () => {
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
          ).source
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
      test(`selects leaf ${ONTOLOGY_ID_HUMAN_CHILD}, all siblings currently selected`, () => {
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
          ).source
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
      test(`selects leaf ${ONTOLOGY_ID_HUMAN_CARNEGIE_CS1}, all siblings and aunt currently selected`, () => {
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
          ).source
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
      test(`selects node ${ONTOLOGY_ID_HUMAN_EMBRYONIC_HUMAN}, nothing currently selected`, () => {
        const nextFilters = buildNextOntologyCategoryFilters(
          CATEGORY_FILTER_ID.DEVELOPMENT_STAGE,
          ONTOLOGY_ID_HUMAN_EMBRYONIC_HUMAN,
          [], // No filters
          CATEGORY_VALUE_KEYS,
          (
            CATEGORY_FILTER_CONFIGS_BY_ID[
              CATEGORY_FILTER_ID.DEVELOPMENT_STAGE
            ] as CuratedOntologyCategoryFilterConfig
          ).source
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
      test(`selects node ${ONTOLOGY_ID_HUMAN_EMBRYONIC_HUMAN}, child currently selected`, () => {
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
          ).source
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
      test(`selects node ${ONTOLOGY_ID_HUMAN_EMBRYONIC_HUMAN}, sibling currently selected`, () => {
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
          ).source
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
      test(`de-selects leaf ${ONTOLOGY_ID_HUMAN_CARNEGIE_CS1}, ancestor currently selected`, () => {
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
          ).source
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

/**
 * Return the category value IDs for the given array of multi-panel views.
 * @param views - Array of select category value views to map category value IDs from.
 */
function listCategoryValueIds(
  views: SelectCategoryValueView[]
): CategoryValueId[] {
  return views.map((view) => view.categoryValueId);
}

/**
 * Update the selected and selected partial values for the category value view with the given key.
 * @param valueViews - Map of select category value views keyed by category value ID.
 * @param categoryValueId - ID of category value view to update.
 * @param selected - True if category value view is to be updated to selected.
 * @param selectedPartial - True if category value view is to be updated to selected partial.
 */
function updateCategoryValueViewSelected(
  valueViews: Map<CategoryValueId, SelectCategoryValueView>,
  categoryValueId: CategoryValueId,
  selected: boolean,
  selectedPartial: boolean
) {
  const valueView = {
    ...valueViews.get(categoryValueId),
  } as SelectCategoryValueView;
  valueView.selected = selected;
  valueView.selectedPartial = selectedPartial;
  valueViews.set(categoryValueId, valueView);
}

/**
 * Update the selected and selected partial values for the category value view with the given key.
 * @param multiPanelUIState - Complete multi-panel UI state.
 * @param categoryFilterId - The category filter to update.
 * @param selected - True if category value view is to be updated to selected.
 * @param selectedPartial - True if category value view is to be updated to selected partial.
 */
function updateCategoryFilterUIStateSelected(
  multiPanelUIState: MultiPanelUIState,
  categoryFilterId: CATEGORY_FILTER_ID,
  selected: CategoryValueId[],
  selectedPartial: CategoryValueId[]
) {
  const categoryFilterUIState = multiPanelUIState.get(categoryFilterId);
  if (!categoryFilterUIState) {
    console.error(
      `MultiPanelCategoryFilterUIState not found for ${categoryFilterId}`
    );
    return;
  }
  const updatedCategoryFilterUIState = {
    ...categoryFilterUIState,
  };
  updatedCategoryFilterUIState.selected = selected;
  updatedCategoryFilterUIState.selectedPartial = selectedPartial;
  multiPanelUIState.set(categoryFilterId, updatedCategoryFilterUIState);
}

/**
 * Update the selected and selected partial values for the category value with the given key.
 * @param selectCategoryValues - Map of select category value keyed by category value ID.
 * @param categoryValueId - ID of category value view to update.
 * @param selected - True if category value is to be updated to selected.
 * @param selectedPartial - True if category value is to be updated to selected partial.
 */
function updateSelectCategoryValueSelected(
  selectCategoryValues: Map<CategoryValueId, SelectCategoryValue>,
  categoryValueId: CategoryValueId,
  selected: boolean,
  selectedPartial: boolean
) {
  const selectCategoryValue = {
    ...selectCategoryValues.get(categoryValueId),
  } as SelectCategoryValue;
  selectCategoryValue.selected = selected;
  selectCategoryValue.selectedPartial = selectedPartial;
  selectCategoryValues.set(categoryValueId, selectCategoryValue);
}
