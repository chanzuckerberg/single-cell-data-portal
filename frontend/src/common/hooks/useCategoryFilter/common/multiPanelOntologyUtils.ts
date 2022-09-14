import { Filters, Row } from "react-table";
import { track } from "src/common/analytics";
import { EMPTY_ARRAY } from "src/common/constants/utils";
import {
  buildNextSelectCategoryFilters,
  buildSelectCategoryValueViews,
} from "src/common/hooks/useCategoryFilter/common/selectUtils";
import {
  CATEGORY_FILTER_CONFIGS_BY_ID,
  LABEL_SUFFIX_NON_SPECIFIC,
} from "src/components/common/Filter/common/constants";
import {
  Categories,
  CategoryFilter,
  CategoryFilterConfig,
  CategoryFilterPanelConfig,
  CategoryValueId,
  CATEGORY_FILTER_ID,
  KeyedSelectCategoryValue,
  MultiPanelCategoryFilterUIState,
  MultiPanelOntologyCategoryView,
  MultiPanelOntologyFilterConfig,
  MultiPanelSelectedUIState,
  MultiPanelUINode,
  MultiPanelUIState,
  OntologyDescendants,
  OntologyPanelCategoryView,
  ON_FILTER_SOURCE,
  OrFilterPrefix,
  SelectCategoryValue,
  SelectCategoryValueView,
} from "src/components/common/Filter/common/entities";
import { listOntologyTreeIds } from "src/components/common/Filter/common/utils";
import { getCategoryFilter, sortCategoryValueViews } from "./utils";

/**
 * Utils specific to multi-panel category filters.
 */

/**
 * Model of a category value's selected state. This is an "intermediate" object that is used during recursive
 * selected/partially selected calculations.
 */
interface PartialSelectedCheck {
  isMaskedByChild: boolean;
  isImpliedSelectedByChildren: boolean;
}

/**
 * Update label of views that also appear in parent panels.
 * @param panelCategoryValueViews - Views to display in the current panel.
 * @param parentUINodes - All UI nodes in parent panels
 */
function applyNonSpecificLabel(
  panelCategoryValueViews: SelectCategoryValueView[],
  parentUINodes: MultiPanelUINode[]
) {
  // Collect the category values of every parent.
  const parentCategoryValueIds = parentUINodes.map(
    (parentUINode: MultiPanelUINode) => parentUINode.categoryValueId
  );

  // Check if any view has a label collision (e.g. "blood" in the organ panel and "blood" in the tissue panel) and if
  // so, update the view to include ", non-specific" terminology in its label.
  panelCategoryValueViews.forEach(
    (panelCategoryValueView: SelectCategoryValueView) => {
      if (
        isLabelCollision(
          panelCategoryValueView.categoryValueId,
          parentCategoryValueIds
        )
      ) {
        panelCategoryValueView.label = `${panelCategoryValueView.label}${LABEL_SUFFIX_NON_SPECIFIC}`;
      }
    }
  );
}

/**
 * Build filter value with an inferred prefix. Used by ontology-aware filter categories that require filtering by both
 * inferred and explicit values (e.g. tissue).
 * @param ontologyTermId - The term ID to prefix with inferred.
 * @returns String containing explicit prefix prepended to the given ontology term ID.
 */
export function buildExplicitOntologyTermId(ontologyTermId: string): string {
  return `${OrFilterPrefix.EXPLICIT}:${ontologyTermId}`;
}

/**
 * Build filter value with an explicit prefix. Used by ontology-aware filter categories that require filtering by both
 * inferred and explicit values (e.g. tissue).
 * @param ontologyTermId - The term ID to prefix with inferred.
 * @returns String containing inferred prefix prepended to the given ontology term ID.
 */
export function buildInferredOntologyTermId(ontologyTermId: string): string {
  return `${OrFilterPrefix.INFERRED}:${ontologyTermId}`;
}

/**
 * Build view model of multi-panel category.
 * @param config - Config model of ontology category.
 * @param categoryValuesByValue - Internal filter model of single or multiselect category.
 * @param multiPanelUIState - Current set of category values that the user has selected in multi-panel category filters.
 * @param ontologyTermLabelsById - Set of ontology term labels keyed by term ID, used to determine labels for ontology
 * @returns Multi-panel category view model.
 */
export function buildMultiPanelCategoryView(
  config: MultiPanelOntologyFilterConfig,
  categoryValuesByValue: KeyedSelectCategoryValue,
  multiPanelUIState: MultiPanelUIState,
  ontologyTermLabelsById: Map<string, string>
): MultiPanelOntologyCategoryView {
  const { categoryFilterId } = config;
  const categoryFilterUIState = multiPanelUIState.get(categoryFilterId);
  if (!categoryFilterUIState) {
    console.log(
      `Multi-panel category filter state not found for category ${categoryFilterId}`
    );
    return {
      categoryFilterId: categoryFilterId,
      label: config.label,
      panels: EMPTY_ARRAY,
      selectedViews: EMPTY_ARRAY,
    };
  }

  // Build value view models for each panel, making sure cross-panel filtering is applied.
  const ontologyPanelCategoryViews = buildOntologyPanelCategoryViews(
    config,
    categoryValuesByValue,
    categoryFilterUIState,
    ontologyTermLabelsById
  );

  // Determine set of selected values for this multi-panel category
  const allCategoryValueViews = ontologyPanelCategoryViews
    .map(
      (ontologyPanelCategoryView: OntologyPanelCategoryView) =>
        ontologyPanelCategoryView.views
    )
    .flat();
  const selectedViews = listMultiPanelSelectedViews(
    allCategoryValueViews,
    categoryFilterUIState
  );

  // Build view model of multi-panel category.
  return {
    categoryFilterId: categoryFilterId,
    label: config.label,
    panels: ontologyPanelCategoryViews,
    selectedViews: selectedViews,
  };
}

/**
 * Build up react-table filters model from the given multi-panel UI state. Do not add partially selected values as
 * only the current set of filters applied by react-table for this category filter should be returned.
 * @param multiPanelUIState - ID of category filter to filters.
 * @returns Set of selected values in the format expected by react-table.
 */
export function buildMultiPanelCurrentFilters<T extends Categories>(
  multiPanelUIState: MultiPanelUIState
): Filters<T> {
  return [...multiPanelUIState.keys()].reduce(
    (accum: CategoryFilter[], categoryFilterId: CATEGORY_FILTER_ID) => {
      // Don't add category to filters if it has no selected values.
      const categoryFilterUIState = multiPanelUIState.get(categoryFilterId);
      if (!categoryFilterUIState) {
        console.log(
          `Category filter UI state not found for category filter ${categoryFilterId}`
        );
        return accum; // Error state
      }
      const { selected: selectedCategoryValues } = categoryFilterUIState;
      if (!selectedCategoryValues || !selectedCategoryValues.length) {
        return accum;
      }

      // Otherwise, this category has selected values: add them!
      accum.push({
        id: categoryFilterId,
        value: selectedCategoryValues,
      });
      return accum;
    },
    [] as CategoryFilter[]
  );
}

/**
 * Build serializable model of selected and partially selected state of each multi-panel category filters. Used by
 * callees of hook to save filter state to local storage.
 * @param multiPanelUIState - Current set of category values that the user has selected in multi-panel category filters.
 * @returns Serializable model of selected and partially selected state of each multi-panel category filter.
 */
export function buildMultiPanelSelectedUIState(
  multiPanelUIState: MultiPanelUIState
) {
  return [...multiPanelUIState.keys()].reduce((accum, categoryFilterId) => {
    const categoryFilterUIState = multiPanelUIState.get(categoryFilterId);
    accum[categoryFilterId] = {
      // Adding chaining (and nullish coalescing) here to satisfy "possibly undefined" errors but as we are reducing
      // over the keys in the multi-panel UI state, these values are actually guaranteed to be defined.
      selected: categoryFilterUIState?.selected ?? [],
      selectedPartial: categoryFilterUIState?.selectedPartial ?? [],
    };
    return accum;
  }, {} as MultiPanelSelectedUIState);
}

/**
 * Build up the base UI model for each category filter.
 * @param originalRows - Original result set before filtering.
 * @param categoryFilterIds - Set of category IDs to include for this filter instance.
 * @param initialMultiPanelSelectedUIState - Selected state of category to set as an initial state.
 * @returns Map of UI nodes keyed by category value, facilitates easy lookups of UI node state.
 */
export function buildMultiPanelUIState<T extends Categories>(
  originalRows: Row<T>[],
  categoryFilterIds: Set<CATEGORY_FILTER_ID>,
  initialMultiPanelSelectedUIState: MultiPanelSelectedUIState
): MultiPanelUIState {
  return Array.from(categoryFilterIds.values()).reduce(
    (
      accum: Map<CATEGORY_FILTER_ID, MultiPanelCategoryFilterUIState>,
      categoryFilterId: CATEGORY_FILTER_ID
    ) => {
      // Ignore categories that are not multi-panels.
      const config = CATEGORY_FILTER_CONFIGS_BY_ID[categoryFilterId];
      if (!isMultiPanelCategoryFilterConfig(config)) {
        return accum;
      }

      // Determine the set of values for each panel.
      const categoryValueIdsByPanel = keyCategoryValueIdsByPanel(
        config,
        originalRows
      );

      // Build up parent/children relationships for each category value.
      const uiHierarchyByCategoryValue = buildUINodesByCategoryValueId(
        categoryValueIdsByPanel,
        config.descendants
      );

      // Determine the initial selected state for this category.
      const categorySelectedUIState =
        initialMultiPanelSelectedUIState?.[categoryFilterId];

      accum.set(categoryFilterId, {
        selected: categorySelectedUIState?.selected ?? [],
        selectedPartial: categorySelectedUIState?.selectedPartial ?? [],
        uiNodesByCategoryValueId: uiHierarchyByCategoryValue,
      });

      return accum;
    },
    new Map<CATEGORY_FILTER_ID, MultiPanelCategoryFilterUIState>()
  );
}

/**
 * Build updated set of selected filters for the given single or multiselect category and the selected category values.
 * If this is a remove action from a tag, remove selected descendants as well, if applicable.
 * @param config - Configuration model of selected category.
 * @param selectedValue - Selected category value ID to use as selected value.
 * @param source - Location where select/remove was triggered, either filter menu or selected tag.
 * @param categoryFilterUIState - UI model of partially/selected values in the category filter.
 * @param currentFilters - Set of selected values in categor filter with no overrides/filtering/partial logic applied.
 */
function buildNextMultiPanelCategoryFilters<T extends Categories>(
  config: MultiPanelOntologyFilterConfig,
  selectedValue: CategoryValueId,
  source: ON_FILTER_SOURCE,
  categoryFilterUIState: MultiPanelCategoryFilterUIState,
  currentFilters: Filters<T>
): CategoryValueId[] {
  // Determine selected set of values for this category filter based on the current selected values for this multi-panel
  // category filter; toggle current selected values.
  const selectedCategoryFilters = buildNextSelectCategoryFilters(
    config,
    selectedValue,
    currentFilters
  );

  // If the source of the filter action wasn't a tag, we can simply toggle the value and exit here.
  if (source !== ON_FILTER_SOURCE.TAG) {
    return selectedCategoryFilters;
  }

  // If the source of the filter action was a selected tag, we know a remove action has occurred. Clear out any
  // children of the removed tag unless the children themselves have another that is partially/selected.
  return onRemoveMultiPanelCategoryValueTag(
    selectedValue,
    selectedCategoryFilters,
    categoryFilterUIState,
    config.descendants
  );
}

/**
 * Build view models for each panel in the multi-panel category.
 * @param config - Config model of ontology category.
 * @param categoryValuesByValue - Internal filter model of multi-panel category.
 * @param categoryFilterUIState - UI model of partially/selected values in the multi-panel category.
 * @param ontologyTermLabelsById - Set of ontology term labels keyed by term ID, used to determine labels for ontology.
 * @returns View model of multi-panel category.
 */
function buildOntologyPanelCategoryViews(
  config: MultiPanelOntologyFilterConfig,
  categoryValuesByValue: KeyedSelectCategoryValue,
  categoryFilterUIState: MultiPanelCategoryFilterUIState,
  ontologyTermLabelsById: Map<string, string>
): OntologyPanelCategoryView[] {
  const { descendants, panels: panelConfigs } = config;
  return panelConfigs.reduce((accum, panelConfig, panelIndex) => {
    // Determine the indices of the panels that are parents to the current panel.
    const parentPanelIndices = listParentPanelIndices(
      panelIndex,
      panelConfigs.length
    );

    // Determine if this panel is filtered by checking if there are any selected values in this panel's parents.
    const selectedParentUINodes = listSelectedParentUINodes(
      parentPanelIndices,
      categoryFilterUIState,
      descendants
    );

    // If panel is filtered (that is, there are selected values in this panel's parents), only include values in this
    // panel that are descendants of the selected parent values. If there are no selected values in parent panels,
    // include all values in this panel.
    const panelUINodes = listPanelUINodes(panelIndex, categoryFilterUIState);
    const includeUINodes = selectedParentUINodes.length
      ? filterPanelUINodes(
          panelUINodes,
          selectedParentUINodes,
          categoryFilterUIState,
          descendants
        )
      : panelUINodes;

    // Collect the select category values for the values to be included.
    const selectCategoryValues = includeUINodes
      .map((uiNode: MultiPanelUINode) =>
        categoryValuesByValue.get(uiNode.categoryValueId)
      )
      .filter(
        (
          selectCategoryValue: SelectCategoryValue | undefined
        ): selectCategoryValue is SelectCategoryValue => !!selectCategoryValue
      );

    // Build views for each selected category value to be included.
    const panelViews: SelectCategoryValueView[] = buildSelectCategoryValueViews(
      config,
      selectCategoryValues,
      ontologyTermLabelsById
    );

    // Apply non-specific labels to values that appear in this panel as well as parent panels.
    if (panelConfig.sourceKind === "EXPLICIT_ONLY") {
      // Find all UI nodes of parents.
      const parentUINodes = listPanelsUINodes(
        parentPanelIndices,
        categoryFilterUIState
      );

      // Update label collisions.
      applyNonSpecificLabel(panelViews, parentUINodes);
    }

    // Sort views by label.
    panelViews.sort(sortCategoryValueViews("label"));

    // Build panel view.
    accum.push({
      isSearchMultiselect: panelConfig.searchKind === "SEARCH_MULTI_SELECT",
      label: panelConfig.label,
      views: panelViews,
    });

    return accum;
  }, [] as OntologyPanelCategoryView[]);
}

/**
 * Build up the base UI nodes for each category value. This base model is updated as the UI state changes.
 * @param categoryValueIdsByPanel - Category values grouped by panel.
 * @param descendants - Map of descendants keyed by ancestor.
 * @returns Map of UI nodes keyed by category value, facilitates easy lookups of UI node state.
 */
export function buildUINodesByCategoryValueId(
  categoryValueIdsByPanel: CategoryValueId[][],
  descendants: OntologyDescendants
): Map<CategoryValueId, MultiPanelUINode> {
  return categoryValueIdsByPanel.reduce(
    (
      uiAccum: Map<CategoryValueId, MultiPanelUINode>,
      categoryValueIds: CategoryValueId[],
      index: number
    ) => {
      // Add all values in this panel to the map of values.
      categoryValueIds.forEach((categoryValueId) =>
        uiAccum.set(categoryValueId, {
          categoryValueId: categoryValueId,
          panelIndex: index,
          uiChildren: [],
          uiParents: [],
        })
      );

      // Build parent child relationships for each value in this panel.
      const parentCategoryValueIdsByPanel = categoryValueIdsByPanel.slice(
        0,
        index
      );
      categoryValueIds.forEach((categoryValueId) => {
        linkParentsAndChildren(
          categoryValueId,
          parentCategoryValueIdsByPanel,
          uiAccum,
          descendants
        );
      });

      return uiAccum;
    },
    new Map<CategoryValueId, MultiPanelUINode>()
  );
}

/**
 * Determine the UI nodes to include for display. If a parent panel has selected (or possibly partially selected)
 * values, then the UI nodes for this panel need to be limited to descendants of the parent selected values.
 * @param panelUINodes - All UI nodes in the panel.
 * @param selectedParentUINodes - UI nodes in parent panels that are selected.
 * @param categoryFilterUIState - UI model of partially/selected values in the multi-panel category.
 * @param descendants - Map of descendants keyed by ancestor.
 * @returns An array of UI nodes for the panel that are to be included for display.
 */
function filterPanelUINodes(
  panelUINodes: MultiPanelUINode[],
  selectedParentUINodes: MultiPanelUINode[],
  categoryFilterUIState: MultiPanelCategoryFilterUIState,
  descendants: OntologyDescendants
): MultiPanelUINode[] {
  return panelUINodes.filter((panelUINode: MultiPanelUINode) => {
    const { categoryValueId } = panelUINode;

    // If the node is selected or partially selected, it is always displayed.
    if (
      categoryFilterUIState.selected.includes(categoryValueId) ||
      categoryFilterUIState.selectedPartial.includes(categoryValueId)
    ) {
      return true;
    }

    // Otherwise, only display values that are descendants of selected values.
    return selectedParentUINodes.some(
      (selectedParentUINode: MultiPanelUINode) =>
        isDescendant(
          categoryValueId,
          selectedParentUINode.categoryValueId,
          descendants
        )
    );
  });
}

/**
 * Return the UI children for the given node. UI children are initialized to an empty array on hook init but require
 * truthy checking due to model being backed by a Map. Encapsulate the truthy checking here.
 * @param categoryValueId - Category value to return UI children of.
 * @param uiNodesByCategoryValueId - UI model containing parent/child relationships.
 */
function getUIChildren(
  categoryValueId: CategoryValueId,
  uiNodesByCategoryValueId: Map<CategoryValueId, MultiPanelUINode>
): CategoryValueId[] {
  return uiNodesByCategoryValueId.get(categoryValueId)?.uiChildren ?? [];
}

/**
 * Return the UI parents for the given node. UI parents are initialized to an empty array on hook init but require
 * truthy checking due to model being backed by a Map. Encapsulate the truthy checking here.
 * @param categoryValueId - Category value to return UI parents of.
 * @param uiNodesByCategoryValueId - UI model containing parent/child relationships.
 */
function getUIParents(
  categoryValueId: CategoryValueId,
  uiNodesByCategoryValueId: Map<CategoryValueId, MultiPanelUINode>
): CategoryValueId[] {
  return uiNodesByCategoryValueId.get(categoryValueId)?.uiParents ?? [];
}

/**
 * Determine if there is a ancestor/descendant relationship between the given category values. A descendant can
 * also be an explicit ontology term ID.
 * @param descendantCategoryValueId - Category value to check if it's a descendant.
 * @param ancestorCategoryValueId - Category value to check if it's an ancestor.
 * @param descendants - Map of descendants keyed by ancestor.
 * @returns True if there is an ancestor/descendant relationship or inferred/explicit relationship between the given
 * category IDs.
 */
function isDescendant(
  descendantCategoryValueId: CategoryValueId,
  ancestorCategoryValueId: CategoryValueId,
  descendants: OntologyDescendants
): boolean {
  // Check ancestor/descendant relationship.
  const ontologyTermId = removeOntologyTermIdPrefix(descendantCategoryValueId);
  const ancestorOntologyTermId = removeOntologyTermIdPrefix(
    ancestorCategoryValueId
  );
  const isDescendantTerm = (descendants[ancestorOntologyTermId] ?? []).includes(
    ontologyTermId
  );

  // Check inferred/explicit relationship.
  const isExplicitTerm = isExplicitTermOfInferredTerm(
    descendantCategoryValueId,
    ancestorCategoryValueId
  );

  return isDescendantTerm || isExplicitTerm;
}

/**
 * Determine if there is an inferred/explicit relationship between the given category values. For example, "E:UBERON:X"
 * is an explicit version of "I:UBERON:X". Specifically, "blood, non-specific" is an explicit version of "blood".
 * @param explicitCategoryValueId - Category value to check if it's an explicit version of the inferred category value.
 * @param inferredCategoryValueId - Category value to check if it's an inferred version of the explicit category value.
 * @returns True if there is an inferred/explicit relationship between the given category values.
 */
function isExplicitTermOfInferredTerm(
  explicitCategoryValueId: CategoryValueId,
  inferredCategoryValueId: CategoryValueId
): boolean {
  const ontologyTermId = removeOntologyTermIdPrefix(explicitCategoryValueId);
  const inferredOntologyTermId = removeOntologyTermIdPrefix(
    inferredCategoryValueId
  );
  return (
    isExplicitOntologyTermId(explicitCategoryValueId) &&
    !isExplicitOntologyTermId(inferredCategoryValueId) &&
    ontologyTermId === inferredOntologyTermId
  );
}

/**
 * Check if the given category value is an explicit value and not an inferred value.
 * @param categoryValueId - Category value to check.
 * @returns True if category value is an explicit value.
 */
function isExplicitOntologyTermId(categoryValueId: CategoryValueId): boolean {
  const prefix = removeOntologyTermId(categoryValueId);
  return prefix === OrFilterPrefix.EXPLICIT;
}

/**
 * Determine if the given category values are in the same panel.
 * @param categoryValueId - First category value to check panel index of.
 * @param compareToCategoryValueId - Second category value to check panel index of.
 * @param uiNodesByCategoryValueId - Map of category value ID to UI parents, UI children and panel index.
 * @returns True if the given category values are in the same panel.
 */
function isInSamePanel(
  categoryValueId: CategoryValueId,
  compareToCategoryValueId: CategoryValueId,
  uiNodesByCategoryValueId: Map<CategoryValueId, MultiPanelUINode>
): boolean {
  const categoryValueUINode = uiNodesByCategoryValueId.get(categoryValueId);
  if (!categoryValueUINode) {
    return false;
  }
  const compareToUINode = uiNodesByCategoryValueId.get(
    compareToCategoryValueId
  );
  if (!compareToUINode) {
    return false;
  }
  return categoryValueUINode.panelIndex === compareToUINode.panelIndex;
}

/**
 * Determine if the category value ID also exists in the parent set of category value IDs.
 * @param categoryValueId - Category value to check for against parent category value IDs.
 * @param parentCategoryValueIds - All category values in the parent panels.
 * @returns True if category value is also present in a parent panel.
 */
function isLabelCollision(
  categoryValueId: CategoryValueId,
  parentCategoryValueIds: CategoryValueId[]
): boolean {
  return parentCategoryValueIds.some((parentCategoryValueId: CategoryValueId) =>
    isExplicitTermOfInferredTerm(categoryValueId, parentCategoryValueId)
  );
}

/**
 * Determine if a category value should be marked as selected or partially selected.
 * @param categoryValueId - The category value we are working on.
 * @param selectedCategoryValueIds - All explicitly selected category value ids in all panels.
 * @param uiNodesByCategoryValueId - Map of categoryValueID to UI parents and children.
 * @returns True if category value is partially selected.
 */
function isPartiallySelectedCategoryValue(
  categoryValueId: CategoryValueId,
  selectedCategoryValueIds: CategoryValueId[],
  uiNodesByCategoryValueId: Map<CategoryValueId, MultiPanelUINode>
): boolean {
  return isPartiallySelectedCategoryValueInternal(
    categoryValueId,
    selectedCategoryValueIds,
    uiNodesByCategoryValueId
  ).isMaskedByChild;
}

/**
 * Recursive function used when determining if a category value should be marked as selected or partially selected.
 * @param categoryValueId - The category value we are working on.
 * @param selectedCategoryValueIds - All explicitly selected category value ids in all panels.
 * @param uiNodesByCategoryValueId - Map of categoryValueID to UI parents and children.
 * @returns Model of selected state for the given category value.
 */
function isPartiallySelectedCategoryValueInternal(
  categoryValueId: CategoryValueId,
  selectedCategoryValueIds: CategoryValueId[],
  uiNodesByCategoryValueId: Map<CategoryValueId, MultiPanelUINode>
): PartialSelectedCheck {
  let isMaskedByChild;
  let isImpliedSelectedByChildren;

  // Get this node's UI children
  const uiChildren = getUIChildren(categoryValueId, uiNodesByCategoryValueId);

  // If you have no children, you are neither implied by your children nor partial (masked by a child).
  if (!uiChildren.length) {
    isMaskedByChild = false;
    isImpliedSelectedByChildren = false;
    return { isImpliedSelectedByChildren, isMaskedByChild };
  }

  // You have children so collect and check them.
  const childrenCheckByCategoryValueId = uiChildren.reduce((accum, uiChild) => {
    const childCheck = isPartiallySelectedCategoryValueInternal(
      uiChild,
      selectedCategoryValueIds,
      uiNodesByCategoryValueId
    );

    accum.set(uiChild, childCheck);

    return accum;
  }, new Map<CategoryValueId, PartialSelectedCheck>());

  // Evaluate your children...

  // Determine your isImpliedSelectedByChildren - you are implied selected if all of your children are either selected
  // or implied selected.
  isImpliedSelectedByChildren = uiChildren.every((uiChild) => {
    return isSelectedOrImplied(
      uiChild,
      selectedCategoryValueIds,
      childrenCheckByCategoryValueId
    );
  });

  // Determine your isMaskedByChild - you are masked by a child (partially selected) if some but not all of your
  // children are selected or implied selected.
  const isAnySelectedOrImplied = uiChildren.some((uiChild) => {
    return isSelectedOrImplied(
      uiChild,
      selectedCategoryValueIds,
      childrenCheckByCategoryValueId
    );
  });

  const isEverySelectedOrImplied = uiChildren.every((uiChild) => {
    return isSelectedOrImplied(
      uiChild,
      selectedCategoryValueIds,
      childrenCheckByCategoryValueId
    );
  });

  const isAnyMasked = uiChildren.some((uiChild) => {
    return isMasked(uiChild, childrenCheckByCategoryValueId);
  });

  isMaskedByChild =
    isAnyMasked || (isAnySelectedOrImplied && !isEverySelectedOrImplied);

  return { isImpliedSelectedByChildren, isMaskedByChild };
} /* FIN */

/**
 * Returns true if the UI node has a descendant node that is selected.
 * @param uiChild - Category value to check if masked.
 * @param childrenChecksByCategoryValueId - Model of partially selected states keyed by category value.
 * @returns True if UI child is masked by a selected child value.
 */
function isMasked(
  uiChild: CategoryValueId,
  childrenChecksByCategoryValueId: Map<CategoryValueId, PartialSelectedCheck>
): boolean {
  const childCheck = childrenChecksByCategoryValueId.get(uiChild);
  if (!childCheck) {
    console.log(`Child check not found for ${uiChild}`);
    return false; // Error state - can't determine if masked; return false.
  }
  return childCheck.isMaskedByChild;
}

/**
 * Determine if the given category config is a multi-panel category.
 * @param config - Config model of category, either an ontology category config or a base category config.
 * @returns True if category config is for a multi-panel category.
 */
export function isMultiPanelCategoryFilterConfig(
  config: CategoryFilterConfig
): config is MultiPanelOntologyFilterConfig {
  return config.viewKind === "MULTI_PANEL";
}

/**
 * Check if all children of the given UI node are either selected or selected implied. (e.g. "blood" in the organ panel
 * selected implied if all blood descendants in the tissue panel are selected).
 * @param uiChild - Category value to check if masked.
 * @param selectedCategoryValueIds - All explicitly selected category value ids in all panels.
 * @param childrenChecksByCategoryValueId - Model of partially selected states keyed by category value.
 * @returns True if UI child is masked by a selected child value.
 */
function isSelectedOrImplied(
  uiChild: CategoryValueId,
  selectedCategoryValueIds: CategoryValueId[],
  childrenChecksByCategoryValueId: Map<CategoryValueId, PartialSelectedCheck>
): boolean {
  const childSelected = selectedCategoryValueIds.includes(uiChild);
  const childCheck = childrenChecksByCategoryValueId.get(uiChild);
  if (!childCheck) {
    console.log(`Child check not found for ${uiChild}`);
    return false; // Error state - can't determine if selected/implied; return false.
  }
  const childImplied = childCheck.isImpliedSelectedByChildren;
  return childSelected || childImplied;
}

/**
 * Determine if the given category value is to be included as a selected tag.
 * @param categoryValueId - Category value to check if it can be included as a selected tag.
 * @param selected - Array of category values that are selected.
 * @param selectedPartial - Array of category values that are partially selected.
 * @param uiNodesByCategoryValueId - Map of categoryValueID to UI parents and children.
 * @returns True if category value is to be included as a selected tag.
 */
function isSelectedViewTagVisible(
  categoryValueId: CategoryValueId,
  selected: CategoryValueId[],
  selectedPartial: CategoryValueId[],
  uiNodesByCategoryValueId: Map<CategoryValueId, MultiPanelUINode>
): boolean {
  // Grab the parents for this selected view.
  const uiParents = getUIParents(categoryValueId, uiNodesByCategoryValueId);

  // Grab the selected parents.
  const selectedParents = uiParents.filter((uiParent: CategoryValueId) =>
    selected.includes(uiParent)
  );

  // If there are no selected parents, check all ancestors to see if they are included. For example, hema system,
  // blood, non-specific, umbilical cord blood, venous blood.
  if (!selectedParents.length) {
    return uiParents.every((uiParent: CategoryValueId) =>
      isSelectedViewTagVisible(
        uiParent,
        selected,
        selectedPartial,
        uiNodesByCategoryValueId
      )
    );
  }

  // If any parent is partially selected, selected view is visible. The view is visible if any parent is partially
  // selected.
  return selectedParents.some((uiParent) => selectedPartial.includes(uiParent));
}

/**
 * Group all values from the original rows for the given multi-panel category filter and group by panel.
 * @param config - Configuration model of selected category.
 * @param originalRows - Original result set before filtering.
 * @returns Array of arrays, one outer array for each panel in the multi-panel category filter, and in inner array
 * containing the category values for the panel.
 */
export function keyCategoryValueIdsByPanel<T extends Categories>(
  config: MultiPanelOntologyFilterConfig,
  originalRows: Row<T>[]
): CategoryValueId[][] {
  const { filterOnKey, panels: panelConfigs } = config;
  return panelConfigs.reduce(
    (accum: CategoryValueId[][], panelConfig: CategoryFilterPanelConfig) => {
      // Determine the set of values for curated ontology panels.
      if (panelConfig.sourceKind === "CURATED") {
        const categoryValueIds = [
          ...listOntologyTreeIds(panelConfig.source),
        ].map((ontologyTermId) => buildInferredOntologyTermId(ontologyTermId));
        accum.push(categoryValueIds);
        return accum;
      }

      // Otherwise, build up the set of values for this panel from the original rows: only include explicit values.
      const categoryValueIds = originalRows.reduce(
        (explicitAccum, originalRow) => {
          // Adding type assertion here: curated categories are backed by x_ancestors fields, which are string arrays.
          // @ts-expect-error -- resolve mismatch between Categories and FilterKey.
          (originalRow.original[filterOnKey] as string[]).forEach(
            (categoryValueId) => {
              if (isExplicitOntologyTermId(categoryValueId as string)) {
                explicitAccum.add(categoryValueId);
              }
            }
          );
          return explicitAccum;
        },
        new Set<CategoryValueId>()
      );
      accum.push([...categoryValueIds]);
      return accum;
    },
    [] as CategoryValueId[][]
  );
}

/**
 * Create parent/child relationships between the given values. This is a not an ancestor/descendant relations but a UI-
 * specific parent/child relationship. For example, "blood, non-specific" is a descendant of both "hematopoietic system"
 * and "blood" but is only a child of "blood".
 * @param categoryValueId - Category value to build parent/child relationships for.
 * @param parentCategoryValueIdsByPanel - Category values of the panel that is an ancestor panel of the category
 * value's panel.
 * @param uiNodesByCategoryValueId - UI model to update with parent/child relationships.
 * @param descendants - Map of descendants keyed by ancestor.
 */
function linkParentsAndChildren(
  categoryValueId: CategoryValueId,
  parentCategoryValueIdsByPanel: CategoryValueId[][],
  uiNodesByCategoryValueId: Map<CategoryValueId, MultiPanelUINode>,
  descendants: OntologyDescendants
) {
  // Link parent and children between this panel and the parent panel if:
  // 1. Child value is a descendant of the parent value and,
  // 2. Child value is not a descendant of any children already added to parent.
  parentCategoryValueIdsByPanel.forEach((parentCategoryValueIds) => {
    parentCategoryValueIds.forEach((panelCategoryValueId) => {
      // Check if value is a descendant of panel value.
      const isDescendantOfPanelValue = isDescendant(
        categoryValueId,
        panelCategoryValueId,
        descendants
      );

      // If value isn't a descendant, there's no parent child relationship to link here.
      if (!isDescendantOfPanelValue) {
        return;
      }

      // Value is a descendant, check it is not "blocked". That is, the value is not a descendant of any children of
      // the panel value.
      // Note: the logic here possibly might need revisiting if more than three panels are added. For example, is it
      // possible that the children list is ever incomplete at this point?
      const panelCategoryValueUIChildren = getUIChildren(
        panelCategoryValueId,
        uiNodesByCategoryValueId
      );
      const isDescendantOfUIChild = panelCategoryValueUIChildren.some(
        (panelCategoryValueUIChild) => {
          // There are no descendants of explicit values; value can't be a descendant if child is explicit.
          if (isExplicitOntologyTermId(panelCategoryValueUIChild)) {
            return false;
          }

          // Check if value is a descendant of panel children.
          return isDescendant(
            categoryValueId,
            panelCategoryValueUIChild,
            descendants
          );
        }
      );

      // Only add value if it's not "blocked".
      if (!isDescendantOfUIChild) {
        linkParentAndChild(
          categoryValueId,
          panelCategoryValueId,
          uiNodesByCategoryValueId
        );
      }
    });
  });
}

/**
 * Create parent/child relationship between the given values.
 * @param categoryValueId - Category value to add as a child.
 * @param parentCategoryValueId - Category value to add as a parent.
 * @param uiNodesByCategoryValueId - UI model to update with parent/child relationships.
 */
function linkParentAndChild(
  categoryValueId: CategoryValueId,
  parentCategoryValueId: CategoryValueId,
  uiNodesByCategoryValueId: Map<CategoryValueId, MultiPanelUINode>
) {
  getUIChildren(parentCategoryValueId, uiNodesByCategoryValueId).push(
    categoryValueId
  );
  getUIParents(categoryValueId, uiNodesByCategoryValueId).push(
    parentCategoryValueId
  );
}

/**
 * Determine the partially selected category values from the given selected values.
 * @param selectedCategoryValueIds - All explicitly selected category value ids in all panels.
 * @param uiNodesByCategoryValueId - Map of categoryValueID to UI parents and children.
 * @return Array of category values that are partially selected.
 */
export function listPartiallySelectedCategoryValues(
  selectedCategoryValueIds: CategoryValueId[],
  uiNodesByCategoryValueId: Map<CategoryValueId, MultiPanelUINode>
): CategoryValueId[] {
  return selectedCategoryValueIds.reduce((accum, selectedCategoryValueId) => {
    if (
      isPartiallySelectedCategoryValue(
        selectedCategoryValueId,
        selectedCategoryValueIds,
        uiNodesByCategoryValueId
      )
    ) {
      accum.push(selectedCategoryValueId);
    }
    return accum;
  }, [] as CategoryValueId[]);
}

/**
 * List the selected values in the given parent panels.
 * @param parentPanelIndices - Array containing the indices of parent panels.
 * @param categoryFilterUIState - UI model of partially/selected values in the multi-panel category.
 * @param descendants - Map of descendants keyed by ancestor.
 * @returns Array of selected UI nodes in the parent panels.
 */
function listSelectedParentUINodes(
  parentPanelIndices: number[],
  categoryFilterUIState: MultiPanelCategoryFilterUIState,
  descendants: OntologyDescendants
): MultiPanelUINode[] {
  // Collect the selected values in the parent panels.
  const selectedParentUINodes = [
    ...categoryFilterUIState.uiNodesByCategoryValueId.values(),
  ].filter(
    (uiNode: MultiPanelUINode) =>
      categoryFilterUIState.selected.includes(uiNode.categoryValueId) &&
      parentPanelIndices.includes(uiNode.panelIndex)
  );

  // If any selected value has a descendant that is selected, remove the value from the set of selected values. For
  // example, if "digestive system" and "tongue" are both selected, we only want to include "tongue" as a selected value
  // otherwise the panel will show descendants for both "digestive system" and "tongue". Do not exclude descendants that
  // are in the same panel, for example "leukocyte" and "T cell": we always want to execute an "or" within a panel.
  return selectedParentUINodes.filter(
    (selectedParentUINode: MultiPanelUINode) =>
      !selectedParentUINodes.some((otherSelectedParentUINode) => {
        return (
          !isInSamePanel(
            otherSelectedParentUINode.categoryValueId,
            selectedParentUINode.categoryValueId,
            categoryFilterUIState.uiNodesByCategoryValueId
          ) &&
          isDescendant(
            otherSelectedParentUINode.categoryValueId,
            selectedParentUINode.categoryValueId,
            descendants
          )
        );
      })
  );
}

/**
 * Determine the indices of the panels that are parent to the given panel, if any.
 * @param panelIndex - Index of the panel to list parent panel indices for.
 * @param panelCount - Total number of panels.
 * @returns An array containing the indices of the parent panels.
 */
function listParentPanelIndices(
  panelIndex: number,
  panelCount: number
): number[] {
  return [...Array(panelCount).keys()].slice(0, panelIndex);
}

/**
 * List the parent nodes of the given panel.
 * @param panelIndex - Index of the panel to list parent nodes of.
 * @param categoryFilterUIState - UI model of partially/selected values in the multi-panel category.
 * @returns An array of parent UI nodes.
 */
function listPanelUINodes(
  panelIndex: number,
  categoryFilterUIState: MultiPanelCategoryFilterUIState
): MultiPanelUINode[] {
  return [...categoryFilterUIState.uiNodesByCategoryValueId.values()].filter(
    (uiNode: MultiPanelUINode) => uiNode.panelIndex === panelIndex
  );
}

/**
 * List the parent nodes of the given panels.
 * @param panelIndices - Indices of the panels to list parent nodes of.
 * @param categoryFilterUIState - UI model of partially/selected values in the multi-panel category.
 * @returns An array of parent UI nodes.
 */
function listPanelsUINodes(
  panelIndices: number[],
  categoryFilterUIState: MultiPanelCategoryFilterUIState
): MultiPanelUINode[] {
  return panelIndices.reduce(
    (accum: MultiPanelUINode[], panelIndex: number) => {
      accum.push(...listPanelUINodes(panelIndex, categoryFilterUIState));
      return accum;
    },
    [] as MultiPanelUINode[]
  );
}

/**
 * Build the list of selected values to be displayed as selected (blue) tags. Selected values must be "rolled up" to
 * fully-selected ancestors.
 * @param categoryValueViews - All view models in the panel.
 * @param categoryFilterUIState - UI model of partially/selected values in the category filter.
 * @returns An array of possibly "rolled-up" value views to display as selected tags.
 */
export function listMultiPanelSelectedViews(
  categoryValueViews: SelectCategoryValueView[],
  categoryFilterUIState: MultiPanelCategoryFilterUIState
): SelectCategoryValueView[] {
  const { selected, selectedPartial, uiNodesByCategoryValueId } =
    categoryFilterUIState;
  // Check if we can add any views to the select set, used to display selected tags.
  return categoryValueViews.filter(
    (categoryValueView: SelectCategoryValueView) =>
      categoryValueView.selected &&
      isSelectedViewTagVisible(
        categoryValueView.categoryValueId,
        selected,
        selectedPartial,
        uiNodesByCategoryValueId
      )
  );
}

/**
 * Handle select of multi-panel value: build and set next set of filters for this category. Track selected select value.
 * @param config - Configuration model of selected category.
 * @param selectedValue - Selected category value ID to use as selected value.
 * @param selectedLabel - Display value of selected value; required for tracking.
 * @param source - Location where select/remove was triggered, either filter menu or selected tag.
 * @param multiPanelUIState - Current set of category values that the user has selected on the UI for all multi-
 * panel category filters.
 */
export function onFilterMultiPanelCategory(
  config: MultiPanelOntologyFilterConfig,
  selectedValue: CategoryValueId,
  selectedLabel: string,
  source: ON_FILTER_SOURCE,
  multiPanelUIState: MultiPanelUIState
): [MultiPanelUIState, CategoryValueId[]] {
  // Model the selected values for this category filter in a react-table filters format.
  const currentFilters = buildMultiPanelCurrentFilters(multiPanelUIState);

  // Grab the UI model backing this category.
  const { categoryFilterId } = config;
  const categoryFilterUIState = multiPanelUIState.get(categoryFilterId);
  if (!categoryFilterUIState) {
    return [multiPanelUIState, []]; // Error state
  }

  // Track selected category and value.
  trackMultiPanelCategoryValueSelected(
    config,
    selectedValue,
    selectedLabel,
    currentFilters,
    categoryFilterUIState
  );

  // Determine selected set of values for this category filter based on the current selected values; toggle current
  // selected values. If this is a remove action from a tag, remove descendants of fully (not partially) selected
  // ancestors.
  const nextFilters = buildNextMultiPanelCategoryFilters(
    config,
    selectedValue,
    source,
    categoryFilterUIState,
    currentFilters
  );

  // Determine the selected values which are now partially selected, if any. For example, if "digestive system" and
  // "tongue" are both selected in the UI, internally we record "digestive system" as partially selected and "tongue"
  // as selected.
  const { uiNodesByCategoryValueId } = categoryFilterUIState;
  const selectedPartial = listPartiallySelectedCategoryValues(
    nextFilters,
    uiNodesByCategoryValueId
  );

  // Determine the selected values to pass to react-table to filter the rows by applying restrictions across selected
  // filter values. For example, even if "digestive system" and "tongue" are both selected in the UI, we only want to
  // pass the more restrictive "tongue" to react-table as we only want to show rows that match "tongue".
  const overriddenNextFilters = overrideSelectedParents(
    nextFilters,
    config.descendants,
    uiNodesByCategoryValueId
  );

  // Update internal multi-panel filter state with the updated, toggled set of selected filters for this category as
  // well as the partially (overridden) selected values.
  const nextMultiPanelUIState = new Map(multiPanelUIState);
  nextMultiPanelUIState.set(categoryFilterId, {
    ...categoryFilterUIState,
    selected: nextFilters,
    selectedPartial,
  });

  return [nextMultiPanelUIState, overriddenNextFilters];
}

/**
 * Handle remove of a category value from a multi-panel from a selected tag. If a category value is selected (and
 * not partially selected) then remove any selected descendant terms. Also allow tags to "unravel" such as if
 * hematopoietic system and blood non-specific are selected and blood non-specific is removed, hematopoietic system
 * remains selected.
 * @param removedCategoryValueId - The selected value to remove.
 * @param selectedCategoryValueIds - The set of selected values after the selected value has been toggled.
 * @param categoryFilterUIState - UI model of selected values for the multi-panel category filter.
 * @param descendants - Map of descendants keyed by ancestor.
 * @returns An array of selected values with descendant-related removal rules applied.
 */
export function onRemoveMultiPanelCategoryValueTag(
  removedCategoryValueId: CategoryValueId,
  selectedCategoryValueIds: CategoryValueId[],
  categoryFilterUIState: MultiPanelCategoryFilterUIState,
  descendants: OntologyDescendants
): CategoryValueId[] {
  const filteredSelectedCategoryFilters = new Set<CategoryValueId>(
    selectedCategoryValueIds
  );

  // Grab all selected descendants of the removed value.
  const selectedDescendants = selectedCategoryValueIds.filter(
    (selectedCategoryValueId) =>
      isDescendant(
        selectedCategoryValueId,
        removedCategoryValueId,
        descendants
      ) &&
      !isInSamePanel(
        selectedCategoryValueId,
        removedCategoryValueId,
        categoryFilterUIState.uiNodesByCategoryValueId
      )
  );

  // Remove the selected descendant unless the descendant has another parent that is partially/selected.
  const { selected, selectedPartial, uiNodesByCategoryValueId } =
    categoryFilterUIState;

  selectedDescendants.forEach((selectedDescendant) => {
    const uiParents = getUIParents(
      selectedDescendant,
      uiNodesByCategoryValueId
    );

    // If there are no parents, remove the descendant from the selected list.
    if (!uiParents || !uiParents.length) {
      filteredSelectedCategoryFilters.delete(selectedDescendant);
    } else {
      // Don't delete selected descendant if any parent is partially/selected
      const isAnyParentSelected = uiParents
        .filter(
          (uiParent) =>
            uiParent !== removedCategoryValueId &&
            !selectedDescendants.includes(uiParent)
        )
        .some(
          (uiParent) =>
            selected.includes(uiParent) || selectedPartial.includes(uiParent)
        );
      if (!isAnyParentSelected) {
        filteredSelectedCategoryFilters.delete(selectedDescendant);
      }
    }
  });

  return [...filteredSelectedCategoryFilters];
}

/**
 * Determine the selected values to pass to react-table to filter the rows by applying restrictions across selected
 * filter values. For example, even if "digestive system" and "tongue" are both selected in the UI, we only want to
 * pass the more restrictive "tongue" to react-table as we only want to show rows that match "tongue".
 * @param selectedValues - Set of currently selected values.
 * @param descendants - Map of descendants keyed by ancestor.
 * @param uiNodesByCategoryValueId - Map of categoryValueID to UI parents and children.
 * @returns The "effective" set of selected terms to pass to react-table.
 */
export function overrideSelectedParents(
  selectedValues: CategoryValueId[],
  descendants: OntologyDescendants,
  uiNodesByCategoryValueId: Map<CategoryValueId, MultiPanelUINode>
): CategoryValueId[] {
  const selectedOntologyTermIds = selectedValues.map((selectedValue) =>
    removeOntologyTermIdPrefix(selectedValue)
  );
  return selectedValues.reduce((accum, selectedValue: CategoryValueId) => {
    const selectedOntologyId = removeOntologyTermIdPrefix(selectedValue);

    // If the selected value is an explicit value, always include it in the selected set of values.
    if (isExplicitOntologyTermId(selectedValue)) {
      accum.push(selectedValue);
      return accum;
    }

    // Handle self case. Check if the corresponding explicit value is selected. For example, "eye" organ (which is
    // inferred) is a parent of "eye" tissue (which is explicit). In this case, the explicit "eye" tissue takes
    // precedence over the "eye" organ.
    const explicitSelectedValue =
      buildExplicitOntologyTermId(selectedOntologyId);
    if (selectedValues.includes(explicitSelectedValue)) {
      return accum;
    }

    // Otherwise, if any descendant of the selected value is selected, do not include the selected value in the set of
    // selected values. Ignore descendants that are in the same panel as the selected value (for example, "leukocyte"
    // and "T cell"); we want to always execute an "or" within a panel so both ancestor/descendant selected values must
    // be included in this case.
    const descendantsOfSelected = (
      descendants[selectedOntologyId] ?? []
    ).filter((descendant) => {
      return !isInSamePanel(
        selectedValue,
        // Explicit values are handled above; we only want to check for inferred terms.
        buildInferredOntologyTermId(descendant),
        uiNodesByCategoryValueId
      );
    });

    // Check if any descendant is selected.
    const isAnyDescendantSelected = descendantsOfSelected.some((descendant) =>
      selectedOntologyTermIds.includes(descendant)
    );
    if (!isAnyDescendantSelected) {
      accum.push(selectedValue);
    }

    return accum;
  }, [] as CategoryValueId[]);
}

/**
 * Returns the prefix of an ontology term ID that has been prefixed with an inferred or explicit identifier.
 * @param prefixedOntologyTermId - Ontology term ID with either an inferred or explicit prefix.
 * @return The prefix of the given ontology term ID that has been marked as inferred or explicit.
 */
function removeOntologyTermId(prefixedOntologyTermId: string): OrFilterPrefix {
  const [prefix] = splitOntologyTermIdAndPrefix(prefixedOntologyTermId);
  return prefix;
}

/**
 * Returns the core ontology term ID that has been prefixed with an inferred or explicit identifier.
 * @param prefixedOntologyTermId - Ontology term ID with either an inferred or explicit prefix.
 * @return Core ontology term ID.
 */
export function removeOntologyTermIdPrefix(
  prefixedOntologyTermId: string
): string {
  const [, ontologyTermId] = splitOntologyTermIdAndPrefix(
    prefixedOntologyTermId
  );
  return ontologyTermId;
}

/**
 * Returns the prefix and core ontology term ID of an ontology term ID that has been prefixed with an inferred or
 * explicit identifier.
 * @param prefixedOntologyTermId - Ontology term ID with either an inferred or explicit prefix.
 * @return An array containing the prefix and the core ontology term ID.
 */
function splitOntologyTermIdAndPrefix(
  prefixedOntologyTermId: string
): [OrFilterPrefix, string] {
  return prefixedOntologyTermId.split(/:(.*)/s) as [OrFilterPrefix, string];
}

/**
 * Track select of the given multi-panel category and category value.
 * @param config - Configuration model of selected category.
 * @param categoryValueId - Selected category value (e.g. "hematopoietic system").
 * @param selectedLabel - Display value of selected value; used as event payload.
 * @param filters - Current set of selected category values.
 * @param categoryFilterUIState - UI model of partially/selected values in the category filter.
 */
function trackMultiPanelCategoryValueSelected<T extends Categories>(
  config: MultiPanelOntologyFilterConfig,
  categoryValueId: CategoryValueId,
  selectedLabel: string,
  filters: Filters<T>,
  categoryFilterUIState: MultiPanelCategoryFilterUIState
) {
  const { categoryFilterId, panels: panelConfigs } = config;

  // Determine the panel the value was selected from.
  const uiNode =
    categoryFilterUIState.uiNodesByCategoryValueId.get(categoryValueId);
  if (!uiNode) {
    console.log(`UI node not found for ${categoryValueId}.`);
    return; // Error state.
  }

  // Get the event for the panel.
  const { analyticsEvent, analyticsPayloadKey } =
    panelConfigs[uiNode.panelIndex];

  // No tracking if event and payload key aren't specified on multi-panel category configs.
  if (!analyticsEvent || !analyticsPayloadKey) {
    return;
  }

  // Only track the select (and not deselect) of category value.
  const categoryFilters = new Set(
    getCategoryFilter(categoryFilterId, filters)?.value as CategoryValueId
  );
  if (!categoryFilters.has(categoryValueId)) {
    // Build up payload for tracking event and send.
    track(analyticsEvent, { [analyticsPayloadKey]: selectedLabel });
  }
}
