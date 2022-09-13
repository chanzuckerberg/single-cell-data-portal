import { Filters } from "react-table";
import { track } from "src/common/analytics";
import { getCategoryFilter } from "src/common/hooks/useCategoryFilter/common/utils";
import { TOOLTIP_CATEGORY_DISABLED } from "src/components/common/Filter/common/constants";
import {
  Categories,
  CategorySet,
  CategoryValueId,
  CATEGORY_FILTER_ID,
  CuratedOntologyCategoryFilterConfig,
  FilterState,
  KeyedSelectCategoryValue,
  OntologyCategoryTreeNodeView,
  OntologyCategoryTreeView,
  OntologyCategoryView,
  OntologyNode,
  OntologyTermSet,
  ONTOLOGY_VIEW_KEY,
  ONTOLOGY_VIEW_LABEL,
  ORGANISM,
} from "src/components/common/Filter/common/entities";
import {
  findOntologyNodeById,
  findOntologyParentNode,
  getOntologySpeciesKey,
} from "src/components/common/Filter/common/utils";

/**
 * Utils specific to curated ontology category filters.
 */

/**
 * Add all descendents to the set of selected values. Do not add IDs of descendents that are specified in the
 * ontology view but not present in the full set of data as this will filter the data by values that do not exist
 * (resulting in an empty set).
 * @param ontologyNodes - Nodes to add to selected values.
 * @param selectedCategoryValues - The current set of selected values.
 * @param categoryKeyValues - Original, full set of values for this category.
 */
function addOntologyDescendents(
  ontologyNodes: OntologyNode[],
  selectedCategoryValues: Set<CategoryValueId>,
  categoryKeyValues: Set<CategoryValueId>
) {
  // For each node, add self and possibly children.
  ontologyNodes.forEach((ontologyNode) => {
    const { ontology_term_id: ontologyId } = ontologyNode;
    // Only add self if memeber of original set of values for this category.
    if (categoryKeyValues.has(ontologyId)) {
      selectedCategoryValues.add(ontologyId);
    }
    // Add children, if any.
    if (ontologyNode.children) {
      addOntologyDescendents(
        ontologyNode.children,
        selectedCategoryValues,
        categoryKeyValues
      );
    }
  });
}

/**
 * Build view model of a curated ontology category such as development stage.
 * @param config - Config model of a curated ontology category.
 * @param categoryValueByValue - Internal filter model of ontology category.
 * @param filterState - Categories, category value and their counts with the current filter applied. Required when
 * @param ontologyTermLabelsById - Set of ontology term labels keyed by term ID, used to determine labels for ontology.
 * checking enabled state of view that is dependent on the state of another category.
 * @returns Ontology view model.
 */
export function buildCuratedOntologyCategoryView(
  config: CuratedOntologyCategoryFilterConfig,
  categoryValueByValue: KeyedSelectCategoryValue,
  filterState: FilterState,
  ontologyTermLabelsById: Map<string, string>
): OntologyCategoryView {
  const {
    categoryFilterId,
    isLabelVisible,
    isSearchable,
    isZerosVisible,
    label,
    source,
  } = config;

  // Build tree view models (e.g. individual tree structures for displaying different ontologies (e.g. human vs mouse
  // vs other for development stage, or just tissues for tissue).
  const treeViews = Object.keys(source).reduce(
    (accum, ontologyViewKey: string) => {
      const ontologyNodes = source[ontologyViewKey as ONTOLOGY_VIEW_KEY];
      if (!ontologyNodes) {
        return accum; // Error state - ignore species view.
      }

      // Handle special cases where species is to be excluded.
      if (
        categoryFilterId === CATEGORY_FILTER_ID.DEVELOPMENT_STAGE &&
        !isDevelopmentStageSpeciesVisible(
          filterState,
          ontologyViewKey as ONTOLOGY_VIEW_KEY
        )
      ) {
        return accum;
      }

      // Build view model for each node.
      const childrenViews = ontologyNodes.map((ontologyNode) =>
        buildCuratedOntologyCategoryValueView(
          ontologyNode,
          categoryValueByValue,
          ontologyTermLabelsById
        )
      );

      // Calculate partial selected states.
      childrenViews.forEach((childView) =>
        markOntologySelectedPartialViews(childView)
      );

      // Determine the set of selected values that should be included as tags. Adding this as a flat set to prevent
      // UI code from having to recursively determine selected set when rendering tags.
      const selectedViews = childrenViews.reduce((accum, childView) => {
        listOntologySelectedViews(childView, accum);
        return accum;
      }, new Set<OntologyCategoryTreeNodeView>());

      const viewLabel =
        ONTOLOGY_VIEW_LABEL[
          ontologyViewKey as keyof typeof ONTOLOGY_VIEW_LABEL
        ];

      accum.push({
        children: childrenViews,
        label: isLabelVisible ? viewLabel : undefined,
        selectedViews: [...selectedViews.values()],
      });

      return accum;
    },
    [] as OntologyCategoryTreeView[]
  );

  // Build up the ontology category view model.
  const ontologyView: OntologyCategoryView = {
    categoryFilterId,
    isSearchable,
    isZerosVisible,
    label,
    views: treeViews,
  };

  // Check if ontology category is disabled.
  if (isOntologyCategoryViewDisabled(ontologyView)) {
    ontologyView.isDisabled = true;
    ontologyView.tooltip = TOOLTIP_CATEGORY_DISABLED;
  }

  return ontologyView;
}

/**
 * Build view model of node of ontology tree to be displayed as a value in an ontology menu. This is specific to curated
 * ontology category filters such as development stage.
 * @param ontologyNode - Ontology node to build view model for.
 * @param categoryValueByValue - Internal filter model of ontology category.
 * @param ontologyTermLabelsById - Set of ontology term labels keyed by term ID, used to determine labels for ontology.
 * @returns Ontology view model.
 */
export function buildCuratedOntologyCategoryValueView(
  ontologyNode: OntologyNode,
  categoryValueByValue: KeyedSelectCategoryValue,
  ontologyTermLabelsById: Map<string, string>
): OntologyCategoryTreeNodeView {
  const { ontology_term_id: categoryValueId } = ontologyNode;
  const categoryValue = categoryValueByValue.get(categoryValueId);

  // If there's no corresponding category for this node, create a basic model with 0 count and unselected.
  if (!categoryValue) {
    return {
      categoryValueId,
      count: 0,
      label: ontologyNode.label,
      selected: false,
      selectedPartial: false,
    };
  }

  // Build up base view model.
  const view = {
    categoryValueId,
    count: categoryValue.count,
    label: ontologyTermLabelsById.get(categoryValueId) ?? categoryValueId,
  };

  // If this ontology node is a leaf, add its selected value and return.
  if (!ontologyNode.children) {
    return {
      ...view,
      selected: categoryValue.selected,
      selectedPartial: false,
    };
  }

  // Otherwise, build view models for child nodes.
  const children = ontologyNode.children.map((childNode) =>
    buildCuratedOntologyCategoryValueView(
      childNode,
      categoryValueByValue,
      ontologyTermLabelsById
    )
  );

  // Build up view model for this node, including children.
  return {
    ...view,
    children,
    selected: categoryValue.selected,
    selectedPartial: false,
  };
}

/**
 * Build updated set of selected filters for the given ontology tree category and the selected category value.
 * @param categoryFilterId - Category ID (i.e. "development stage") of selected category value.
 * @param selectedValue - Selected category value ID (e.g. "HsapDv:0000003") to update selected state of.
 * @param filters - Current set of selected category values.
 * @param categoryKeyValues - Original, full set of values for this category.
 * @param mask - View model of ontology for this category.
 * @returns Array of selected category values for the given category.
 */
export function buildNextOntologyCategoryFilters<T extends Categories>(
  categoryFilterId: CATEGORY_FILTER_ID,
  selectedValue: CategoryValueId,
  filters: Filters<T>,
  categoryKeyValues: Set<CategoryValueId>,
  mask: OntologyTermSet
): CategoryValueId[] {
  // Grab the current selected values for the category.
  const categoryFilters = new Set(
    getCategoryFilter(categoryFilterId, filters)?.value as CategoryValueId[]
  );

  // Find the selected and parent node, if any, for the selected value.
  const ontologySpeciesKey = getOntologySpeciesKey(selectedValue);
  const ontologyRootNodes = mask[ontologySpeciesKey];
  if (!ontologyRootNodes) {
    return [...categoryFilters.values()]; // Error state - ontology does not exist.
  }
  const selectedOntologyNode = findOntologyNodeById(
    ontologyRootNodes,
    selectedValue
  );
  if (!selectedOntologyNode) {
    return [...categoryFilters.values()]; // Error state - ontology node with given ID does not exist.
  }
  const parentNode = findOntologyParentNode(
    ontologyRootNodes,
    selectedOntologyNode
  );

  // Toggle selected state of selected category value.
  if (categoryFilters.has(selectedValue)) {
    // Selected value is already in the set of selected values, remove it.
    categoryFilters.delete(selectedValue);

    // Also remove any descendents from the selected set.
    if (selectedOntologyNode.children) {
      removeOntologyDescendents(selectedOntologyNode.children, categoryFilters);
    }

    // Reevaluate selected state of parent. If all children were previously selected and now only some are selected,
    // then parent should no longer be selected.
    if (parentNode) {
      handleOntologyChildRemoved(
        ontologyRootNodes,
        parentNode,
        categoryFilters
      );
    }
  } else {
    // Add selected value to selected set.
    categoryFilters.add(selectedValue);

    // Add all descendents of selected value, if any.
    if (selectedOntologyNode.children) {
      addOntologyDescendents(
        selectedOntologyNode.children,
        categoryFilters,
        categoryKeyValues
      );
    }

    // Reevaluate selected state of parent. If all children are now selected, then parent should also be selected.
    if (parentNode) {
      handleOntologyChildAdded(
        ontologyRootNodes,
        parentNode,
        categoryFilters,
        categoryKeyValues
      );
    }
  }

  return [...categoryFilters.values()];
}

/**
 * Reevaluate selected state of parent. It's possible all children of this parent are now selected; if so, add parent
 * to set of selected values and then reevaluate the parent's parent selected state.
 * @param ontologyRootNodes - Top-level nodes in ontology tree.
 * @param ontologyNode - Node to reevaluate selected state of.
 * @param selectedCategoryValues - The current set of selected values.
 * @param categoryKeyValues - Original, full set of values for this category.
 */
function handleOntologyChildAdded(
  ontologyRootNodes: OntologyNode[],
  ontologyNode: OntologyNode,
  selectedCategoryValues: Set<CategoryValueId>,
  categoryKeyValues: Set<CategoryValueId>
) {
  // Check if all children of the parent node are selected, only including children values that are present in the
  // original result set.
  const isEveryChildSelected = isEveryChildNodeSelected(
    ontologyNode,
    selectedCategoryValues,
    categoryKeyValues
  );

  // Add parent if all children are selected.
  if (isEveryChildSelected) {
    selectedCategoryValues.add(ontologyNode.ontology_term_id);

    // Node's parent, if any, must now be reevaluated.
    const parentNode = findOntologyParentNode(ontologyRootNodes, ontologyNode);
    if (parentNode) {
      handleOntologyChildAdded(
        ontologyRootNodes,
        parentNode,
        selectedCategoryValues,
        categoryKeyValues
      );
    }
  }
}

/**
 * Reevaluate selected state of node. If a child has been removed, then a parent can no longer be considered selected;
 * remove node from selected set and reevaluate the parent's parent selected state.
 * @param ontologyRootNodes - Top-level nodes in ontology tree.
 * @param ontologyNode - Node to reevaluate selected state of.
 * @param selectedCategoryValues - The current set of selected values.
 */
function handleOntologyChildRemoved(
  ontologyRootNodes: OntologyNode[],
  ontologyNode: OntologyNode,
  selectedCategoryValues: Set<CategoryValueId>
) {
  const { ontology_term_id: parentId } = ontologyNode;
  if (selectedCategoryValues.has(parentId)) {
    // Parent is currently selected, remove it from the selected set.
    selectedCategoryValues.delete(parentId);

    // Node's parent, if any, must now be reevaluated.
    const parentNode = findOntologyParentNode(ontologyRootNodes, ontologyNode);
    if (parentNode) {
      handleOntologyChildRemoved(
        ontologyRootNodes,
        parentNode,
        selectedCategoryValues
      );
    }
  }
}

/**
 * Development stage species is only visible if:
 * 1. There are no selected organisms or,
 * 2. The given species is selected.
 * @param filterState - Categories, category value and their counts with the current filter applied. Required to
 * determine if development stfage species should be visible.
 * @param speciesKey - The species to check if a corresponding organism has been selected for.
 * @returns True if given species is to be displayed.
 */
function isDevelopmentStageSpeciesVisible(
  filterState: FilterState,
  speciesKey: ONTOLOGY_VIEW_KEY
) {
  // Find the current selected values for organism.
  const organismCategoryValues = filterState[
    CATEGORY_FILTER_ID.ORGANISM
  ] as KeyedSelectCategoryValue;
  const selectedOrganisms = [...organismCategoryValues.values()]
    .filter((selectCategoryValue) => selectCategoryValue.selected)
    .map((selectCategoryValue) => selectCategoryValue.categoryValueId);

  // If no organisms are selected, all species can be displayed.
  if (selectedOrganisms.length === 0) {
    return true;
  }

  // Otherwise this species is only visible if it's selected.
  if (speciesKey === ONTOLOGY_VIEW_KEY.HsapDv) {
    return selectedOrganisms.includes(ORGANISM.HOMO_SAPIENS);
  }
  if (speciesKey === ONTOLOGY_VIEW_KEY.MmusDv) {
    return selectedOrganisms.includes(ORGANISM.MUS_MUSCULUS);
  }
  // Check the "other" case where any species other than human and mouse must be selected.
  return (
    selectedOrganisms.filter(
      (organism) =>
        organism !== ORGANISM.HOMO_SAPIENS && organism !== ORGANISM.MUS_MUSCULUS
    ).length > 0
  );
}

/*
 * Check if all children of the parent node are selected, only including children values that are present in the
 * original result set.
 * @param ontologyNode - Node to check selected state of children.
 * @param selectedCategoryValues - The current set of selected values.
 * @param categoryKeyValues - Original, full set of values for this category.
 */
function isEveryChildNodeSelected(
  ontologyNode: OntologyNode,
  selectedCategoryValues: Set<CategoryValueId>,
  categoryKeyValues: Set<CategoryValueId>
): boolean {
  // Check if all children of the parent node are selected, only including children values that are present in the
  // original result set.
  const childrenIds =
    ontologyNode.children
      ?.map((child) => child.ontology_term_id)
      .filter((childId) => categoryKeyValues.has(childId)) ?? [];
  return Boolean(
    childrenIds.length &&
      childrenIds.every((childId) => selectedCategoryValues.has(childId))
  );
}

/**
 * Returns true if ontology category is disabled, that is, there are no species ontology trees or species
 * ontology trees with values that have a count.
 * @param categoryView - Ontology category view to check enabled/disabled state of.
 * @returns True when there are no species or the count of children for each species is 0.
 */
function isOntologyCategoryViewDisabled(
  categoryView: OntologyCategoryView
): boolean {
  const { views } = categoryView;
  if (!views || views.length === 0) {
    return true;
  }
  return !views.some((s) => s.children.some((child) => child.count > 0));
}

/**
 * Determine the set of selected values (views) that form the basis of the set of selected tags. Adding this as a flat
 * set to prevent UI code from having to recursively determine selected set.
 * @param view - View model of ontology category value.
 * @param selectedSet - Selected set of cateogry value keys.
 */
function listOntologySelectedViews(
  view: OntologyCategoryTreeNodeView,
  selectedSet: Set<OntologyCategoryTreeNodeView>
) {
  // If the view is selected it can be included in the selected set of tags.
  if (view.selected) {
    selectedSet.add(view);
    return;
  }

  // If the view is partially selected, check its children to see which ones can be included in the set of selected tags.
  if (view.selectedPartial && view.children) {
    view.children.forEach((childView) =>
      listOntologySelectedViews(childView, selectedSet)
    );
  }

  return selectedSet;
}

/**
 * Handle select of ontology value: build and set next set of filters for this category. Track selected ontology value.
 * @param config - Configuration model of selected category.
 * @param selectedValue - Selected category value ID to use as selected value.
 * @param filters - Current set of selected category values (values) or ranges keyed by category (id).
 * @param categorySet - Original, unfiltered sets of category values keyed by their category.
 * @returns The selected values with the latest selection applied.
 */
export function onFilterCuratedOntologyCategory<T extends Categories>(
  config: CuratedOntologyCategoryFilterConfig,
  selectedValue: CategoryValueId,
  filters: Filters<T>,
  categorySet: CategorySet
) {
  const { categoryFilterId, source } = config;

  // Track selected category and value.
  trackCuratedOntologyCategoryValueSelected(config, selectedValue, filters);

  // Build and set next set of filters for this category.
  return buildNextOntologyCategoryFilters(
    categoryFilterId,
    selectedValue,
    filters,
    categorySet[categoryFilterId] as Set<CategoryValueId>,
    source
  );
}

/**
 * Update partially selected state of views. If parent isn't selected but some of its children are selected or partially
 * selected, mark parent as partially selected
 * @param view - View model of ontology category value.
 */
function markOntologySelectedPartialViews(view: OntologyCategoryTreeNodeView) {
  // Mark any children as partially selected.
  view.children?.forEach((childView) => {
    markOntologySelectedPartialViews(childView);
  });

  // If this view isn't selected, check if any of it's children are either selected or partially selected.
  if (view.selected) {
    return;
  }
  view.selectedPartial =
    view.children?.some(
      (childView) => childView.selected || childView.selectedPartial
    ) ?? false;
}

/**
 * Remove all descendents from the set of selected values.
 * @param ontologyNodes - Nodes to remove from selected set of values.
 * @param selectedCategoryValues - The current set of selected values.
 */
function removeOntologyDescendents(
  ontologyNodes: OntologyNode[],
  selectedCategoryValues: Set<CategoryValueId>
) {
  // For each node, remove self and children.
  ontologyNodes.forEach((ontologyNode) => {
    selectedCategoryValues.delete(ontologyNode.ontology_term_id);
    if (ontologyNode.children) {
      removeOntologyDescendents(ontologyNode.children, selectedCategoryValues);
    }
  });
}

/**
 * Track select of the given ontology category and category value.
 * @param config - Configuration model of selected category.
 * @param categoryValueKey - Selected category value key (e.g. "HsapDv:0000003").
 * @param filters - Current set of selected category values (values) or ranges keyed by category (id).
 */
function trackCuratedOntologyCategoryValueSelected<T extends Categories>(
  config: CuratedOntologyCategoryFilterConfig,
  categoryValueKey: CategoryValueId,
  filters: Filters<T>
) {
  const { analyticsEvent, categoryFilterId, source } = config;

  // No tracking if event isn't specified on category config.
  if (!analyticsEvent) {
    return;
  }

  // Only track the select (and not deselect) of category value.
  const categoryFilters = new Set(
    getCategoryFilter(categoryFilterId, filters)?.value as CategoryValueId[]
  );
  if (!categoryFilters.has(categoryValueKey)) {
    // Grab the analytics event for this category.

    // Find the node for the selected value.
    const ontologySpeciesKey = getOntologySpeciesKey(categoryValueKey);
    const ontologyRootNodes = source[ontologySpeciesKey];
    if (!ontologyRootNodes) {
      return; // Error state - ontology does not exist.
    }
    const selectedOntologyNode = findOntologyNodeById(
      ontologyRootNodes,
      categoryValueKey
    );
    if (!selectedOntologyNode) {
      return; // Error state - ontology node with given ID does not exist.
    }

    // Build up payload for tracking event and send.
    track(analyticsEvent, {
      label: selectedOntologyNode.label,
      ontologyTermId: selectedOntologyNode.ontology_term_id,
    });
  }
}
