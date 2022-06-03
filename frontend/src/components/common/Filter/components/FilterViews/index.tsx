import {
  CATEGORY_KEY,
  OnFilterFn,
  OntologyCategoryTreeView,
  OnUpdateSearchValueFn,
} from "src/components/common/Filter/common/entities";
import FilterView, {
  MAX_DISPLAYABLE_LIST_ITEMS,
} from "src/components/common/Filter/components/FilterViews/components/FilterView";
import {
  VIEW_LIST_ITEM_HEIGHT,
  VIEW_LIST_SUBHEADER_HEIGHT,
} from "src/components/common/Filter/components/FilterViews/components/FilterView/style";
import { ViewsMenu } from "src/components/common/Filter/components/FilterViews/style";

export const enum CATEGORY_VIEWS_QUANTIFIER {
  NON_SINGLETON = "NON_SINGLETON",
  SINGLETON = "SINGLETON",
}

interface Props {
  categoryKey: CATEGORY_KEY;
  isSearchable: boolean;
  isZerosVisible: boolean;
  onFilter: OnFilterFn;
  onUpdateSearchValue: OnUpdateSearchValueFn;
  views: OntologyCategoryTreeView[];
}

export default function FilterViews({
  categoryKey,
  isSearchable,
  isZerosVisible,
  onFilter,
  onUpdateSearchValue,
  views,
}: Props): JSX.Element {
  const viewsToDisplay = views.filter(
    (s) => s.children && s.children.filter((child) => child.count).length > 0
  );
  const viewsQuantifier =
    viewsToDisplay.length === 1
      ? CATEGORY_VIEWS_QUANTIFIER.SINGLETON
      : CATEGORY_VIEWS_QUANTIFIER.NON_SINGLETON;
  const viewListMaxHeight = calculateViewListMaxHeight(
    viewsToDisplay,
    viewsQuantifier
  );
  return (
    <ViewsMenu>
      {viewsToDisplay.map(({ label, children }, i) => (
        <FilterView
          categoryKey={categoryKey}
          key={`${label || "view"}-${i}`}
          label={label}
          isSearchable={isSearchable}
          isZerosVisible={isZerosVisible}
          onFilter={onFilter}
          onUpdateSearchValue={onUpdateSearchValue}
          showViewDivider={i !== 0}
          values={children}
          viewListMaxHeight={viewListMaxHeight}
        />
      ))}
    </ViewsMenu>
  );
}

/**
 * Returns the maximum possible pixel height of any view list.
 * @param views - Tree view model of ontology view.
 * @param viewsQuantifier - Singleton or non-singleton view.
 * @return maximum possible view list height in pixels.
 */
function calculateViewListMaxHeight(
  views: OntologyCategoryTreeView[],
  viewsQuantifier: CATEGORY_VIEWS_QUANTIFIER
): number {
  return views.reduce((acc, { label }) => {
    const listMaxHeight =
      (MAX_DISPLAYABLE_LIST_ITEMS[viewsQuantifier] + 0.5) *
        VIEW_LIST_ITEM_HEIGHT +
      VIEW_LIST_SUBHEADER_HEIGHT * (label ? 1 : 0);
    return Math.max(acc, listMaxHeight);
  }, 0);
}
