import {
  CATEGORY_KEY,
  OnFilterFn,
  OntologyCategoryTreeNodeView,
  OntologyCategoryTreeView,
} from "src/components/common/Filter/common/entities";
import FilterView from "src/components/common/Filter/components/FilterViews/components/FilterView";
import { MAX_DISPLAYABLE_LIST_ITEMS } from "src/components/common/Filter/components/FilterViews/components/FilterView/style";
import { ViewsMenu } from "src/components/common/Filter/components/FilterViews/style";

interface Props {
  categoryKey: CATEGORY_KEY;
  onFilter: OnFilterFn;
  views: OntologyCategoryTreeView[];
}

export default function FilterViews({
  categoryKey,
  onFilter,
  views,
}: Props): JSX.Element {
  const viewsToDisplay = views.filter(
    (s) => s.children && s.children.filter((child) => child.count).length > 0
  );
  return (
    <ViewsMenu>
      {viewsToDisplay.map(({ label, children }, i) => (
        <FilterView
          categoryKey={categoryKey}
          key={`${label || "view"}-${i}`}
          label={label}
          onFilter={onFilter}
          scrollable={countViews(children) > MAX_DISPLAYABLE_LIST_ITEMS}
          showViewDivider={i !== 0}
          values={children}
        />
      ))}
    </ViewsMenu>
  );
}

/**
 * Returns a count of all views in the given view tree.
 * @param children
 * @param increment
 * @returns number of all views in the given view tree.
 */
function countViews(
  children: OntologyCategoryTreeNodeView[],
  increment = 0
): number {
  return children.reduce((acc, child) => {
    acc++;
    if (child.children) {
      return countViews(child.children, acc);
    }
    return acc;
  }, increment);
}
