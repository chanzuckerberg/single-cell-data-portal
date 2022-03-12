import {
  CATEGORY_KEY,
  OnFilterFn,
  OntologyCategorySpeciesView,
  OntologyCategoryValueView,
} from "src/components/common/Filter/common/entities";
import FilterPanel from "src/components/common/Filter/components/FilterMultiPanel/components/FilterPanel";
import { MAX_DISPLAYABLE_LIST_ITEMS } from "src/components/common/Filter/components/FilterMultiPanel/components/FilterPanel/style";
import { MultiPanelSelector } from "src/components/common/Filter/components/FilterMultiPanel/style";

interface Props {
  categoryKey: CATEGORY_KEY;
  onFilter: OnFilterFn;
  species: OntologyCategorySpeciesView[];
}

export default function FilterMultiPanel({
  categoryKey,
  onFilter,
  species,
}: Props): JSX.Element {
  const speciesToDisplay = species.filter(
    (s) => s.children && s.children.length > 0
  );
  return (
    <MultiPanelSelector>
      {speciesToDisplay.map(({ label, children }, i) => (
        <FilterPanel
          categoryKey={categoryKey}
          key={label}
          label={label}
          onFilter={onFilter}
          scrollable={countViews(children) > MAX_DISPLAYABLE_LIST_ITEMS}
          showPanelDivider={i !== 0}
          values={children}
        />
      ))}
    </MultiPanelSelector>
  );
}

/**
 * Returns a count of all views in the given view tree.
 * @param children
 * @param increment
 * @returns number of species-specific ontology tree nodes.
 */
function countViews(
  children: OntologyCategoryValueView[],
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
