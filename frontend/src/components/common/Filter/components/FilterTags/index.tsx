import { Tag } from "@blueprintjs/core";
import {
  CategoryValueView,
  OnFilterFn,
} from "src/common/hooks/useCategoryFilter";
import { CATEGORY_KEY } from "src/components/common/Filter/common/entities";
import { SelectedTags } from "./style";

interface Props {
  categoryKey: CATEGORY_KEY;
  onFilter: OnFilterFn;
  selectedValues: CategoryValueView[];
}

export default function FilterTags({
  categoryKey,
  onFilter,
  selectedValues,
}: Props): JSX.Element | null {
  return selectedValues.length ? (
    <SelectedTags>
      {selectedValues.map(({ key }) => (
        <Tag
          key={key}
          large
          minimal
          multiline
          onRemove={() => onFilter(categoryKey, key)}
        >
          {key}
        </Tag>
      ))}
    </SelectedTags>
  ) : null;
}
