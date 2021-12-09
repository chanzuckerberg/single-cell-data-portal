import { Tag } from "@blueprintjs/core";
import {
  CategoryValueView,
  CATEGORY_KEY,
  OnFilterFn,
} from "src/components/common/Filter/common/entities";
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
      {selectedValues.map(({ key, label }) => (
        <Tag
          key={key}
          large
          minimal
          multiline
          onRemove={() => onFilter(categoryKey, key)}
        >
          {label}
        </Tag>
      ))}
    </SelectedTags>
  ) : null;
}
