import { Icon, TagFilter } from "@czi-sds/components";
import { SelectedTags } from "./style";

type OnRemoveFn = () => void;

export interface CategoryTag {
  label: string | [string, string]; // Single value ("label") or range tuple ("label[0] &ndash; label[1]")
  onRemove: OnRemoveFn;
}

interface Props {
  tags?: CategoryTag[];
}

export default function FilterTags({ tags }: Props): JSX.Element | null {
  return tags && tags.length ? (
    <SelectedTags>
      {tags.map(({ label, onRemove }, i) => (
        <TagFilter
          key={isLabelRange(label) ? label.join("") : `${label}${i}`}
          deleteIcon={<Icon sdsIcon="xMark" sdsSize="xs" sdsType="button" />}
          label={isLabelRange(label) ? `${label[0]} - ${label[1]}` : label}
          onClick={onRemove}
          onDelete={onRemove}
        />
      ))}
    </SelectedTags>
  ) : null;
}

/**
 * Returns true if the given label is multi value.
 * @param label - Label to check if it is a single string or a string array.
 * @returns True if label is multi value.
 */
function isLabelRange(label: string | string[]): label is [string, string] {
  return Array.isArray(label);
}
