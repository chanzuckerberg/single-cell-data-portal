import { Tag } from "@blueprintjs/core";
import React from "react";
import { SelectedTags } from "./style";

type OnRemoveFn = () => void;

export interface CategoryTag {
  label: string;
  onRemove: OnRemoveFn;
}

interface Props {
  rangeCategoryTag?: CategoryTag;
  selectCategoryTags?: CategoryTag[];
}

export default function FilterTags({
  rangeCategoryTag,
  selectCategoryTags,
}: Props): JSX.Element | null {
  const tags = rangeCategoryTag ? [rangeCategoryTag] : selectCategoryTags;
  return tags && tags.length ? (
    <SelectedTags>
      {tags.map(({ label, onRemove }) => (
        <Tag key={label} large minimal multiline onRemove={onRemove}>
          {label}
        </Tag>
      ))}
    </SelectedTags>
  ) : null;
}
