import { Tag } from "@blueprintjs/core";
import React from "react";
import { SelectedTags } from "./style";

type OnRemoveFn = () => void;

export interface CategoryTag {
  label: string;
  onRemove: OnRemoveFn;
}

interface Props {
  tags?: CategoryTag[];
}

export default function FilterTags({ tags }: Props): JSX.Element | null {
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
