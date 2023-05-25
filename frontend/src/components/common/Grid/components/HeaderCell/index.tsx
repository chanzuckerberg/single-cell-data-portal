import React, { ReactNode } from "react";
import {
  Header,
  HeaderCell as Cell,
} from "src/components/common/Grid/components/HeaderCell/style";
import { ALIGNMENT } from "src/components/common/Grid/common/entities";

interface Props {
  alignment?: ALIGNMENT;
  label: ReactNode;
  tag?: ReactNode;
}

export default function HeaderCell({
  alignment = ALIGNMENT.LEFT,
  label,
  tag,
}: Props): JSX.Element {
  return (
    <Cell alignment={alignment}>
      <Header>{label}</Header>
      {tag}
    </Cell>
  );
}
