import React from "react";
import { PLURALIZED_METADATA_LABEL } from "src/common/constants/metadata";
import NTag from "src/components/common/Grid/components/NTag";

interface Props {
  label: PLURALIZED_METADATA_LABEL;
  values: string[];
}

const MAX_DISPLAYABLE_VALUES = 2;

export default function NTagCell({ label, values }: Props): JSX.Element {
  return values.length <= MAX_DISPLAYABLE_VALUES ? (
    <>
      {values.map((v: string) => (
        <div key={v}>{v}</div>
      ))}
    </>
  ) : (
    <NTag label={label} values={values} />
  );
}
