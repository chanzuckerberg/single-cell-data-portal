import React from "react";
import { PLURALIZED_METADATA_LABEL } from "src/common/constants/metadata";
import NTagCell, {
  TagValue,
} from "src/components/common/Grid/components/NTagCell";

interface Props {
  label: PLURALIZED_METADATA_LABEL;
  values: string[];
}

const DISEASE_CATEGORY_VALUE = {
  NORMAL: "normal",
};

export default function DiseaseCell({ label, values }: Props): JSX.Element {
  const diseases = mapDiseasesToTagValues(values);
  return <NTagCell label={label} values={diseases} />;
}

/**
 * Diseases are mapped to an array of tag values to facilitate the permanent display of the disease "normal"
 * above the remaining +ntag values.
 * @param values - String array of diseases.
 * @returns TagValues - Diseases where the value "normal" is pinned.
 */
function mapDiseasesToTagValues(values: string[]): TagValue[] {
  return values.map((value) => {
    let pinned = false;
    if (value === DISEASE_CATEGORY_VALUE.NORMAL) {
      pinned = true;
    }
    return { label: value, pinned };
  });
}
