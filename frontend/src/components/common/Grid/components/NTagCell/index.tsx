import { Box } from "@mui/material";
import React from "react";
import { PLURALIZED_METADATA_LABEL } from "src/common/constants/metadata";
import NTag from "src/components/common/Grid/components/NTag";

type PartitionedCellValues = [string | undefined, string[], string[]];

export interface TagValue {
  label: string;
  pinned: boolean;
}

interface Props {
  label: PLURALIZED_METADATA_LABEL;
  values: string[] | TagValue[];
}

const MAX_DISPLAYABLE_VALUES = 2;

export default function NTagCell({ label, values }: Props): JSX.Element {
  const [pinned, unpinned, allValues] = partitionCellValues(values);
  return allValues.length <= MAX_DISPLAYABLE_VALUES ? (
    <>
      {/* display all values */}
      {allValues.map((v: string) => (
        <div key={v}>{v}</div>
      ))}
    </>
  ) : (
    <>
      {/* optionally pinned value */}
      {!!pinned && <Box mb={2}>{pinned}</Box>}
      {/* +ntag values */}
      <NTag label={label} values={unpinned} />
    </>
  );
}

/**
 * Determine if the given cell value is a tag value and not a string.
 * @param value - Cell value, either a string or a tag value.
 * @returns True if the given cell value is a tag value.
 */
function isTagValue(value: string | TagValue): value is TagValue {
  return (value as TagValue).label !== undefined;
}

/**
 * Split cell values into a single pinned value, and arrays of non-pinned values and all cell values.
 * The partition facilitates the rendering of a cell with the following conditions:
 * - all cell values for permanent display of two or less values, or
 * - an (optional) pinned value for permanent display above the +ntag values, and
 * - non-pinned values for +ntag values.
 * @param values - Cell values, either a string array or an array of tag values.
 * @returns Tuple containing a single pinned value, a string array of non-pinned values,
 * and a string array of all cell values.
 */
function partitionCellValues(
  values: string[] | TagValue[]
): PartitionedCellValues {
  const partitionedValues: PartitionedCellValues = [undefined, [], []];
  // @ts-expect-error --- suppressing TS limitation (see https://github.com/microsoft/TypeScript/issues/36390 and
  // https://github.com/microsoft/TypeScript/issues/44593#issuecomment-861630080).
  // reduce() is not callable on array union types.
  return values.reduce(
    (accum: PartitionedCellValues, value: string | TagValue) => {
      const [, unpinned, allValues] = accum;
      if (!isTagValue(value)) {
        unpinned.push(value);
        allValues.push(value);
        return accum;
      }
      if (value.pinned) {
        accum[0] = value.label;
      } else {
        unpinned.push(value.label);
      }
      allValues.push(value.label);
      return accum;
    },
    partitionedValues
  );
}
