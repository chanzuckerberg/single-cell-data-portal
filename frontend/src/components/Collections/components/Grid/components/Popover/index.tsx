import { PluralizedMetadataLabel } from "src/common/constants/metadata";
import NTag from "src/components/common/Grid/components/NTag";
import { FieldValues } from "src/components/common/Grid/components/NTag/style";
import { LeftAlignedDetailsCell } from "../Row/common/style";

interface Props {
  label: PluralizedMetadataLabel;
  values: string[];
}

export default function Popover({ label, values }: Props): JSX.Element {
  return (
    <LeftAlignedDetailsCell>
      {values.length <= 2 ? (
        <FieldValues>
          {values[0]}
          <br />
          {values[1]}
        </FieldValues>
      ) : (
        <NTag label={label} values={values} />
      )}
    </LeftAlignedDetailsCell>
  );
}
