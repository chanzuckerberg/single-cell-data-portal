import { DataPortalProps } from "../../common/types/entities";

const VERSION_ONES = ["1.0.0", "1.1.0"];
const VERSION_TWOS = ["2.0.0"];

export function checkValidVersion(corporaProps: DataPortalProps): boolean {
  if (!corporaProps) return false;

  const { version, schema_version: schemaVersion } = corporaProps;

  const isValidVersionOne = VERSION_ONES.includes(
    version?.corpora_schema_version,
  );
  const isValidVersionTwo = VERSION_TWOS.includes(schemaVersion);

  return isValidVersionOne || isValidVersionTwo;
}
