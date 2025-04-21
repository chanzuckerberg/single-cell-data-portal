import { Caption } from "./style";
import { Link } from "@czi-sds/components";

const DISCOVER_API_URL = "https://api.cellxgene.cziscience.com/curation/ui/#/";
const SCHEMA_URL =
  "https://chanzuckerberg.github.io/single-cell-curation/latest-schema.html";

export default function CopyCaption(): JSX.Element {
  return (
    <Caption>
      <p>
        This download link permanently references <b>this version</b> of the
        dataset. If this dataset is updated, a new download link will be created
        that permanently references the next version of this dataset.
      </p>
      <p>
        Individual datasets and their versions may also be downloaded
        programmatically using the{" "}
        <Link
          href={DISCOVER_API_URL}
          rel="noreferrer noopener"
          sdsSize="xs"
          sdsStyle="default"
          target="_blank"
        >
          Discover API
        </Link>
        . The{" "}
        <Link
          href={SCHEMA_URL}
          rel="noreferrer noopener"
          sdsSize="xs"
          sdsStyle="default"
          target="_blank"
        >
          dataset schema
        </Link>{" "}
        describes the required metadata embedded in all datasets submitted to CZ
        CELLxGENE Discover.{" "}
      </p>
    </Caption>
  );
}
