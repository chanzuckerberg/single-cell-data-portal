import { Caption } from "./style";
import { DATASET_ASSET_FORMAT } from "src/common/entities";
import { Link } from "@czi-sds/components";

const DISCOVER_API_URL = "https://api.cellxgene.cziscience.com/curation/ui/#/";
const SCHEMA_URL =
  "https://chanzuckerberg.github.io/single-cell-curation/latest-schema.html";
const SEURAT_SCHEMA_URL =
  "https://github.com/chanzuckerberg/single-cell-curation/blob/main/schema/5.1.0/seurat_encoding.md";

interface Props {
  selectedFormat: DATASET_ASSET_FORMAT | "";
}

export default function CopyCaption({ selectedFormat }: Props): JSX.Element {
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
        {selectedFormat === DATASET_ASSET_FORMAT.RDS && (
          <>
            All datasets are automatically converted to a{" "}
            <Link
              href={SEURAT_SCHEMA_URL}
              rel="noreferrer noopener"
              sdsSize="xs"
              sdsStyle="default"
              target="_blank"
            >
              Seurat v5 object
            </Link>
            .
          </>
        )}
      </p>
    </Caption>
  );
}
