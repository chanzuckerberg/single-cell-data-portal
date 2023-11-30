import { Caption } from "./style";
import { Link } from "@czi-sds/components";
import { DATASET_ASSET_FORMAT } from "src/common/entities";

const DISCOVER_API_URL = "https://api.cellxgene.cziscience.com/curation/ui/#/";
const SCHEMA_URL =
  "https://github.com/chanzuckerberg/single-cell-curation/blob/main/schema/3.1.0/schema.md";
const SEURAT_SCHEMA_URL =
  "https://github.com/chanzuckerberg/single-cell-curation/blob/main/schema/3.1.0/seurat_encoding.md";

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
        <Link href={DISCOVER_API_URL} rel="noreferrer noopener" target="_blank">
          Discover API
        </Link>
        . The{" "}
        <Link href={SCHEMA_URL} rel="noreferrer noopener" target="_blank">
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
              target="_blank"
            >
              Seurat V4 object
            </Link>
            .
          </>
        )}
      </p>
    </Caption>
  );
}
