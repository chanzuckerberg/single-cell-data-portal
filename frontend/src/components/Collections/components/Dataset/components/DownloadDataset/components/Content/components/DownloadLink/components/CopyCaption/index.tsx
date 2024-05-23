import { Caption, StyledLink } from "./style";
import { DATASET_ASSET_FORMAT } from "src/common/entities";

const DISCOVER_API_URL = "https://api.cellxgene.cziscience.com/curation/ui/#/";
const SCHEMA_URL =
  "https://github.com/chanzuckerberg/single-cell-curation/blob/main/schema/5.0.0/schema.md";
const SEURAT_SCHEMA_URL =
  "https://github.com/chanzuckerberg/single-cell-curation/blob/main/schema/5.0.0/seurat_encoding.md";

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
        <StyledLink
          href={DISCOVER_API_URL}
          rel="noreferrer noopener"
          target="_blank"
        >
          Discover API
        </StyledLink>
        . The{" "}
        <StyledLink href={SCHEMA_URL} rel="noreferrer noopener" target="_blank">
          dataset schema
        </StyledLink>{" "}
        describes the required metadata embedded in all datasets submitted to CZ
        CELLxGENE Discover.{" "}
        {selectedFormat === DATASET_ASSET_FORMAT.RDS && (
          <>
            All datasets are automatically converted to a{" "}
            <StyledLink
              href={SEURAT_SCHEMA_URL}
              rel="noreferrer noopener"
              target="_blank"
            >
              Seurat v5 object
            </StyledLink>
            .
          </>
        )}
      </p>
    </Caption>
  );
}
