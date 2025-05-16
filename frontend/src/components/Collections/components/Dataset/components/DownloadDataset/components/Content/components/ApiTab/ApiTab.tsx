import Name from "../BrowserTab/components/Name";
import DownloadLink from "../BrowserTab/components/DownloadLink";
import { FormControl, FormLabel } from "@mui/material";
import { Link } from "@czi-sds/components";

interface ApiTabsProps {
  name: string;
  censusCopyText: string;
  discoverCopyText: string;
}

const DISCOVER_API_URL = "https://api.cellxgene.cziscience.com/curation/ui/#/";
const CENSUS_API_URL =
  "https://chanzuckerberg.github.io/cellxgene-census/index.html";

export const ApiTab = ({
  name,
  censusCopyText,
  discoverCopyText,
}: ApiTabsProps) => {
  return (
    <>
      <Name name={name} />

      <div>
        <FormControl>
          <FormLabel>Census Api</FormLabel>
          <DownloadLink
            copyText={censusCopyText}
            captionText={
              <p>
                <Link
                  href={CENSUS_API_URL}
                  rel="noreferrer noopener"
                  sdsSize="xs"
                  sdsStyle="default"
                  target="_blank"
                >
                  Census API
                </Link>{" "}
                provides efficient computational tooling to access, query, and
                analyze all single-cell RNA data from CZ CELLxGENE Discover. You
                can interact with the data through TileDB-SOMA, or get slices in
                AnnData, Seurat, or SingleCellExperiment objects.
              </p>
            }
            handleAnalytics={() => {
              return; // handleAnalytics(EVENTS.DOWNLOAD_CENSUS_API_LINK_CLICKED);
            }}
          />
        </FormControl>
      </div>

      <div>
        <FormControl>
          <FormLabel>Discover Api</FormLabel>
          <DownloadLink
            copyText={`https://api.cellxgene.cziscience.com/curation/v1/${discoverCopyText}`}
            captionText={
              <p>
                <Link
                  href={DISCOVER_API_URL}
                  rel="noreferrer noopener"
                  sdsSize="xs"
                  sdsStyle="default"
                  target="_blank"
                >
                  Discover API
                </Link>{" "}
                provides programmatic access to individual datasets and their
                versions.
              </p>
            }
            handleAnalytics={() => {
              return; // handleAnalytics(EVENTS.DOWNLOAD_CENSUS_API_LINK_CLICKED);
            }}
          />
        </FormControl>
      </div>
    </>
  );
};
