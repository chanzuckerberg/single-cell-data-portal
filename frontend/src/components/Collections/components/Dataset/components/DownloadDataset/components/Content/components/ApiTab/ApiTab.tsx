import Name from "../BrowserTab/components/Name";
import CopyCodeBlock from "../CopyCodeBlock";
import { FormControl, FormLabel } from "@mui/material";
import { Link } from "@czi-sds/components";
import {
  ActiveTab,
  CENSUS_API_URL,
  DISCOVER_API_URL,
  NOTEBOOK_URLS,
} from "../../utils";
import { EVENTS } from "src/common/analytics/events";

interface ApiTabsProps {
  name: string;
  tab: ActiveTab.PythonApi | ActiveTab.RApi;
  censusCopyText: string;
  handleAnalytics: (event: EVENTS) => void;
  discoverCopyText: string;
}

export const ApiTab = ({
  name,
  tab,
  censusCopyText,
  discoverCopyText,
  handleAnalytics,
}: ApiTabsProps) => {
  return (
    <>
      <Name name={name} />
      <FormControl>
        <FormLabel>Census Api</FormLabel>
        <CopyCodeBlock
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
            handleAnalytics(EVENTS.DOWNLOAD_DATA_COPY);
          }}
        />
      </FormControl>
      <FormControl>
        <FormLabel>Discover Api</FormLabel>
        <CopyCodeBlock
          copyText={discoverCopyText}
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
              versions. See{" "}
              <Link
                href={NOTEBOOK_URLS[tab]}
                rel="noreferrer noopener"
                sdsSize="xs"
                sdsStyle="default"
                target="_blank"
              >
                this notebook
              </Link>{" "}
              for an example.
            </p>
          }
          handleAnalytics={() => {
            handleAnalytics(EVENTS.DOWNLOAD_DATA_COPY);
          }}
        />
      </FormControl>
    </>
  );
};
