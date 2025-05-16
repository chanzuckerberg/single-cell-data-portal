import { DialogLoader as Loader } from "src/components/Datasets/components/DownloadDataset/style";
import { DATASET_ASSET_FORMAT } from "src/common/entities";
import DownloadLink from "./components/DownloadLink";
import DataFormat from "./components/DataFormat";
import Details from "./components/Details";
import Name from "./components/Name";
import { DownloadLinkType } from "../../index";
import { EVENTS } from "src/common/analytics/events";
import { Link } from "@czi-sds/components";

interface BrowserTabProps {
  name: string;
  isError: boolean;
  isLoading: boolean;
  availableFormats: Set<DATASET_ASSET_FORMAT>;
  handleChange: (format: DATASET_ASSET_FORMAT) => void;
  downloadLinks: DownloadLinkType[];
  formatsToDownload: DATASET_ASSET_FORMAT[];
  handleAnalytics: (event: EVENTS, dataFormat?: DATASET_ASSET_FORMAT) => void;
  selectedFormats: DATASET_ASSET_FORMAT[];
  isDownloadLinkLoading: boolean;
}

const DISCOVER_API_URL = "https://api.cellxgene.cziscience.com/curation/ui/#/";
const SCHEMA_URL =
  "https://chanzuckerberg.github.io/single-cell-curation/latest-schema.html";

export const BrowserTab = ({
  name,
  isError,
  isLoading,
  availableFormats,
  handleChange,
  downloadLinks,
  formatsToDownload,
  handleAnalytics,
  selectedFormats,
  isDownloadLinkLoading,
}: BrowserTabProps) => {
  const copyText = downloadLinks.reduce((acc, download) => {
    if (!formatsToDownload.includes(download.filetype)) {
      return acc;
    }
    return acc + download.downloadURL + "\n";
  }, "");
  return (
    <>
      {isError && <div>Dataset download is currently not available.</div>}
      {isLoading && <Loader sdsStyle="minimal" />}
      {!isError && !isLoading && (
        <>
          <Name name={name} />
          <DataFormat
            availableFormats={availableFormats}
            handleChange={handleChange}
            selectedFormats={selectedFormats}
            downloadLinks={downloadLinks}
            isDisabled={downloadLinks.length === 0}
          />
          <Details
            downloadPreview={
              <DownloadLink
                copyText={copyText}
                plural={downloadLinks.length > 1}
                handleAnalytics={() =>
                  handleAnalytics(EVENTS.DOWNLOAD_DATA_COPY)
                }
                captionText={
                  <>
                    <p>
                      This download link permanently references{" "}
                      <b>this version</b> of the dataset. If this dataset is
                      updated, a new download link will be created that
                      permanently references the next version of this dataset.
                    </p>
                    <p>
                      Individual datasets and their versions may also be
                      downloaded programmatically using the{" "}
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
                      describes the required metadata embedded in all datasets
                      submitted to CZ CELLxGENE Discover.{" "}
                    </p>
                  </>
                }
              />
            }
            isLoading={isDownloadLinkLoading}
            hasDownloadLinks={downloadLinks.length > 0}
            selected={selectedFormats.length > 0}
          />
        </>
      )}
    </>
  );
};
