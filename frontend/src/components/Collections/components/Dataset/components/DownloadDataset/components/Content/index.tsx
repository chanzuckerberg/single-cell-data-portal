import React, { FC, useEffect, useMemo, useState } from "react";
import { Dataset, DATASET_ASSET_FORMAT } from "src/common/entities";
import DownloadLink from "./components/DownloadLink";
import DataFormat from "./components/DataFormat";
import Details from "./components/Details";
import Name from "./components/Name";
import { track } from "src/common/analytics";
import { EVENTS } from "src/common/analytics/events";
import { DialogActions, DialogContent, DialogTitle } from "@czi-sds/components";
import { DialogLoader as Loader } from "src/components/Datasets/components/DownloadDataset/style";
import { Button } from "src/components/common/Button";
import {
  POSSIBLE_DOWNLOAD_FORMATS,
  getNotAvailableText,
} from "./components/DataFormat/constants";
import { downloadMultipleFiles, getDownloadLink } from "./utils";
import { Tooltip } from "@czi-sds/components";
interface Props {
  isError?: boolean;
  isLoading?: boolean;
  onClose: () => void;
  name: string;
  dataAssets: Dataset["dataset_assets"];
}

export interface DownloadLinkType {
  filetype: DATASET_ASSET_FORMAT;
  fileSize: number | undefined;
  downloadURL: string | undefined;
}

const Content: FC<Props> = ({
  isError = false,
  isLoading = false,
  onClose,
  name,
  dataAssets,
}) => {
  const [selectedFormats, setSelectedFormats] = useState<
    DATASET_ASSET_FORMAT[]
  >([]);
  const [downloadLinks, setDownloadLinks] = useState<DownloadLinkType[]>([]);
  const [isDownloadLinkLoading, setIsDownloadLinkLoading] =
    useState<boolean>(false);

  const noDownloadLinksAvailable = useMemo(() => {
    return (
      selectedFormats.length > 0 &&
      selectedFormats.every((format) => {
        const link = downloadLinks.find((l) => l.filetype === format);

        return (
          !link || // no matching link
          link.downloadURL === undefined ||
          link.downloadURL === getNotAvailableText(format)
        );
      })
    );
  }, [selectedFormats, downloadLinks]);

  const isDownloadDisabled =
    !selectedFormats.length ||
    isDownloadLinkLoading ||
    isError ||
    isLoading ||
    noDownloadLinksAvailable;

  const availableFormats = useMemo(
    () => new Set(dataAssets.map((dataAsset) => dataAsset.filetype)),
    [dataAssets]
  );

  // Determine formats to download, ensuring ATAC_INDEX includes ATAC_FRAGMENT
  const formatsToDownload = useMemo(() => {
    return selectedFormats.flatMap((format) =>
      format === DATASET_ASSET_FORMAT.ATAC_INDEX
        ? [DATASET_ASSET_FORMAT.ATAC_INDEX, DATASET_ASSET_FORMAT.ATAC_FRAGMENT]
        : [format]
    );
  }, [selectedFormats]);

  // Set selected formats based on available formats and ATAC index
  useEffect(() => {
    const isATACIncomplete = (format: DATASET_ASSET_FORMAT) => {
      if (format !== DATASET_ASSET_FORMAT.ATAC_INDEX) return false;
      const hasIndex = availableFormats.has(DATASET_ASSET_FORMAT.ATAC_INDEX);
      const hasFragment = availableFormats.has(
        DATASET_ASSET_FORMAT.ATAC_FRAGMENT
      );
      return !(hasIndex && hasFragment);
    };
    const allFormats = POSSIBLE_DOWNLOAD_FORMATS.filter(
      ({ format }) => availableFormats.has(format) && !isATACIncomplete(format)
    ).map(({ format }) => format);

    const selectedFormats: DATASET_ASSET_FORMAT[] = [];

    // Default to just H5AD selected if available
    allFormats.includes(DATASET_ASSET_FORMAT.H5AD)
      ? selectedFormats.push(DATASET_ASSET_FORMAT.H5AD)
      : selectedFormats.push(allFormats[0]);

    setSelectedFormats(selectedFormats);
  }, [availableFormats, dataAssets]);

  useEffect(() => {
    const shouldSkipFetching = () => {
      if (selectedFormats.length === 0 || formatsToDownload.length === 0) {
        return true;
      }
      return formatsToDownload.every((format) =>
        downloadLinks.some((link) => link.filetype === format)
      );
    };
    const getAssetsToRetrieve = () => {
      return dataAssets.filter(
        (asset) =>
          !downloadLinks.some((link) => link.filetype === asset.filetype) &&
          formatsToDownload.includes(asset.filetype)
      );
    };
    const replaceMissingURLs = (links: (DownloadLinkType | null)[]) => {
      return links
        .filter((link) => link !== null)
        .map((link) => ({
          ...link,
          downloadURL: link.downloadURL ?? getNotAvailableText(link.filetype),
        }));
    };

    const fetchDownloadLinks = async () => {
      setIsDownloadLinkLoading(true);

      if (shouldSkipFetching()) {
        setIsDownloadLinkLoading(false);
        return;
      }

      const assetsToRetrieve = getAssetsToRetrieve();
      if (assetsToRetrieve.length === 0) {
        setIsDownloadLinkLoading(false);
        console.error("No Download Formats available");
        return;
      }

      try {
        const links = await Promise.all(assetsToRetrieve.map(getDownloadLink));
        const cleanedLinks = replaceMissingURLs(links);
        setDownloadLinks((prev) => [...new Set([...prev, ...cleanedLinks])]);
      } catch (error) {
        console.error("Error fetching download links", error);
      } finally {
        setIsDownloadLinkLoading(false);
      }
    };

    fetchDownloadLinks();
  }, [selectedFormats, dataAssets, downloadLinks, formatsToDownload]);
  /**
   * Tracks dataset download analytics as specified by the custom analytics event.
   * @param event - Custom analytics event.
   * @param dataFormat - Data format (optionally specified with event "DOWNLOAD_DATA_FORMAT_CLICKED").
   */
  const handleAnalytics = (
    event: EVENTS,
    dataFormat?: DATASET_ASSET_FORMAT
  ) => {
    track(event, {
      dataset_name: name,
      data_format: dataFormat || selectedFormats,
    });
  };

  const handleChange = (format: DATASET_ASSET_FORMAT) => {
    if (selectedFormats.includes(format)) {
      setSelectedFormats(selectedFormats.filter((f) => f !== format));
    } else {
      setSelectedFormats([...selectedFormats, format]);
    }
    handleAnalytics(EVENTS.DOWNLOAD_DATA_FORMAT_CLICKED, format);
  };

  return (
    <>
      <DialogTitle title="Download Data" />
      <DialogContent>
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
                  formatsToDownload={formatsToDownload}
                  downloadLinks={downloadLinks}
                  handleAnalytics={() =>
                    handleAnalytics(EVENTS.DOWNLOAD_DATA_COPY)
                  }
                />
              }
              isLoading={isDownloadLinkLoading}
              hasDownloadLinks={downloadLinks.length > 0}
              selected={selectedFormats.length > 0}
            />
          </>
        )}
      </DialogContent>
      <DialogActions>
        <Button
          isAllCaps={false}
          onClick={onClose}
          sdsStyle="minimal"
          sdsType="secondary"
        >
          Cancel
        </Button>
        <Tooltip
          placement="top"
          disableHoverListener={!isDownloadDisabled}
          disableFocusListener={!isDownloadDisabled}
          disableTouchListener={!isDownloadDisabled}
          title={"Select at least one valid format"}
        >
          <span>
            <Button
              data-testid="download-asset-download-button"
              disabled={isDownloadDisabled}
              onClick={() => {
                downloadMultipleFiles(formatsToDownload, downloadLinks);
                handleAnalytics(EVENTS.DOWNLOAD_DATA_COMPLETE);
              }}
              sdsStyle="square"
              sdsType="primary"
            >
              Download
            </Button>
          </span>
        </Tooltip>
      </DialogActions>
    </>
  );
};

export default Content;
