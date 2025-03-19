// // TODO: (smccanny) delete this before merging
/* eslint-disable sonarjs/no-duplicate-string */
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
import { possibleDownloadFormats } from "./components/DataFormat/constants";
import { downloadMultipleFiles, getDownloadLink } from "./utils";

// // TODO: (smccanny) delete this hard coded test array before merging
// const dataAssets = [
//   {
//     created_at: 0,
//     dataset_id: "30a7b283-82db-408f-92fa-238998378eb5",
//     filename: "raw.h5ad",
//     filetype: "RAW_H5AD",
//     id: "14aa6934-b111-4738-90f4-feadc3584e2e",
//     s3_uri:
//       "s3://corpora-data-staging/30a7b283-82db-408f-92fa-238998378eb5/raw.h5ad",
//     updated_at: 0,
//     user_submitted: true,
//   },
//   {
//     created_at: 0,
//     dataset_id: "30a7b283-82db-408f-92fa-238998378eb5",
//     filename: "local.h5ad",
//     filetype: "H5AD",
//     id: "8a7f782f-0ab3-4118-bea3-0706eca08f12",
//     s3_uri:
//       "s3://corpora-data-staging/30a7b283-82db-408f-92fa-238998378eb5/local.h5ad",
//     updated_at: 0,
//     user_submitted: true,
//   },
//   {
//     created_at: 0,
//     dataset_id: "30a7b283-82db-408f-92fa-238998378eb5",
//     filename: "",
//     filetype: "CXG",
//     id: "1ec622f0-f0ef-4a50-983b-bc7774b66dcc",
//     s3_uri:
//       "s3://hosted-cellxgene-staging/30a7b283-82db-408f-92fa-238998378eb5.cxg/",
//     updated_at: 0,
//     user_submitted: true,
//   },
//   {
//     created_at: 0,
//     dataset_id: "30a7b283-82db-408f-92fa-238998378eb5",
//     filename: "local.rds",
//     filetype: "RDS",
//     id: "4dd72160-f924-447f-93f4-ac5253ce4606",
//     s3_uri:
//       "s3://corpora-data-staging/30a7b283-82db-408f-92fa-238998378eb5/local.rds",
//     updated_at: 0,
//     user_submitted: true,
//   },
//   {
//     created_at: 0,
//     dataset_id: "30a7b283-82db-408f-92fa-238998378eb5",
//     filename: "tsv.gz",
//     filetype: "ATAC_INDEX",
//     id: "4dd72160-f924-447f-93f4-ac5253ce4606",
//     s3_uri:
//       "s3://corpora-data-staging/30a7b283-82db-408f-92fa-238998378eb5/local.rds",
//     updated_at: 0,
//     user_submitted: true,
//   },
//   {
//     dataset_id: "30a7b283-82db-408f-92fa-238998378eb5",
//     filename: "fragment",
//     filetype: "ATAC_FRAGMENT",
//     id: "4dd72160-f924-447f-93f4-ac5253ce4606",
//     s3_uri:
//       "s3://corpora-data-staging/30a7b283-82db-408f-92fa-238998378eb5/local.rds",
//   },
// ] as Dataset["dataset_assets"];

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

  const isDownloadDisabled =
    !selectedFormats.length || isDownloadLinkLoading || isError || isLoading;

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

  useEffect(() => {
    const isATACIncomplete = (format: DATASET_ASSET_FORMAT) => {
      if (format !== DATASET_ASSET_FORMAT.ATAC_INDEX) return false;
      const hasIndex = availableFormats.has(DATASET_ASSET_FORMAT.ATAC_INDEX);
      const hasFragment = availableFormats.has(
        DATASET_ASSET_FORMAT.ATAC_FRAGMENT
      );
      return !(hasIndex && hasFragment);
    };

    const selectedFormats = possibleDownloadFormats
      .filter(
        ({ format }) =>
          availableFormats.has(format) && !isATACIncomplete(format)
      )
      .map(({ format }) => format);

    setSelectedFormats(selectedFormats);
  }, [availableFormats, dataAssets]);

  useEffect(() => {
    const fetchDownloadLinks = async () => {
      setIsDownloadLinkLoading(true);
      if (selectedFormats.length === 0) {
        setIsDownloadLinkLoading(false);
        return;
      }
      if (formatsToDownload.length === 0) return;

      // Exit early if all formats already have download links
      if (
        formatsToDownload.every((format) =>
          downloadLinks.some((link) => link.filetype === format)
        )
      ) {
        setIsDownloadLinkLoading(false);
        return;
      }

      // Filter assets for selected formats
      const assets = dataAssets.filter((asset) =>
        formatsToDownload.includes(asset.filetype)
      );
      if (assets.length === 0) {
        setIsDownloadLinkLoading(false);
        throw new Error("No Download Formats available");
      }

      try {
        // Fetch download links and update state
        const newLinks = await Promise.all(assets.map(getDownloadLink));
        const newLinksFiltered = newLinks.filter(
          (link) => link !== null
        ) as DownloadLinkType[];
        setDownloadLinks((prev) => [
          ...new Set([...prev, ...newLinksFiltered]),
        ]);
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
      </DialogActions>
    </>
  );
};

export default Content;
