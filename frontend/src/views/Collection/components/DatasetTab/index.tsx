import { Button, Intent, UL } from "@blueprintjs/core";
import { IconNames } from "@blueprintjs/icons";
import memoize from "lodash/memoize";
import React, { FC, useState } from "react";
import { MutateFunction, useQueryCache } from "react-query";
import { Collection, Dataset } from "src/common/entities";
import {
  useCollection,
  useCollectionUploadLinks,
  useReuploadDataset,
  USE_COLLECTION,
} from "src/common/queries/collections";
import DatasetsGrid from "src/components/Collections/components/Grid/components/DatasetsGrid";
import DropboxChooser, { UploadingFile } from "src/components/DropboxChooser";
import { StyledLink } from "src/views/Collection/common/style";
import { UploadedFiles } from "src/views/Collection/components/ActionButtons";
import DatasetUploadToast from "src/views/Collection/components/DatasetUploadToast";
import EmptyModal from "../EmptyModal";

interface Props {
  collectionID: Collection["id"];
  visibility: Collection["visibility"];
  datasets: Array<Dataset>;
  isRevision: boolean;
}

const DatasetTab: FC<Props> = ({
  collectionID: collectionId,
  visibility,
  datasets,
  isRevision,
}) => {
  const CLI_README_LINK =
    "https://github.com/chanzuckerberg/cellxgene/blob/main/dev_docs/schema_guide.md";

  const [uploadLink] = useCollectionUploadLinks(collectionId, visibility);
  const [reuploadDataset] = useReuploadDataset(collectionId);
  const [uploadedFiles, setUploadedFiles] = useState({} as UploadedFiles);
  const { data: collection } = useCollection({ id: collectionId, visibility });

  const queryCache = useQueryCache();

  const isDatasetPresent =
    datasets?.length > 0 || Object.keys(uploadedFiles).length > 0;

  const invalidateCollectionQuery = memoize(
    () => {
      queryCache.invalidateQueries([USE_COLLECTION, collectionId, visibility]);
    },
    () => collectionId + visibility
  );

  const addNewFile = (
    mutationFunction = uploadLink as MutateFunction<string, unknown, unknown>,
    originalId?: string
  ) => {
    return (newFile: UploadingFile) => {
      if (!newFile.link) return;

      const payload = JSON.stringify({ url: newFile.link, uuid: originalId });
      mutationFunction(
        { collectionId: collectionId, payload },
        {
          onSuccess: (datasetID: Dataset["id"]) => {
            newFile.id = datasetID;
            DatasetUploadToast.show({
              icon: IconNames.TICK,
              intent: Intent.PRIMARY,
              message:
                "Your file is being uploaded which will continue in the background, even if you close this window.",
            });
            setUploadedFiles({ ...uploadedFiles, [newFile.id]: newFile });
            queryCache.invalidateQueries(USE_COLLECTION);
          },
        }
      );
    };
  };

  return (
    <>
      {isDatasetPresent ? (
        <DatasetsGrid
          visibility={visibility}
          accessType={collection?.access_type}
          datasets={datasets}
          uploadedFiles={uploadedFiles}
          invalidateCollectionQuery={invalidateCollectionQuery}
          isRevision={isRevision}
          onUploadFile={addNewFile}
          reuploadDataset={reuploadDataset}
        />
      ) : (
        <EmptyModal
          title="No datasets uploaded"
          content={
            <div>
              Before you begin uploading dataset files:
              <UL>
                <li>
                  You must validate your dataset locally. We provide a local CLI
                  script to do this.{" "}
                  <StyledLink href={CLI_README_LINK}>Learn More</StyledLink>
                </li>
                <li>
                  We only support adding datasets in the h5ad format at this
                  time.
                </li>
              </UL>
            </div>
          }
          button={
            <DropboxChooser onUploadFile={addNewFile()}>
              <Button
                intent={Intent.PRIMARY}
                outlined
                text={"Add Dataset from Dropbox"}
              />
            </DropboxChooser>
          }
        />
      )}
    </>
  );
};
export default DatasetTab;
