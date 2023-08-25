import { Button, Intent, UL } from "@blueprintjs/core";
import { IconNames } from "@blueprintjs/icons";
import { FC, useCallback, useState } from "react";
import { useQueryClient } from "react-query";
import { ACCESS_TYPE, Collection, Dataset } from "src/common/entities";
import {
  USE_COLLECTION,
  useCollection,
  useCollectionUploadLinks,
  useReuploadDataset,
} from "src/common/queries/collections";
import { isTombstonedCollection } from "src/common/utils/typeGuards";
import { CollectionDatasetsGrid } from "src/components/Collection/components/CollectionDatasetsGrid/style";
import DropboxChooser, { UploadingFile } from "src/components/DropboxChooser";
import { StyledLink } from "src/views/Collection/common/style";
import Toast from "src/views/Collection/components/Toast";
import EmptyModal from "../EmptyModal";
import { UploadedFiles } from "src/views/Collection/components/CollectionActions/components/AddButton";

interface Props {
  collectionId: Collection["id"];
  visibility: Collection["visibility"];
  datasets: Array<Dataset>;
  isRevision: boolean;
}

const DatasetTab: FC<Props> = ({
  collectionId,
  visibility,
  datasets,
  isRevision,
}) => {
  const CLI_README_LINK =
    "https://github.com/chanzuckerberg/single-cell-curation/blob/main/readme.md";

  const { mutateAsync: uploadLink } = useCollectionUploadLinks(collectionId);
  const { mutateAsync: reuploadDataset } = useReuploadDataset(collectionId);
  const [uploadedFiles, setUploadedFiles] = useState({} as UploadedFiles);
  const { data: collection } = useCollection({ id: collectionId });

  const queryClient = useQueryClient();

  const invalidateCollectionQuery = useCallback(() => {
    queryClient.invalidateQueries([USE_COLLECTION, collectionId]);
  }, [collectionId, queryClient]);

  if (isTombstonedCollection(collection)) return null;

  const hasWriteAccess = collection?.access_type === ACCESS_TYPE.WRITE;

  const isDatasetPresent =
    datasets?.length > 0 || Object.keys(uploadedFiles).length > 0;

  const addNewFile = (mutationFunction = uploadLink, originalId?: string) => {
    return (newFile: UploadingFile) => {
      if (!newFile.link) return;

      const payload = JSON.stringify({ id: originalId, url: newFile.link });
      mutationFunction(
        { collectionId: collectionId, payload },
        {
          onSuccess: (datasetID: Dataset["id"]) => {
            newFile.id = datasetID;
            Toast.show({
              icon: IconNames.TICK,
              intent: Intent.PRIMARY,
              message:
                "Your file is being uploaded which will continue in the background, even if you close this window.",
            });
            setUploadedFiles({ ...uploadedFiles, [newFile.id]: newFile });
          },
        }
      );
    };
  };

  return (
    <>
      {isDatasetPresent ? (
        <CollectionDatasetsGrid
          accessType={collection?.access_type}
          collectionId={collectionId}
          datasets={datasets}
          invalidateCollectionQuery={invalidateCollectionQuery}
          isRevision={isRevision}
          onUploadFile={addNewFile}
          reuploadDataset={reuploadDataset}
          uploadedFiles={uploadedFiles}
          visibility={visibility}
        />
      ) : (
        <EmptyModal
          title="No datasets uploaded"
          content={
            hasWriteAccess ? (
              <div>
                Before you begin uploading dataset files:
                <UL>
                  <li>
                    You must validate your dataset locally. We provide a local
                    CLI script to do this.{" "}
                    <StyledLink href={CLI_README_LINK}>Learn More</StyledLink>
                  </li>
                  <li>
                    We only support adding datasets in the h5ad format at this
                    time.
                  </li>
                </UL>
              </div>
            ) : undefined
          }
          button={
            hasWriteAccess ? (
              <DropboxChooser onUploadFile={addNewFile()}>
                <Button
                  intent={Intent.PRIMARY}
                  outlined
                  text={"Add Dataset from Dropbox"}
                />
              </DropboxChooser>
            ) : undefined
          }
        />
      )}
    </>
  );
};
export default DatasetTab;
