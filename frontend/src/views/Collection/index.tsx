import { Button, H3, Intent } from "@blueprintjs/core";
import { IconNames } from "@blueprintjs/icons";
import { RouteComponentProps } from "@reach/router";
import { memoize } from "lodash-es";
import React, { FC, useState } from "react";
import { useQueryCache } from "react-query";
import { Dataset, VISIBILITY_TYPE } from "src/common/entities";
import {
  useCollection,
  useCollectionUploadLinks,
  USE_COLLECTION,
} from "src/common/queries/collections";
import DeleteCollection from "src/components/Collections/components/DeleteCollection";
import DatasetsGrid from "src/components/Collections/components/Grid/components/DatasetsGrid";
import PublishCollection from "src/components/Collections/components/PublishCollection";
import { UploadingFile } from "src/components/DropboxChooser";
import { ViewGrid } from "../globalStyle";
import ActionButtons, { UploadedFiles } from "./components/ActionButtons";
import DatasetUploadToast from "./components/DatasetUploadToast";
import EmptyDatasets from "./components/EmptyDatasets";
import {
  CollectionButtons,
  CollectionInfo,
  DatasetContainer,
  Description,
  LinkContainer,
} from "./style";
import { getIsPublishable, renderLinks } from "./utils";

interface RouteProps {
  id?: string;
}

export type Props = RouteComponentProps<RouteProps>;

const Collection: FC<Props> = ({ id = "" }) => {
  const isPrivate = window.location.pathname.includes("/private");

  const visibility = isPrivate
    ? VISIBILITY_TYPE.PRIVATE
    : VISIBILITY_TYPE.PUBLIC;

  const [selected, setSelected] = useState<Dataset["id"]>("");

  const [uploadedFiles, setUploadedFiles] = useState({} as UploadedFiles);

  const queryCache = useQueryCache();

  const { data: collection, isError } = useCollection(id, visibility);

  const [uploadLink] = useCollectionUploadLinks(id, visibility);

  const addNewFile = (newFile: UploadingFile) => {
    if (!newFile.link) return;

    const payload = JSON.stringify({ url: newFile.link });
    uploadLink(
      { collectionId: id, payload },
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

  if (!collection || isError) {
    return null;
  }

  const datasets = collection.datasets;

  const isDatasetPresent =
    datasets?.length > 0 || Object.keys(uploadedFiles).length > 0;

  const isPublishable = getIsPublishable(datasets);

  const invalidateCollectionQuery = memoize(
    () => {
      queryCache.invalidateQueries([USE_COLLECTION, id, visibility]);
    },
    () => id + visibility
  );

  return (
    <ViewGrid>
      <CollectionInfo>
        <H3>{collection.name}</H3>
        <Description>{collection.description}</Description>
        <LinkContainer>{renderLinks(collection.links)}</LinkContainer>
      </CollectionInfo>

      <CollectionButtons>
        {isPrivate && <DeleteCollection id={id} />}
        <Button
          intent={Intent.PRIMARY}
          minimal
          outlined
          disabled
          text="Share"
        />
        {isPrivate && (
          <PublishCollection isPublishable={isPublishable} id={id} />
        )}
      </CollectionButtons>
      <DatasetContainer>
        {isDatasetPresent ? (
          <DatasetsGrid
            datasets={datasets}
            uploadedFiles={uploadedFiles}
            invalidateCollectionQuery={invalidateCollectionQuery}
            onSelect={setSelected}
          />
        ) : (
          <EmptyDatasets onUploadFile={addNewFile} />
        )}
      </DatasetContainer>
      {isDatasetPresent && (
        <ActionButtons
          collectionId={collection?.id}
          selectedDatasetId={selected}
          visibility={visibility}
          addNewFile={addNewFile}
        />
      )}
    </ViewGrid>
  );
};

export default Collection;
