import { Button, Classes, H3, Intent } from "@blueprintjs/core";
import { IconNames } from "@blueprintjs/icons";
import { RouteComponentProps } from "@reach/router";
import { memoize } from "lodash-es";
import React, { FC, useState } from "react";
import { useQueryCache } from "react-query";
import {
  COLLECTION_LINK_TYPE_OPTIONS,
  Dataset,
  Link,
  VISIBILITY_TYPE,
} from "src/common/entities";
import { hasAssets } from "src/common/modules/datasets/selectors";
import {
  useCollection,
  useCollectionUploadLinks,
  USE_COLLECTION,
} from "src/common/queries/collections";
import { getUrlHost } from "src/common/utils/getUrlHost";
import DownloadDataset from "src/components/Collections/components/Dataset/components/DownloadDataset";
import DatasetsGrid from "src/components/Collections/components/Grid/components/DatasetsGrid";
import DeleteDataset from "src/components/Collections/components/Grid/components/Row/DatasetRow/components/DeleteDataset";
import DropboxChooser, { UploadingFile } from "src/components/DropboxChooser";
import { ViewGrid } from "../globalStyle";
import { StyledLink } from "./common/style";
import DatasetUploadToast from "./components/DatasetUploadToast";
import EmptyDatasets from "./components/EmptyDatasets";
import {
  CollectionInfo,
  DatasetContainer,
  Description,
  LinkContainer,
  StyledDiv,
} from "./style";

interface RouteProps {
  id?: string;
}

export type Props = RouteComponentProps<RouteProps>;

const renderLinks = (links: Link[]) => {
  return links?.map(({ link_url: url, link_type: type }) => {
    const linkTypeOption = COLLECTION_LINK_TYPE_OPTIONS[type];

    if (!linkTypeOption) return null;

    const urlHost = getUrlHost(url);

    const { text } = linkTypeOption;

    if (!urlHost) return null;

    return (
      <React.Fragment key={`${type}+${url}`}>
        <span className={Classes.TEXT_MUTED}>{text}</span>
        <StyledLink href={url}>{urlHost}</StyledLink>
      </React.Fragment>
    );
  });
};
export interface UploadedFiles {
  [datasetID: string]: UploadingFile;
}

const Collection: FC<Props> = ({ id = "" }) => {
  const isPrivate = window.location.pathname.includes("/private");
  const visibility = isPrivate
    ? VISIBILITY_TYPE.PRIVATE
    : VISIBILITY_TYPE.PUBLIC;

  const [selected, setSelected] = useState<Dataset["id"]>("");

  const [uploadedFiles, setUploadedFiles] = useState({} as UploadedFiles);

  const queryCache = useQueryCache();

  const { data: collection, isError } = useCollection(id, visibility);

  const [mutate] = useCollectionUploadLinks(id, visibility);

  const addNewFile = (newFile: UploadingFile) => {
    if (!newFile.link) return;

    const payload = JSON.stringify({ url: newFile.link });
    mutate(
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

  const isDatasetPresent =
    collection.datasets?.length > 0 || Object.keys(uploadedFiles).length > 0;

  const invalidateCollectionQuery = memoize(
    () => {
      queryCache.invalidateQueries([USE_COLLECTION, id, visibility]);
    },
    () => id + visibility
  );

  const selectedDataset = getSelectedDataset({
    datasets: collection.datasets,
    selectedId: selected,
  });

  return (
    <ViewGrid>
      <CollectionInfo>
        <H3>{collection.name}</H3>
        <Description>{collection.description}</Description>
        <LinkContainer>{renderLinks(collection.links)}</LinkContainer>
      </CollectionInfo>
      <DatasetContainer>
        {isDatasetPresent ? (
          <DatasetsGrid
            datasets={collection.datasets}
            uploadedFiles={uploadedFiles}
            invalidateCollectionQuery={invalidateCollectionQuery}
            onSelect={setSelected}
          />
        ) : (
          <EmptyDatasets onUploadFile={addNewFile} />
        )}
      </DatasetContainer>
      {isDatasetPresent && (
        <StyledDiv>
          <DropboxChooser onUploadFile={addNewFile}>
            <Button intent={Intent.PRIMARY} outlined>
              Add
            </Button>
          </DropboxChooser>
          <DownloadDataset
            isDisabled={!hasAssets(selectedDataset)}
            name={selectedDataset?.name || ""}
            dataAssets={selectedDataset?.dataset_assets || []}
            Button={DownloadButton}
          />
          <DeleteDataset id={selected} collectionId={collection?.id} />
        </StyledDiv>
      )}
    </ViewGrid>
  );
};

function DownloadButton({ ...props }) {
  return (
    <Button intent={Intent.PRIMARY} outlined {...props}>
      Download
    </Button>
  );
}

function getSelectedDataset({
  selectedId,
  datasets,
}: {
  selectedId: string;
  datasets: Dataset[];
}) {
  return datasets.find((dataset) => dataset.id === selectedId);
}

export default Collection;
