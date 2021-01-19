import { Button, Classes, H3, H6, Intent, Popover } from "@blueprintjs/core";
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
import {
  useCollection,
  useCollectionUploadLinks,
  USE_COLLECTION,
} from "src/common/queries/collections";
import { useDeleteDataset } from "src/common/queries/datasets";
import { getUrlHost } from "src/common/utils/getUrlHost";
import DatasetsGrid from "src/components/Collections/components/Grid/components/DatasetsGrid";
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
  StyledAlert,
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

  const [selected, setSelected] = useState<Dataset["id"]>();

  const [uploadedFiles, setUploadedFiles] = useState({} as UploadedFiles);

  const queryCache = useQueryCache();

  const { data: collection, isError } = useCollection(id, visibility);

  const [mutate] = useCollectionUploadLinks(id, visibility);

  const [isAlertOpen, openAlert] = useState(false);

  const [deleteDataset] = useDeleteDataset(collection?.id);

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

  const datasetPresent =
    collection.datasets?.length > 0 || Object.keys(uploadedFiles).length > 0;

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
      <DatasetContainer>
        {datasetPresent ? (
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
      {datasetPresent && (
        <StyledDiv>
          <DropboxChooser onUploadFile={addNewFile}>
            <Button intent={Intent.PRIMARY} outlined>
              Add
            </Button>
          </DropboxChooser>
          <Popover>
            <Button intent={Intent.PRIMARY} outlined>
              Download
            </Button>
          </Popover>
          <Button
            icon={IconNames.TRASH}
            minimal
            onClick={() => {
              if (selected) openAlert(!isAlertOpen);
            }}
          ></Button>
        </StyledDiv>
      )}
      <StyledAlert
        cancelButtonText={"Cancel"}
        confirmButtonText={"Delete Dataset"}
        intent={Intent.DANGER}
        isOpen={isAlertOpen}
        onCancel={() => {
          openAlert(!isAlertOpen);
        }}
        onConfirm={() => {
          deleteDataset(selected);
          openAlert(!isAlertOpen);
        }}
      >
        <H6>Are you sure you want to delete this dataset?</H6>
        <p>You cannot undo this action</p>
      </StyledAlert>
    </ViewGrid>
  );
};

export default Collection;
