import { Button, Classes, H3, Intent } from "@blueprintjs/core";
import { IconNames } from "@blueprintjs/icons";
import { RouteComponentProps } from "@reach/router";
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
import { getUrlHost } from "src/common/utils/getUrlHost";
import DatasetsGrid from "src/components/Collections/components/DatasetsGrid";
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

const Collection: FC<Props> = ({ id = "" }) => {
  const isPrivate = window.location.pathname.includes("/private");
  const visibility = isPrivate
    ? VISIBILITY_TYPE.PRIVATE
    : VISIBILITY_TYPE.PUBLIC;

  const [uploadedFiles, setUploadedFiles] = useState(
    new Map<Dataset["id"], UploadingFile>()
  );

  const queryCache = useQueryCache();

  const { data: collection, isError } = useCollection(id, visibility);

  const [mutate] = useCollectionUploadLinks(id, visibility);

  const addNewFile = (newFile: UploadingFile) => {
    if (!newFile.link || newFile.id) return;

    const payload = JSON.stringify({ url: newFile.link });
    mutate(
      { collectionId: id, payload },
      {
        onSuccess: (data) => {
          newFile.id = data;
          if (!newFile.id) return;
          DatasetUploadToast.show({
            icon: IconNames.TICK,
            intent: Intent.PRIMARY,
            message:
              "Your file is being uploaded which will continue in the background, even if you close this window.",
          });
          setUploadedFiles(
            new Map(
              Array.from(uploadedFiles.entries()).concat([
                [newFile.id, newFile],
              ])
            )
          );
          queryCache.invalidateQueries(USE_COLLECTION);
        },
      }
    );
  };

  if (!collection || isError) return null;

  console.log(uploadedFiles);

  return (
    <ViewGrid>
      <CollectionInfo>
        <H3>{collection.name}</H3>
        <Description>{collection.description}</Description>
        <LinkContainer>{renderLinks(collection.links)}</LinkContainer>
      </CollectionInfo>

      <DatasetContainer>
        {collection?.datasets?.length > 0 ? (
          <DatasetsGrid
            datasets={collection.datasets}
            uploadedFiles={uploadedFiles}
          />
        ) : (
          <EmptyDatasets onUploadFile={addNewFile} />
        )}
      </DatasetContainer>
      {collection?.datasets?.length > 0 && (
        <StyledDiv>
          <DropboxChooser onUploadFile={addNewFile}>
            <Button intent={Intent.PRIMARY} outlined>
              Add
            </Button>
          </DropboxChooser>
          <Button intent={Intent.PRIMARY} outlined>
            Download
          </Button>
          <Button icon={IconNames.TRASH} minimal></Button>
        </StyledDiv>
      )}
    </ViewGrid>
  );
};

export default Collection;
