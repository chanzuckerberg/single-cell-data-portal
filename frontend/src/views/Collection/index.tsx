import { Classes, H3 } from "@blueprintjs/core";
import { RouteComponentProps } from "@reach/router";
import React, { FC, useEffect, useState } from "react";
import {
  COLLECTION_LINK_TYPE_OPTIONS,
  Link,
  VISIBILITY_TYPE,
} from "src/common/entities";
import {
  useCollection,
  useCollectionUploadLinks,
} from "src/common/queries/collections";
import { getUrlHost } from "src/common/utils/getUrlHost";
import { ViewGrid } from "../globalStyle";
import { StyledLink } from "./common/style";
import EmptyDatasets from "./components/EmptyDatasets";
import {
  CollectionInfo,
  DatasetContainer,
  Description,
  LinkContainer,
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

  const [uploadLink, setUploadLink] = useState("");

  const { data: collection, isError } = useCollection(id, visibility);

  const [mutate] = useCollectionUploadLinks(id, visibility);

  useEffect(() => {
    if (!uploadLink) return;

    const payload = JSON.stringify({ url: uploadLink });

    mutate({ collectionId: id, payload });
  }, [uploadLink, mutate, id]);

  if (!collection || isError) return null;

  return (
    <ViewGrid>
      <CollectionInfo>
        <H3>{collection.name}</H3>
        <Description>{collection.description}</Description>
        <LinkContainer>{renderLinks(collection.links)}</LinkContainer>
      </CollectionInfo>

      <DatasetContainer>
        {
          // eslint-disable-next-line no-constant-condition
          collection?.datasets?.length > 0 &&
          // eslint-disable-next-line sonarjs/no-redundant-boolean
          false ? /* DATASETS VIEW GOES HERE */ null : (
            <EmptyDatasets onSelectUploadLink={setUploadLink} />
          )
        }
      </DatasetContainer>
    </ViewGrid>
  );
};

export default Collection;
