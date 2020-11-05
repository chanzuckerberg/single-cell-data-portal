import { Button, Classes, H3, H4, Intent, UL } from "@blueprintjs/core";
import { RouteComponentProps } from "@reach/router";
import React, { FC } from "react";
import { COLLECTION_LINK_TYPE_OPTIONS, Link } from "src/common/entities";
import { useCollection } from "src/common/queries/collections";
import {
  CenterAlignedDiv,
  CollectionInfo,
  DatasetContainer,
  Description,
  LinkContainer,
  StyledLink,
  ViewGrid,
} from "./style";

interface RouteProps {
  id?: string;
}

export type Props = RouteComponentProps<RouteProps>;

const getDomain = (url: string): string | null => {
  let result;

  try {
    result = new URL(url);
  } catch {
    return null;
  }

  return result.host;
};

const RenderEmptyDatasets = () => {
  return (
    <CenterAlignedDiv>
      <H4>No datasets uploaded</H4>
      <div>
        Before you begin uploading dataset files:
        <UL>
          <li>
            You must validate your dataset locally. We provide a local CLI
            script to do this. <StyledLink>Learn More</StyledLink>
          </li>
          <li>
            We only support adding datasets in the h5ad format at this time.
          </li>
        </UL>
      </div>
      <Button
        intent={Intent.PRIMARY}
        outlined
        text={"Add Dataset from Dropbox"}
      />
    </CenterAlignedDiv>
  );
};

const renderLinks = (links: Link[]) => {
  return links.map(({ url, type }) => {
    const linkTypeOption = COLLECTION_LINK_TYPE_OPTIONS[type];

    if (!linkTypeOption) return null;

    const domain = getDomain(url);

    const { text } = linkTypeOption;

    if (!domain) return null;

    return (
      <React.Fragment key={`${type}+${url}`}>
        <span className={Classes.TEXT_MUTED}>{text}</span>
        <StyledLink href={url}>{domain}</StyledLink>
      </React.Fragment>
    );
  });
};

const Collection: FC<Props> = ({ id }) => {
  const { data: collection } = useCollection(id ?? "");
  if (!collection) return null;

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
          collection.datasets.length > 0 &&
          // eslint-disable-next-line sonarjs/no-redundant-boolean
          false ? /* DATASETS VIEW GOES HERE */ null : (
            <RenderEmptyDatasets />
          )
        }
      </DatasetContainer>
    </ViewGrid>
  );
};

export default Collection;
