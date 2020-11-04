import { Button, Classes, H3, H4, Intent, UL } from "@blueprintjs/core";
import { IconNames } from "@blueprintjs/icons";
import { RouteComponentProps } from "@reach/router";
import psl from "psl";
import React, { FC } from "react";
import { Link, LINK_TYPE } from "src/common/entities";
import { useCollection } from "src/common/queries/collections";
import {
  CenterAlignedDiv,
  CollectionButtons,
  CollectionInfo,
  DatasetContainer,
  Description,
  LinkContainer,
  StyledButton,
  ViewGrid,
} from "./style";

interface RouteProps {
  id?: string;
}

export type Props = RouteComponentProps<RouteProps>;

const getDomain = (url: string): string | null => {
  let hostname;
  if (url.indexOf("//") > -1) {
    hostname = url.split("/")[2];
  } else {
    hostname = url.split("/")[0];
  }
  const parsedURL = psl.parse(hostname);
  if (parsedURL.error) {
    console.error(parsedURL.input, parsedURL.error);
    return null;
  }
  return parsedURL.domain;
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
            script to do this. Learn More
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
    if (!Object.values(LINK_TYPE).includes(type)) return null;
    const domain = getDomain(url);
    if (domain)
      return (
        <React.Fragment key={`${type}+${url}`}>
          <span className={Classes.TEXT_MUTED}>{type}</span>
          <a href={url}>{domain}</a>
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
      <CollectionButtons>
        <StyledButton
          icon={IconNames.TRASH}
          text={"Delete"}
          minimal
          intent={Intent.DANGER}
        />
        <StyledButton text={"Share"} outlined intent={Intent.PRIMARY} />
        <StyledButton text={"Publish"} intent={Intent.PRIMARY} />
      </CollectionButtons>
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
