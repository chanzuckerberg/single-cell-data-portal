import React, { ReactNode } from "react";
import { ROUTES } from "src/common/constants/routes";
import LinkCell from "src/components/common/Grid/components/LinkCell";
import { SubTitle } from "src/components/common/Grid/components/SubTitle";
import { Title } from "src/components/common/Grid/components/Title";

interface Props {
  children?: ReactNode;
  collectionId?: string;
  collectionName?: string;
  name: string;
}

export default function DatasetNameCell({
  children,
  collectionId = "",
  collectionName = "",
  name,
}: Props): JSX.Element {
  const url = ROUTES.COLLECTION.replace(":id", collectionId);
  return (
    <>
      <Title>{name}</Title>
      {!!collectionName &&
        (collectionId ? (
          <SubTitle>
            <LinkCell
              data-testid="collection-link"
              url={url}
              value={collectionName}
            />
          </SubTitle>
        ) : (
          <SubTitle>{collectionName}</SubTitle>
        ))}
      {children || null}
    </>
  );
}
