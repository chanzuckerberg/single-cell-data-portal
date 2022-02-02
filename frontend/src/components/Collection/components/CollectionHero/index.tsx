import { H3 } from "@blueprintjs/core";
import React, { ReactElement } from "react";
import { CollectionHero as Hero } from "src/components/Collection/components/CollectionHero/style";

interface Props {
  CollectionActionButtons: ReactElement;
  collectionName: string;
  DeleteCollectionButton: ReactElement;
  shouldShowPrivateWriteAction: boolean;
  shouldShowPublicWriteAction: boolean;
}

export default function CollectionHero({
  CollectionActionButtons,
  collectionName,
  DeleteCollectionButton,
  shouldShowPrivateWriteAction,
  shouldShowPublicWriteAction,
}: Props): JSX.Element {
  return (
    <Hero>
      <H3 data-test-id="collection-name">{collectionName}</H3>
      {shouldShowPrivateWriteAction && CollectionActionButtons}
      {shouldShowPublicWriteAction && DeleteCollectionButton}
    </Hero>
  );
}
