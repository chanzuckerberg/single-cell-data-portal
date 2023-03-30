import { Intent } from "@blueprintjs/core";
import React from "react";
import { CollectionRevisionCallout } from "src/components/Collection/components/CollectionRevisionStatusCallout/style";

interface Props {
  isRevisionDifferent: boolean;
}

export default function CollectionRevisionStatusCallout({
  isRevisionDifferent,
}: Props): JSX.Element {
  return (
    <CollectionRevisionCallout intent={Intent.PRIMARY} icon={null}>
      <span data-testid="revision-status">
        {isRevisionDifferent
          ? "This collection has changed since you last published it."
          : "This is a private revision of a public collection."}
      </span>
    </CollectionRevisionCallout>
  );
}
