import React from "react";
import { CollectionRevisionCallout } from "src/components/Collection/components/CollectionRevisionStatusCallout/style";
import Link from "next/link";
import { TextLink } from "./style";
import { ACCESS_TYPE, Collection } from "src/common/entities";

interface Props {
  collection: Collection;
}

export default function CollectionRevisionStatusCallout({
  collection,
}: Props): JSX.Element | null {
  const { access_type, revision_of, revising_in } = collection;
  return access_type === ACCESS_TYPE.WRITE && (revision_of || revising_in) ? (
    <CollectionRevisionCallout dismissible={false} sdsType="primary">
      <span data-testid="revision-status">
        {!!revising_in && (
          <span>
            This public collection has a pending revision.{" "}
            {/* (thuang): use `legacyBehavior` prop, since `<TextLink />` is `<a />` */}
            <Link href={`/collections/${revising_in}`} legacyBehavior passHref>
              <TextLink href="passHref">Continue Revision</TextLink>
            </Link>
          </span>
        )}
        {!!revision_of && (
          <span>
            This is a private revision of a published collection.{" "}
            {/* (thuang): use `legacyBehavior` prop, since `<TextLink />` is `<a />` */}
            <Link href={`/collections/${revision_of}`} legacyBehavior passHref>
              <TextLink href="passHref">Open Published Collection</TextLink>
            </Link>
          </span>
        )}
      </span>
    </CollectionRevisionCallout>
  ) : null;
}
