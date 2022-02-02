import Link from "next/link";
import React from "react";
import {
  CollectionMetadata as Metadata,
  MetadataLabel,
  MetadataValue,
} from "src/components/Collection/components/CollectionMetadata/style";

export interface CollectionMetadataLink {
  dataTestId: string;
  label: string;
  url: string;
  value: string;
}

interface Props {
  collectionMetadata: CollectionMetadataLink[];
}

export default function CollectionMetadata({
  collectionMetadata,
}: Props): JSX.Element {
  return (
    <Metadata>
      {collectionMetadata.map(({ dataTestId, label, url, value }) => (
        <React.Fragment key={label}>
          <MetadataLabel>{label}</MetadataLabel>
          <Link href={url} passHref>
            <MetadataValue
              data-test-id={dataTestId}
              href="passHref"
              rel="noopener"
              target="_blank"
            >
              {value}
            </MetadataValue>
          </Link>
        </React.Fragment>
      ))}
    </Metadata>
  );
}
