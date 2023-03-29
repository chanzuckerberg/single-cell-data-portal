import Link from "next/link";
import React from "react";
import {
  CollectionMetadata as Metadata,
  MetadataLabel,
  MetadataValue,
} from "src/components/Collection/components/CollectionMetadata/style";

export interface CollectionMetadataLink {
  label: string;
  testId: string;
  url: string;
  value: string;
}

interface Props {
  collectionMetadataLinks: CollectionMetadataLink[];
}

export default function CollectionMetadata({
  collectionMetadataLinks,
}: Props): JSX.Element {
  return (
    <Metadata>
      {collectionMetadataLinks.map(({ label, testId, url, value }, i) => (
        <React.Fragment key={`${value}${i}`}>
          <MetadataLabel>{label}</MetadataLabel>
          <Link href={url} passHref>
            <MetadataValue
              data-testid={testId}
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
