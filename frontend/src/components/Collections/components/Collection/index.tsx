import React, { FC } from "react";
import { VISIBILITY_TYPE } from "src/common/entities";
import { useCollection } from "src/common/queries/collections";
import Dataset from "../Dataset";

interface Props {
  id: string;
  visibility: VISIBILITY_TYPE;
}

const Collection: FC<Props> = ({ id, visibility }) => {
  const { data: collection } = useCollection(id, visibility);

  if (!collection?.datasets) return null;

  return (
    <>
      {collection.datasets
        .map((dataset) => ({
          dataset,
          links: collection.links,
        }))
        .map(({ dataset, links }) => {
          return <Dataset key={dataset.id} dataset={dataset} links={links} />;
        })}
    </>
  );
};

export default Collection;
