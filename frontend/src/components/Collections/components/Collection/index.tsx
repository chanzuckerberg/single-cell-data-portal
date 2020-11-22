import React, { FC } from "react";
import { useCollection, VISIBILITY } from "src/common/queries/collections";
import Dataset from "../Dataset";

interface Props {
  id: string;
  visibility: VISIBILITY;
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
