import React, { FC } from "react";
import { useCollection } from "src/common/queries/collections";
import Dataset from "../Dataset";

interface Props {
  id: string;
}

const Collection: FC<Props> = ({ id }) => {
  const { data: collection } = useCollection(id);

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
