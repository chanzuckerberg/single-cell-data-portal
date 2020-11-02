import React, { FC } from "react";
import { useCollections } from "src/common/queries/collections";
import Collection from "./components/Collection";
import Heading from "./components/Heading";

const Collections: FC = () => {
  const { isFetching, data: collections } = useCollections();

  if (isFetching && !collections) return <div>Loading collections...</div>;

  if (!collections) return <div>Sorry, we could not find any collections</div>;

  return (
    <>
      <h1>Datasets</h1>
      <p>
        The cellxgene data portal is a repository of public, explorable
        single-cell datasets. If you have a public dataset that you would like
        hosted on this site, please let us know at{" "}
        <a href="mailto:cellxgene@chanzuckerberg.com">
          cellxgene@chanzuckerberg.com
        </a>
        .
      </p>
      <Heading />
      {collections?.map((collection) => (
        <Collection id={collection.id} key={collection.id} />
      ))}
    </>
  );
};

export default Collections;
