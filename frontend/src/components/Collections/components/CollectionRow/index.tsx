import React, { FC } from "react";
import { useCollection } from "src/common/queries/collections";

interface Props {
  id: string;
}

const CollectionRow: FC<Props> = ({ id }) => {
  const { data: collection } = useCollection(id);

  if (!collection?.datasets) return null;

  return (
    <tr>
      <td>{collection.name}</td>
      <td>{collection.organs?.join(", ")}</td>
      <td>{collection.assays?.join(", ")}</td>
      <td>{collection.species?.join(", ")}</td>
      <td>{collection.cell_count}</td>
      <td>{collection.status}</td>
    </tr>
  );
};

export default CollectionRow;
