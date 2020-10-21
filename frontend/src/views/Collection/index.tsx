import React, { FC } from "react";

interface Props {
  id: string;
}

const Collection: FC<Props> = ({ id }) => {
  return <>Collection {id}</>;
};

export default Collection;
