import { RouteComponentProps } from "@reach/router";
import React, { FC } from "react";

interface RouteProps {
  id: string;
}

export type Props = RouteComponentProps<RouteProps>;

const Collection: FC<Props> = ({ id }) => {
  return <>Collection {id}</>;
};

export default Collection;
