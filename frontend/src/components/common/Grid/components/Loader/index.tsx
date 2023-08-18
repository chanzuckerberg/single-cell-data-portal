import React from "react";
import { LoadingIndicator } from "@czi-sds/components";
import { Loader as GridLoader } from "./style";

export default function Loader(): JSX.Element {
  return (
    <GridLoader>
      <LoadingIndicator sdsStyle="tag" />
    </GridLoader>
  );
}
