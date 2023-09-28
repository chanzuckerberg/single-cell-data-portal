import React from "react";
import { LoadingIndicator } from "@czi-sds/components";

interface Props {
  sdsStyle?: "minimal" | "tag";
}

export default function Loader({
  sdsStyle = "tag",
  ...props /* Spread props to allow for data-testid and other Loader props. */
}: Props): JSX.Element {
  return (
    <div {...props}>
      <LoadingIndicator sdsStyle={sdsStyle} />
    </div>
  );
}
