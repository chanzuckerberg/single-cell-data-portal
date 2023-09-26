import React from "react";
import { LoadingIndicator } from "@czi-sds/components";

interface Props {
  className?: string;
  sdsStyle?: "minimal" | "tag";
}

export default function Loader({
  className,
  sdsStyle = "tag",
}: Props): JSX.Element {
  return (
    <div className={className}>
      <LoadingIndicator sdsStyle={sdsStyle} />
    </div>
  );
}
