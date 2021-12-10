import React from "react";

interface Props {
  cellCount: number;
}

export default function CountCell({ cellCount }: Props): JSX.Element {
  const count = cellCount ? cellCount.toLocaleString() : "";
  return <div>{count}</div>;
}
