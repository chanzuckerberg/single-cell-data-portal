import { Button, Intent } from "@blueprintjs/core";
import React from "react";
import { Dataset } from "src/common/entities";

export function DownloadButton({ ...props }) {
  return (
    <Button intent={Intent.PRIMARY} outlined {...props}>
      Download
    </Button>
  );
}

export function getSelectedDataset({
  selectedId,
  datasets,
}: {
  selectedId: string;
  datasets?: Dataset[];
}) {
  return datasets?.find((dataset) => dataset.id === selectedId);
}
