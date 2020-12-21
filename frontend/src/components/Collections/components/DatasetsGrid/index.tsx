import React, { FC } from "react";
import { Dataset } from "src/common/entities";
import {
  CollectionHeaderCell,
  LeftAlignedHeaderCell,
  RightAlignedHeaderCell,
  StyledCollectionsGrid,
} from "src/components/Collections/components/CollectionsGrid/style";
import DatasetRow from "./components/DatasetRow";

interface Props {
  datasets: Dataset[];
}

const DatasetsGrid: FC<Props> = ({ datasets }) => {
  return (
    <StyledCollectionsGrid bordered>
      <thead>
        <tr>
          <CollectionHeaderCell>Dataset</CollectionHeaderCell>
          <LeftAlignedHeaderCell>Tissue</LeftAlignedHeaderCell>
          <LeftAlignedHeaderCell>Assay</LeftAlignedHeaderCell>
          <LeftAlignedHeaderCell>Disease</LeftAlignedHeaderCell>
          <LeftAlignedHeaderCell>Organism</LeftAlignedHeaderCell>
          <RightAlignedHeaderCell>Cell Count</RightAlignedHeaderCell>
        </tr>
      </thead>
      <tbody>
        {datasets?.map((dataset) => (
          <DatasetRow key={dataset.id} dataset={dataset} />
        ))}
      </tbody>
    </StyledCollectionsGrid>
  );
};

export default DatasetsGrid;
