import React, { FC, useState } from "react";
import { Dataset } from "src/common/entities";
import {
  CollectionHeaderCell,
  LeftAlignedHeaderCell,
  RightAlignedHeaderCell,
  StyledCollectionsGrid,
} from "src/components/Collections/components/CollectionsGrid/style";
import { UploadingFile } from "src/components/DropboxChooser";
import DatasetRow from "./components/DatasetRow";

interface Props {
  datasets: Dataset[];
  uploadedFiles: Map<Dataset["id"], UploadingFile>;
}

const DatasetsGrid: FC<Props> = ({ datasets, uploadedFiles }) => {
  const [selected, setSelected] = useState(new Set());

  const handleSelect = (id: string) => {
    const newSelected = new Set(selected);
    if (newSelected.has(id)) newSelected.delete(id);
    else newSelected.add(id);
    setSelected(newSelected);
  };

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
          <DatasetRow
            key={dataset.id}
            dataset={dataset}
            checkHandler={handleSelect}
            file={uploadedFiles.get(dataset.id)}
          />
        ))}
        {/* There is change we'll have to render datasets which yet to have been populated in datasets, but are in uploadedFiles */}
      </tbody>
    </StyledCollectionsGrid>
  );
};

export default DatasetsGrid;
