import React, {
  Dispatch,
  FC,
  SetStateAction,
  useEffect,
  useState,
} from "react";
import { Dataset } from "src/common/entities";
import {
  DatasetHeaderCell,
  LeftAlignedHeaderCell,
  RightAlignedHeaderCell,
  StyledCollectionsGrid,
} from "src/components/Collections/components/Grid/common/style";
import { UploadedFiles } from "src/views/Collection";
import DatasetRow from "../Row/DatasetRow";
interface Props {
  datasets: Dataset[];
  uploadedFiles: UploadedFiles;
  invalidateCollectionQuery: () => void;
  onSelect: Dispatch<SetStateAction<string>>;
}

const DatasetsGrid: FC<Props> = ({
  datasets,
  uploadedFiles,
  invalidateCollectionQuery,
  onSelect,
}) => {
  const [selected, setSelected] = useState<Dataset["id"]>("");

  useEffect(() => {
    onSelect(selected);
  }, [selected, onSelect]);

  const handleSelect = (id: Dataset["id"]) => {
    if (id !== selected) {
      setSelected(id);
    }
  };

  return (
    <StyledCollectionsGrid bordered>
      <thead>
        <tr>
          <DatasetHeaderCell>Dataset</DatasetHeaderCell>
          <LeftAlignedHeaderCell>Tissue</LeftAlignedHeaderCell>
          <LeftAlignedHeaderCell>Assay</LeftAlignedHeaderCell>
          <LeftAlignedHeaderCell>Disease</LeftAlignedHeaderCell>
          <LeftAlignedHeaderCell>Organism</LeftAlignedHeaderCell>
          <RightAlignedHeaderCell>Cell Count</RightAlignedHeaderCell>
          <RightAlignedHeaderCell />
        </tr>
      </thead>
      <tbody>
        {sortByCreatedAtAscending(datasets).map((dataset) => (
          <DatasetRow
            key={dataset.id}
            dataset={dataset}
            file={uploadedFiles[dataset.id]}
            invalidateCollectionQuery={invalidateCollectionQuery}
            onSelect={handleSelect}
            selected={selected}
          />
        ))}
      </tbody>
    </StyledCollectionsGrid>
  );
};

function sortByCreatedAtAscending(datasets: Dataset[]): Dataset[] {
  return datasets?.sort((a, b) => a.created_at - b.created_at) || [];
}

export default DatasetsGrid;
