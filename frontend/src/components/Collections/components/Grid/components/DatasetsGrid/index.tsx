import React, { FC } from "react";
import { Collection, Dataset, VISIBILITY_TYPE } from "src/common/entities";
import {
  DatasetHeaderCell,
  LeftAlignedHeaderCell,
  RightAlignedHeaderCell,
  StyledCollectionsGrid,
} from "src/components/Collections/components/Grid/common/style";
import { Props as ChooserProps } from "src/components/DropboxChooser/index";
import { UploadedFiles } from "src/views/Collection/components/ActionButtons";
import DatasetRow from "../Row/DatasetRow";

interface Props {
  datasets: Dataset[];
  uploadedFiles: UploadedFiles;
  invalidateCollectionQuery: () => void;
  visibility: VISIBILITY_TYPE;
  accessType?: Collection["access_type"];
  isRevision: boolean;
  onUploadFile: ChooserProps["onUploadFile"];
}

const DatasetsGrid: FC<Props> = ({
  datasets,
  uploadedFiles,
  invalidateCollectionQuery,
  visibility,
  accessType,
  isRevision,
}) => {
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
            visibility={visibility}
            accessType={accessType}
            key={dataset.id}
            dataset={dataset}
            file={uploadedFiles[dataset.id]}
            invalidateCollectionQuery={invalidateCollectionQuery}
            isRevision={isRevision}
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
