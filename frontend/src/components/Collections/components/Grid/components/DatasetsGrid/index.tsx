import { FC } from "react";
import { MutateFunction } from "react-query";
import { Collection, Dataset, VISIBILITY_TYPE } from "src/common/entities";
import { ReuploadLink } from "src/common/queries/collections";
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
  onUploadFile: (
    reuploadDataset?: MutateFunction<string, unknown, ReuploadLink>,
    datasetId?: string
  ) => ChooserProps["onUploadFile"];
  reuploadDataset: MutateFunction<string, unknown, ReuploadLink>;
}

const DatasetsGrid: FC<Props> = ({
  datasets,
  uploadedFiles,
  invalidateCollectionQuery,
  visibility,
  accessType,
  isRevision,
  onUploadFile,
  reuploadDataset,
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
        {sortByCellCountDescending(datasets).map((dataset) => (
          <DatasetRow
            visibility={visibility}
            accessType={accessType}
            key={dataset.id}
            dataset={dataset}
            file={uploadedFiles[dataset.id]}
            invalidateCollectionQuery={invalidateCollectionQuery}
            revisionsEnabled={isRevision}
            onUploadFile={onUploadFile(reuploadDataset, dataset.id)}
          />
        ))}
      </tbody>
    </StyledCollectionsGrid>
  );
};

export function sortByCellCountDescending(datasets: Dataset[]): Dataset[] {
  return (
    datasets?.sort((a, b) => (b.cell_count ?? 0) - (a.cell_count ?? 0)) || []
  );
}

export default DatasetsGrid;
