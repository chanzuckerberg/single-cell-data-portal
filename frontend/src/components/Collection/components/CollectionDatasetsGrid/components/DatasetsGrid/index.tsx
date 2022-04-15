/* Copied from src/components/Collections/components/Grid/components/DatasetsGrid and modified as an intermediate
   upgrade to the collection datasets table while keeping the existing core datasets grid outside of the filter feature
   flag untouched. Once filter feature flag is removed, the existing core datasets grid can be deleted and replaced
   with this version. Ideally collection datasets table should be moved to react-table but requires changes to how
   dataset information (statuses etc) is calculated and rolled out. */
import { FC } from "react";
import { UseMutateAsyncFunction } from "react-query";
import { Collection, Dataset, VISIBILITY_TYPE } from "src/common/entities";
import { ReuploadLink } from "src/common/queries/collections";
import DatasetRow from "src/components/Collection/components/CollectionDatasetsGrid/components/Row/DatsetRow";
import { RightAlignCell } from "src/components/common/Grid/components/RightAlignCell";
import { Grid as StyledGrid } from "src/components/common/Grid/style";
import { Props as ChooserProps } from "src/components/DropboxChooser/index";
import { UploadedFiles } from "src/views/Collection/components/ActionButtons";

interface Props {
  className?: string;
  datasets: Dataset[];
  uploadedFiles: UploadedFiles;
  invalidateCollectionQuery: () => void;
  visibility: VISIBILITY_TYPE;
  accessType?: Collection["access_type"];
  isRevision: boolean;
  onUploadFile: (
    reuploadDataset?: UseMutateAsyncFunction<
      unknown,
      unknown,
      ReuploadLink,
      unknown
    >,
    datasetId?: string
  ) => ChooserProps["onUploadFile"];
  reuploadDataset: UseMutateAsyncFunction<
    unknown,
    unknown,
    ReuploadLink,
    unknown
  >;
}

const DatasetsGrid: FC<Props> = ({
  className,
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
    <StyledGrid className={className}>
      <thead>
        <tr>
          <th>Dataset</th>
          <th>Tissue</th>
          <th>Disease</th>
          <th>Assay</th>
          <th>Organism</th>
          <th>
            <RightAlignCell>Cells</RightAlignCell>
          </th>
          <th />
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
    </StyledGrid>
  );
};

export function sortByCellCountDescending(datasets: Dataset[]): Dataset[] {
  return (
    datasets?.sort((a, b) => (b.cell_count ?? 0) - (a.cell_count ?? 0)) || []
  );
}

export default DatasetsGrid;
