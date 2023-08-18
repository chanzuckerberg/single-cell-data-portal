import { FC } from "react";
import { UseMutateAsyncFunction } from "react-query";
import { Collection, Dataset, VISIBILITY_TYPE } from "src/common/entities";
import { ReuploadLink } from "src/common/queries/collections";
import DatasetRow from "src/components/Collection/components/CollectionDatasetsGrid/components/Row/DatsetRow";
import { RightAlignCell } from "src/components/common/Grid/components/RightAlignCell";
import { Grid as StyledGrid } from "src/components/common/Grid/style";
import { Props as ChooserProps } from "src/components/DropboxChooser/index";
import { sortByCellCountDescending } from "./common/util";
import { UploadedFiles } from "src/views/Collection/components/CollectionActions/components/AddButton";

interface Props {
  className?: string;
  collectionId: Collection["id"];
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
  collectionId,
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
            collectionId={collectionId}
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

export default DatasetsGrid;
