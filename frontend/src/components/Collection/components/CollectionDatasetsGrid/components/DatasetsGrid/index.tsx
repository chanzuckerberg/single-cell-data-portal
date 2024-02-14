import { FC } from "react";
import { UseMutateAsyncFunction } from "react-query";
import { Collection, Dataset, VISIBILITY_TYPE } from "src/common/entities";
import { ReuploadLink } from "src/common/queries/collections";
import DatasetRow from "src/components/Collection/components/CollectionDatasetsGrid/components/Row/DatsetRow";
import { RightAlignCell } from "src/components/common/Grid/components/RightAlignCell";
import { Grid as StyledGrid } from "src/components/common/Grid/style";
import { Props as ChooserProps } from "src/components/DropboxChooser/index";
import { UploadedFiles } from "src/views/Collection/components/CollectionActions/components/AddButton";
import { useDragAndDrop } from "src/views/Collection/hooks/useDragAndDrop/useDragAndDrop";
import { Reorder } from "src/views/Collection/hooks/useReorder/common/entities";
import { getDragAndDrop } from "src/views/Collection/hooks/useDragAndDrop/common/utils";

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
  reorder: Reorder;
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
  reorder,
  reuploadDataset,
}) => {
  const dragAndDrop = useDragAndDrop();
  return (
    <StyledGrid className={className}>
      <thead>
        <tr>
          {reorder.isReorder && <th />}
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
        {datasets.map((dataset, index) => (
          <DatasetRow
            accessType={accessType}
            key={dataset.id}
            collectionId={collectionId}
            dataset={dataset}
            dragAndDrop={getDragAndDrop(dragAndDrop, index)}
            file={uploadedFiles[dataset.id]}
            invalidateCollectionQuery={invalidateCollectionQuery}
            onUploadFile={onUploadFile(reuploadDataset, dataset.id)}
            reorder={reorder}
            revisionsEnabled={isRevision}
            visibility={visibility}
          />
        ))}
      </tbody>
    </StyledGrid>
  );
};

export default DatasetsGrid;
