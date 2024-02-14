import { FC } from "react";
import { UseMutateAsyncFunction } from "react-query";
import { Collection, Dataset, VISIBILITY_TYPE } from "src/common/entities";
import { ReuploadLink } from "src/common/queries/collections";
import DatasetRow from "src/components/Collection/components/CollectionDatasetsGrid/components/Row/DatsetRow";
import { RightAlignCell } from "src/components/common/Grid/components/RightAlignCell";
import { Grid as StyledGrid } from "src/components/common/Grid/style";
import { Props as ChooserProps } from "src/components/DropboxChooser/index";
import { UploadedFiles } from "src/views/Collection/components/CollectionActions/components/AddButton";
import { ReorderAction } from "src/views/Collection/hooks/useReorderMode/useReorderMode";
import { useDragAndDrop } from "src/views/Collection/hooks/useDragAndDrop/useDragAndDrop";

interface Props {
  className?: string;
  collectionId: Collection["id"];
  datasets: Dataset[];
  uploadedFiles: UploadedFiles;
  invalidateCollectionQuery: () => void;
  visibility: VISIBILITY_TYPE;
  accessType?: Collection["access_type"];
  isReorder: boolean;
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
  reorderAction: ReorderAction;
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
  isReorder,
  onUploadFile,
  reorderAction,
  reuploadDataset,
}) => {
  const { dragAndDropAction, dragAndDropStyles } = useDragAndDrop();
  return (
    <StyledGrid className={className}>
      <thead>
        <tr>
          {isReorder && <th />}
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
            datasetIndex={index}
            dragAndDropAction={dragAndDropAction}
            dragAndDropStyles={dragAndDropStyles}
            file={uploadedFiles[dataset.id]}
            invalidateCollectionQuery={invalidateCollectionQuery}
            isReorder={isReorder}
            onUploadFile={onUploadFile(reuploadDataset, dataset.id)}
            reorderAction={reorderAction}
            revisionsEnabled={isRevision}
            visibility={visibility}
          />
        ))}
      </tbody>
    </StyledGrid>
  );
};

export default DatasetsGrid;
