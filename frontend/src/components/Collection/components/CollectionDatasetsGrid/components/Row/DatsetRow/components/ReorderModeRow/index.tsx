import { DragEvent, ReactNode, useCallback, useState } from "react";
import { Dataset } from "src/common/entities";
import ReorderCell from "src/components/common/Grid/components/ReorderCell";
import { Row } from "src/components/Collection/components/CollectionDatasetsGrid/components/Row/DatsetRow/components/ReorderModeRow/style";
import { ReorderAction } from "src/views/Collection/hooks/useReorderMode";
import { DragAndDropAction } from "src/views/Collection/hooks/useDragAndDrop/useDragAndDrop";
import { SerializedStyles } from "@emotion/react";
import { OffsetByIndex } from "src/views/Collection/hooks/useDragAndDrop/common/entities";
import {
  SELECTOR_TABLE_BODY,
  SELECTOR_TABLE_ROW,
} from "src/views/Collection/hooks/useDragAndDrop/common/constants";

export enum DRAG_EVENT_TYPE {
  DRAG = "DRAG",
  DRAG_START = "DRAG_START",
}

export type ReorderModeRowProps = Props;

interface Props {
  children: ReactNode | ReactNode[];
  dataset: Dataset;
  datasetIndex: number;
  dragAndDropAction: DragAndDropAction;
  draggingStyles: SerializedStyles;
  reorderAction: ReorderAction;
}

export default function ReorderModeRow({
  children,
  dataset,
  datasetIndex,
  dragAndDropAction,
  draggingStyles,
  reorderAction,
}: Props): JSX.Element {
  const { onReorder } = reorderAction;
  const {
    onDragging,
    onDraggingOver,
    onDropping,
    onEndDragging,
    onStartDragging,
  } = dragAndDropAction;
  const [dragEventType, setDragEventTypeState] = useState<DRAG_EVENT_TYPE>();

  // Row is being dragged.
  const onDrag = useCallback(
    (dragEvent: DragEvent<HTMLTableRowElement>) => {
      dragEvent.preventDefault();
      // Set drag event type state to DRAG as dragging is in progress;
      // now that the representation of the original row is generated,
      // the original row styles are updated to indicate the row is being dragged.
      setDragEventTypeState(DRAG_EVENT_TYPE.DRAG);
      // While dragging, update the direction of the drag event.
      onDragging(dragEvent);
    },
    [onDragging]
  );

  // Row dragging is successful or cancelled.
  const onDragEnd = useCallback(
    (dragEvent: DragEvent<HTMLTableRowElement>) => {
      dragEvent.preventDefault();
      dragEvent.dataTransfer.clearData();
      // Set drag event type state to undefined as dragging is now complete; row styles revert to their default setting.
      setDragEventTypeState(undefined);
      onEndDragging();
    },
    [onEndDragging]
  );

  // Row is dragging over a valid drop target.
  const onDragOver = useCallback(
    (dragEvent: DragEvent<HTMLTableRowElement>, droppingIndex: number) => {
      dragEvent.preventDefault(); // Enables drop container to receive drop events.
      dragEvent.dataTransfer.dropEffect = "move";
      // While dragging over a valid drop target, update the dropping index.
      onDraggingOver(droppingIndex);
    },
    [onDraggingOver]
  );

  // Row dragging has started; row element is referenced for dragging.
  const onDragStart = useCallback(
    (
      dragEvent: DragEvent<HTMLTableRowElement>,
      datasetID: string,
      draggingIndex: number
    ) => {
      dragEvent.dataTransfer.setData("text/plain", datasetID);
      dragEvent.dataTransfer.effectAllowed = "move";
      // Set drag event type state to DRAG_START as dragging has started; row styles are updated to facilitate
      // the image capture of the row being dragged. This image is generated from the drag target, and will represent
      // the draggable row as it is being dragged.
      setDragEventTypeState(DRAG_EVENT_TYPE.DRAG_START);
      // Initiate drag and drop process.
      onStartDragging(draggingIndex, getOffsetByIndex(dragEvent));
    },
    [onStartDragging]
  );

  // Row is dropped on a valid drop target.
  const onDrop = useCallback(
    (dragEvent: DragEvent<HTMLTableRowElement>) => {
      dragEvent.preventDefault();
      // Get the ID of the row being dragged.
      const datasetID = dragEvent.dataTransfer.getData("text");
      if (!datasetID) return; // Drop unsuccessful.
      // Update the order of the datasets.
      onDropping(onReorder);
    },
    [onDropping, onReorder]
  );

  return (
    <Row
      dragEventType={dragEventType}
      draggable
      draggingStyles={draggingStyles}
      onDrag={onDrag}
      onDragEnd={onDragEnd}
      onDragOver={(event) => onDragOver(event, datasetIndex)}
      onDragStart={(event) => onDragStart(event, dataset.id, datasetIndex)}
      onDrop={onDrop}
    >
      {/* Reordering handle */}
      <td>
        <ReorderCell />
      </td>
      {children}
    </Row>
  );
}

/**
 * Returns a map of row height by index.
 * Referenced by dragging styles state to calculate the y-axis translation of the dragging row and dropping rows
 * during drag operation.
 * @param dragEvent - Drag event.
 * @returns row height by index.
 */
function getOffsetByIndex(dragEvent: DragEvent<HTMLElement>): OffsetByIndex {
  const tbodyEl = dragEvent.currentTarget.closest(SELECTOR_TABLE_BODY);
  const rowsEl = tbodyEl?.querySelectorAll(SELECTOR_TABLE_ROW);
  const offsetByIndex = new Map();
  rowsEl?.forEach((rowEl, index) => {
    offsetByIndex.set(index, rowEl.getBoundingClientRect().height);
  });
  return offsetByIndex;
}
