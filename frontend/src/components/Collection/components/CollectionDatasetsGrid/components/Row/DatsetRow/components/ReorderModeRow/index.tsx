import { DragEvent, ReactNode, useCallback, useState } from "react";
import ReorderCell from "src/components/common/Grid/components/ReorderCell";
import { Row } from "src/components/Collection/components/CollectionDatasetsGrid/components/Row/DatsetRow/components/ReorderModeRow/style";
import {
  DragAndDrop,
  ElementHeightByIndex,
} from "src/views/Collection/hooks/useDragAndDrop/common/entities";
import {
  SELECTOR_TABLE_BODY,
  SELECTOR_TABLE_ROW,
} from "src/views/Collection/hooks/useDragAndDrop/common/constants";
import { Reorder } from "src/views/Collection/hooks/useReorder/common/entities";

export enum DRAG_EVENT_TYPE {
  DRAG = "DRAG",
  DRAG_START = "DRAG_START",
}

export type ReorderModeRowProps = Props;

interface Props {
  children: ReactNode | ReactNode[];
  dragAndDrop: DragAndDrop;
  reorder: Reorder;
  testId: string;
}

export default function ReorderModeRow({
  children,
  dragAndDrop,
  reorder,
  testId,
}: Props): JSX.Element {
  const { reorderAction } = reorder;
  const { onReorder } = reorderAction;
  const { datasetIndex, dragAndDropAction } = dragAndDrop;
  const { onDraggingOver, onDropping, onEndDragging, onStartDragging } =
    dragAndDropAction;
  const [dragEventType, setDragEventTypeState] = useState<DRAG_EVENT_TYPE>();

  // Row is being dragged.
  const onDrag = useCallback((dragEvent: DragEvent) => {
    dragEvent.preventDefault();
    // Set drag event type state to DRAG as dragging is in progress;
    // now that the representation of the original row is generated,
    // the original row styles are updated to indicate the row is being dragged.
    setDragEventTypeState(DRAG_EVENT_TYPE.DRAG);
  }, []);

  // Row dragging is successful or cancelled.
  const onDragEnd = useCallback(
    (dragEvent: DragEvent) => {
      dragEvent.preventDefault();
      // Set drag event type state to undefined as dragging is now complete; row styles revert to their default setting.
      setDragEventTypeState(undefined);
      onEndDragging();
    },
    [onEndDragging]
  );

  // Row is dragging over a valid drop target.
  const onDragOver = useCallback(
    (dragEvent: DragEvent, dropTargetIndex: number) => {
      dragEvent.preventDefault(); // Enables drop container to receive drop events.
      dragEvent.dataTransfer.dropEffect = "move";
      // While dragging over a valid drop target, update the dragging state.
      onDraggingOver(dropTargetIndex, dragEvent.clientY);
    },
    [onDraggingOver]
  );

  // Row dragging has started; row element is referenced for dragging.
  const onDragStart = useCallback(
    (dragEvent: DragEvent, dragTargetIndex: number) => {
      dragEvent.dataTransfer.effectAllowed = "move";
      // Set drag event type state to DRAG_START as dragging has started; row styles are updated to facilitate
      // the image capture of the row being dragged. This image is generated from the drag target, and will represent
      // the draggable row as it is being dragged.
      setDragEventTypeState(DRAG_EVENT_TYPE.DRAG_START);
      // Initiate drag and drop process.
      onStartDragging(
        dragTargetIndex,
        dragEvent.clientY,
        getElementHeightByIndex(dragEvent)
      );
    },
    [onStartDragging]
  );

  // Row is dropped on a valid drop target.
  const onDrop = useCallback(
    (dragEvent: DragEvent) => {
      dragEvent.preventDefault();
      // Update the order of the datasets.
      onDropping(onReorder);
    },
    [onDropping, onReorder]
  );

  return (
    <Row
      data-testid={testId}
      dragEventType={dragEventType}
      draggable
      onDrag={onDrag}
      onDragEnd={onDragEnd}
      onDragOver={(event) => onDragOver(event, datasetIndex)}
      onDragStart={(event) => onDragStart(event, datasetIndex)}
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
 * Referenced by dragging styles state to calculate the y-axis translation of the dragged and dragged over rows.
 * @param dragEvent - Drag event.
 * @returns row height keyed by index.
 */
function getElementHeightByIndex(dragEvent: DragEvent): ElementHeightByIndex {
  const tbodyEl = dragEvent.currentTarget.closest(SELECTOR_TABLE_BODY);
  const rowsEl = tbodyEl?.querySelectorAll(SELECTOR_TABLE_ROW);
  const elementHeightByIndex = new Map();
  rowsEl?.forEach((rowEl, index) => {
    elementHeightByIndex.set(index, rowEl.getBoundingClientRect().height);
  });
  return elementHeightByIndex;
}
