import { DragEvent, ReactNode, useCallback, useState } from "react";
import { Dataset } from "src/common/entities";
import ReorderCell from "src/components/common/Grid/components/ReorderCell";
import { Row } from "src/components/Collection/components/CollectionDatasetsGrid/components/Row/DatsetRow/components/ReorderModeRow/style";
import { ReorderAction } from "src/common/hooks/useReorderMode";

export enum DRAG_EVENT_TYPE {
  DRAG = "DRAG",
  DRAG_START = "DRAG_START",
}

export type ReorderModeRowProps = Props;

interface Props {
  children: ReactNode | ReactNode[];
  dataset: Dataset;
  reorderAction: ReorderAction;
}

export default function ReorderModeRow({
  children,
  dataset,
  reorderAction,
}: Props): JSX.Element {
  const { onReorder } = reorderAction;
  const [dragEventType, setDragEventTypeState] = useState<DRAG_EVENT_TYPE>();

  // Row is being dragged.
  const onDrag = useCallback((dragEvent: DragEvent<HTMLTableRowElement>) => {
    dragEvent.preventDefault();
    // Set drag event type state to DRAG as dragging is in progress;
    // now that the representation of the original row is generated,
    // the original row styles are updated to indicate the row is being dragged.
    setDragEventTypeState(DRAG_EVENT_TYPE.DRAG);
  }, []);

  // Row dragging is successful or cancelled.
  const onDragEnd = useCallback((dragEvent: DragEvent<HTMLTableRowElement>) => {
    dragEvent.preventDefault();
    dragEvent.dataTransfer.clearData();
    // Set drag event type state to undefined as dragging is now complete; row styles revert to their default setting.
    setDragEventTypeState(undefined);
  }, []);

  // Row is dragging over a valid drop target.
  const onDragOver = useCallback(
    (dragEvent: DragEvent<HTMLTableRowElement>) => {
      dragEvent.preventDefault(); // Enables drop container to receive drop events.
      dragEvent.dataTransfer.dropEffect = "move";
    },
    []
  );

  // Row dragging has started; row element is referenced for dragging.
  const onDragStart = (
    dragEvent: DragEvent<HTMLTableRowElement>,
    datasetID: string
  ) => {
    dragEvent.dataTransfer.setData("text/plain", datasetID);
    dragEvent.dataTransfer.effectAllowed = "move";
    // Set drag event type state to DRAG_START as dragging has started; row styles are updated to facilitate
    // the image capture of the row being dragged. This image is generated from the drag target, and will represent
    // the draggable row as it is being dragged.
    setDragEventTypeState(DRAG_EVENT_TYPE.DRAG_START);
  };

  // Row is dropped on a valid drop target.
  const onDrop = useCallback(
    (dragEvent: DragEvent<HTMLTableRowElement>, targetDatasetID: string) => {
      dragEvent.preventDefault();
      // Calculate the position of the row being dropped relative to the drop target; the position is either
      // above or below the drop target mid-point.
      const targetEl = dragEvent.target as HTMLTableRowElement;
      const domRect = targetEl.getBoundingClientRect();
      const midPoint = domRect.top + domRect.height / 2;
      const orderPosition = Math.sign(dragEvent.clientY - midPoint);
      // Get the ID of the row being dragged.
      const datasetID = dragEvent.dataTransfer.getData("text");
      // Update the order of the datasets.
      onReorder(datasetID, targetDatasetID, orderPosition);
    },
    [onReorder]
  );

  return (
    <Row
      dragEventType={dragEventType}
      draggable
      onDrag={onDrag}
      onDragEnd={onDragEnd}
      onDragOver={onDragOver}
      onDragStart={(event) => onDragStart(event, dataset.id)}
      onDrop={(event) => onDrop(event, dataset.id)}
    >
      {/* Reordering handle */}
      <td>
        <ReorderCell />
      </td>
      {children}
    </Row>
  );
}
