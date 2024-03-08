import { UseReorder } from "src/views/Collection/hooks/useReorder/useReorder";

export interface Reorder extends UseReorder {
  disabled: boolean;
  startReorder: () => void;
}
