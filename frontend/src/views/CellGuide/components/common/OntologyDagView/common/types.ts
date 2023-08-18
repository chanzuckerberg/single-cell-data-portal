import { CellOntologyTreeResponse as TreeNode } from "src/common/queries/cellGuide";

export interface TreeNodeWithState extends TreeNode {
  isExpanded?: boolean;
  showAllChildren?: boolean;
  x0?: number;
  y0?: number;
  x?: number;
  y?: number;
}

export interface LegendProps {
  xPos: number;
  yPos: number;
}
