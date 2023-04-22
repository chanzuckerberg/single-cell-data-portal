/**
 * Grid column configuration, extends React Table's ColumnInstance.
 */
export interface GridColumnProps {
  alignment?: ALIGNMENT;
  showCountAndTotal?: boolean;
}

/**
 * Table header and data cell alignment.
 */
export enum ALIGNMENT {
  LEFT = "LEFT",
  RIGHT = "RIGHT",
}
