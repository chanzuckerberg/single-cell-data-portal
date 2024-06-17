export interface Organism {
  id: string;
  name: string;
}

export type ExcludeOverlappingCells =
  | "excludeOne"
  | "excludeTwo"
  | "retainBoth";
