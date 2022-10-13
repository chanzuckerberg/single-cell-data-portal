import { Dataset } from "src/common/entities";

export function sortByCellCountDescending(datasets: Dataset[]): Dataset[] {
  return (
    datasets?.sort((a, b) => (b.cell_count ?? 0) - (a.cell_count ?? 0)) || []
  );
}
