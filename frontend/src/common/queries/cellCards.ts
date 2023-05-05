import { allCellTypes } from "src/views/CellCards/common/fixtures";

// This will eventually be replaced by a query to the backend
export const useCellTypes = () => {
  return allCellTypes;
};

export const useCellTypesById = () => {
  const cellTypes = useCellTypes();
  return cellTypes.reduce((acc, curr) => {
    const { id, label } = curr;
    acc[id] = label;
    return acc;
  }, {});
};
