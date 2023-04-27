import { allCellTypes } from "src/views/CellCards/components/CellCardSearchBar/fixture";

// This will eventually be replaced by a query to the backend
export const useCellTypes = () => {
  return allCellTypes;
};

export const useCellTypesById = () => {
  const cellTypes = useCellTypes();
  return cellTypes.reduce((acc, curr) => {
    const [key] = Object.keys(curr);
    const value = curr[key];
    acc[key] = value;
    return acc;
  }, {});
};
