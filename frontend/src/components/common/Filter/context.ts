import { createContext } from "react";

export interface FilterControlContextProps {
  openSpecificFilter: string | null;
  setOpenSpecificFilter: (filter: string | null) => void;
}

const originalContext = {
  openSpecificFilter: null,
  setOpenSpecificFilter: (filter: string | null) => {
    console.log("setOpenSpecificFilter not implemented", filter);
  },
};

export const FilterControlContext =
  createContext<FilterControlContextProps>(originalContext);
