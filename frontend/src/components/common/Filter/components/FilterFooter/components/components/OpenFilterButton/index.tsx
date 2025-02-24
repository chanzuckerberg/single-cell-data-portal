import { LinkButton } from "../../../style";
import { ReactNode, useContext } from "react";
import { FilterControlContext } from "src/components/common/Filter/context";

export function OpenFilterButton({ children }: { children: ReactNode }) {
  const { setOpenSpecificFilter } = useContext(FilterControlContext);
  const openOrganismFilter = () => {
    setOpenSpecificFilter("Organism");
  };
  return <LinkButton onClick={openOrganismFilter}>{children}</LinkButton>;
}
