import { InputDropdownProps as IInputDropdownProps } from "@czi-sds/components";
import { SORT_BY } from "src/views/WheresMyGeneV2/common/types";

export const DEFAULT_INPUT_DROPDOWN_PROPS: Partial<IInputDropdownProps> = {
  sdsStyle: "square",
};

export const CELL_TYPE_OPTIONS = [
  { id: SORT_BY.CELL_ONTOLOGY, name: "Cell Ontology" },
  { id: SORT_BY.H_CLUSTER, name: "Hierarchical" },
];

export const GENE_OPTIONS = [
  { id: SORT_BY.USER_ENTERED, name: "As Entered" },
  { id: SORT_BY.H_CLUSTER, name: "Hierarchical" },
];
