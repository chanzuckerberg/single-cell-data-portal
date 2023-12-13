import { RightSidebarProperties } from "../../../../components/common/RightSideBar";
import { OntologyTerm } from "src/common/queries/wheresMyGene";
import { State } from "../../common/store";

export interface CellInfoBarProps extends RightSidebarProperties {
  cellInfoCellType: Exclude<State["cellInfoCellType"], null>;
  tissueInfo: OntologyTerm;
  generateGeneInfo: (gene: string) => void;
}
