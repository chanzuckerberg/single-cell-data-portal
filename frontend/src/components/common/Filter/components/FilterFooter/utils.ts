import { CellTypeFooter } from "./CellTypeFooter";
import { TissueFooter } from "./TissueFooter";
import { DevelopmentStageFooter } from "./DevelopmentStageFooter";

export const FOOTER_COMPONENT_MAP: Record<string, React.FC<any>> = {
  tissueFooter: TissueFooter,
  cellTypeFooter: CellTypeFooter,
  developmentStageFooter: DevelopmentStageFooter,
};
