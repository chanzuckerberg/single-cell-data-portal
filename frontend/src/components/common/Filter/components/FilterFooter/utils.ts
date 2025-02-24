import { CellTypeFooter } from "./components/CellTypeFooter";
import { TissueFooter } from "./components/TissueFooter";
import { DevelopmentStageFooter } from "./components/DevelopmentStageFooter";

export const FOOTER_COMPONENT_MAP: Record<string, React.FC<any>> = {
  tissueFooter: TissueFooter,
  cellTypeFooter: CellTypeFooter,
  developmentStageFooter: DevelopmentStageFooter,
};
