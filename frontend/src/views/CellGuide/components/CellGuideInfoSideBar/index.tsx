import React, { Dispatch, SetStateAction } from "react";
import { RightSidebarProperties } from "src/components/common/RightSideBar";
import { CellType } from "../common/OntologyDagView/common/types";
import Description from "../CellGuideCard/components/Description";
import MarkerGeneTables from "../CellGuideCard/components/MarkerGeneTables";

export interface CellGuideInfoBarProps extends RightSidebarProperties {
  cellInfoCellType: CellType;
  setGeneInfoGene: React.Dispatch<React.SetStateAction<string | null>>;
  setTooltipContent: Dispatch<
    SetStateAction<{
      title: string;
      element: JSX.Element;
    } | null>
  >;
  skinnyMode: boolean;
  selectedOrganName: string;
  selectedOrganId: string;
  organismName: string;
  selectedGene?: string;
  selectGene: (gene: string) => void;
}

function CellGuideInfoBar({
  cellInfoCellType,
  setTooltipContent,
  setGeneInfoGene,
  selectedOrganName,
  selectedOrganId,
  organismName,
  selectedGene,
  skinnyMode,
  selectGene,
}: CellGuideInfoBarProps): JSX.Element | null {
  return (
    <div>
      <Description
        cellTypeId={cellInfoCellType.cellTypeId}
        cellTypeName={cellInfoCellType.cellTypeName}
        skinnyMode={true}
        inSideBar
      />
      <MarkerGeneTables
        setTooltipContent={setTooltipContent}
        key={cellInfoCellType.cellTypeId}
        cellTypeId={cellInfoCellType.cellTypeId}
        setGeneInfoGene={setGeneInfoGene}
        cellTypeName={cellInfoCellType.cellTypeName}
        skinnyMode={skinnyMode}
        inSideBar
        organName={selectedOrganName}
        organId={selectedOrganId}
        organismName={organismName}
        selectedGene={selectedGene}
        selectGene={selectGene}
      />
    </div>
  );
}

export default React.memo(CellGuideInfoBar);
