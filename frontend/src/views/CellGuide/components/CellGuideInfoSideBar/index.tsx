import React, { Dispatch, SetStateAction } from "react";
import { RightSidebarProperties } from "src/components/common/RightSideBar";
import { CellType } from "../common/OntologyDagView/common/types";
import Description from "../CellGuideCard/components/Description";
import MarkerGeneTables from "../CellGuideCard/components/MarkerGeneTables";
import { StyledButton } from "../CellGuideCard/components/Description/style";
import { useRouter } from "next/router";
import { ROUTES } from "src/common/constants/routes";
import { MarkerGeneTableWrapper } from "./style";
import {
  CELLGUIDE_VIEW_PAGE_SIDEBAR_BUTTON_TEST_ID,
  CELLGUIDE_INFO_SIDEBAR_TEST_ID,
} from "./constants";
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
  const router = useRouter();
  return (
    <div data-testid={CELLGUIDE_INFO_SIDEBAR_TEST_ID}>
      <StyledButton
        data-testid={CELLGUIDE_VIEW_PAGE_SIDEBAR_BUTTON_TEST_ID}
        sdsType="primary"
        sdsStyle="minimal"
        onClick={() => {
          router.push(
            `${ROUTES.CELL_GUIDE}/${cellInfoCellType.cellTypeId.replace(
              ":",
              "_"
            )}`
          );
        }}
      >
        View CellGuide Page
      </StyledButton>
      <Description
        cellTypeId={cellInfoCellType.cellTypeId}
        cellTypeName={cellInfoCellType.cellTypeName}
        skinnyMode={true}
        inSideBar
      />
      <MarkerGeneTableWrapper>
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
      </MarkerGeneTableWrapper>
    </div>
  );
}

export default React.memo(CellGuideInfoBar);
