import React, { Dispatch, SetStateAction } from "react";
import { RightSidebarProperties } from "src/components/common/RightSideBar";
import { CellType } from "../common/OntologyDagView/common/types";
import Description from "../CellGuideCard/components/Description";
import MarkerGeneTables from "../CellGuideCard/components/MarkerGeneTables";
import { ROUTES } from "src/common/constants/routes";
import { MarkerGeneTableWrapper, StyledLink } from "./style";
import {
  CELLGUIDE_VIEW_PAGE_SIDEBAR_BUTTON_TEST_ID,
  CELLGUIDE_INFO_SIDEBAR_TEST_ID,
} from "./constants";
import { track } from "src/common/analytics";
import { EVENTS } from "src/common/analytics/events";
import { useRouter } from "next/router";

export interface CellGuideInfoBarProps extends RightSidebarProperties {
  cellInfoCellType: CellType;
  setGeneInfoGene: React.Dispatch<React.SetStateAction<string | null>>;
  setCellInfoCellType: React.Dispatch<React.SetStateAction<CellType | null>>;
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
  setCellInfoCellType,
}: CellGuideInfoBarProps): JSX.Element | null {
  const router = useRouter();
  return (
    <div data-testid={CELLGUIDE_INFO_SIDEBAR_TEST_ID}>
      <StyledLink
        data-testid={CELLGUIDE_VIEW_PAGE_SIDEBAR_BUTTON_TEST_ID}
        href={`${ROUTES.CELL_GUIDE}/${cellInfoCellType.cellTypeId.replace(
          ":",
          "_"
        )}`}
        onClick={(e) => {
          if (!e.metaKey && !e.ctrlKey) {
            e.preventDefault();
            const href = e.currentTarget.getAttribute("href");
            if (!href) return;
            router.push(href);
          }
          track(EVENTS.CG_VIEW_CELLGUIDE_PAGE_CLICKED, {
            cell_type: cellInfoCellType.cellTypeName,
          });
        }}
      >
        View CellGuide Page
      </StyledLink>
      <Description
        cellTypeId={cellInfoCellType.cellTypeId}
        cellTypeName={cellInfoCellType.cellTypeName}
        skinnyMode={true}
        setTooltipContent={setTooltipContent}
        inSideBar
      />
      <MarkerGeneTableWrapper>
        <MarkerGeneTables
          setTooltipContent={setTooltipContent}
          key={cellInfoCellType.cellTypeId}
          cellTypeId={cellInfoCellType.cellTypeId}
          setGeneInfoGene={setGeneInfoGene}
          setCellInfoCellType={setCellInfoCellType}
          cellTypeName={cellInfoCellType.cellTypeName}
          skinnyMode={skinnyMode}
          cellInfoCellType={cellInfoCellType}
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
