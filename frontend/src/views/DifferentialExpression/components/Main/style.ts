import styled from "@emotion/styled";
import { Drawer } from "@blueprintjs/core";
import {
  fontHeaderXl,
  fontBodyS,
  getColors,
  CommonThemeProps,
  Button,
  fontCapsXxs,
} from "@czi-sds/components";

const LEFT_PANEL_WIDTH = "60vw";
const RIGHT_PANEL_WIDTH = "40vw";

export const TwoPanelLayout = styled.div`
  display: flex;
  flex-direction: row;

  .leftPanel {
    border-right: 1px solid #ccc; // Gray vertical line
    width: ${LEFT_PANEL_WIDTH};
    min-width: fit-content;
  }

  .rightPanel {
    width: ${RIGHT_PANEL_WIDTH};
    min-width: fit-content;
  }
`;

export const Wrapper = styled.div`
  display: flex;
  flex-direction: column;
  margin-left: 120px;
  margin-top: 24px;
`;

export const StepHeader = styled.div`
  ${fontHeaderXl}

  margin-bottom: 14px;
  margin-top: 6px;
`;

export const CellGroupTitle = styled.div`
  ${fontBodyS}
  font-weight: 600;
  margin-bottom: 8px;
  margin-top: 25px;
  display: flex;
  flex-direction: row;
  justify-content: space-between;
`;

export const CopyButtonWrapper = styled.div`
  ${fontCapsXxs}
  font-weight: 600;
  cursor: pointer;

  ${(props: CommonThemeProps) => {
    const colors = getColors(props);
    return `
      color: ${colors?.gray[400]}
    `;
  }}
`;

export const WordPop = styled.span<CommonThemeProps>`
  ${(props) => {
    const colors = getColors(props);
    return `
        color: ${colors?.primary[400]};
    `;
  }}
`;

export const StepSubHeader = styled.div`
  ${fontBodyS}

  margin-bottom: 32px;

  ${(props) => {
    const colors = getColors(props);
    return `
        color: ${colors?.gray[500]};
    `;
  }}
`;

export const RunButton = styled(Button)`
  ${fontBodyS}
  margin-bottom: 50px;
  margin-top: 59px;
  width: 256px;
`;

export const RunButtonWrapper = styled.div`
  display: flex;
  flex-direction: row;
  justify-content: flex-end;
`;
export const StyledSidebarDrawer = styled(Drawer)`
  .bp4-drawer-header {
    box-shadow: none;
  }
`;

export const FlexRow = styled.div`
  display: flex;
  flex-direction: row;
  align-items: center;
  column-gap: 16px;
`;

export const QuerySelectorWrapper = styled.div`
  width: 698px;
`;
