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

export const Wrapper = styled.div`
  height: 200vh;
  width: 100%;
  display: flex;
  flex-direction: column;
  margin-left: 120px;
  margin-top: 44px;
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
