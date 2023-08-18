import { css } from "@emotion/css";
import styled from "@emotion/styled";
import { fontBodyXs, TooltipTable } from "@czi-sds/components";
import { gray500, spacesL, spacesXs } from "src/common/theme";

export const StyledTooltipTable = styled(TooltipTable)`
  display: flex;
  flex-direction: column;
  gap: ${spacesL}px;
  padding-bottom: ${spacesXs}px;

  > div {
    padding-top: 0 !important;
  }

  .MuiTableRow-root {
    vertical-align: baseline;
  }

  .MuiTable-root {
    margin-bottom: 0;
  }

  .MuiTableCell-alignLeft {
    ${fontBodyXs}
    width: 30%;
    white-space: nowrap;
  }

  .MuiTableCell-alignRight {
    color: ${gray500};
  }
`;

export const tooltipCss = css`
  margin: 0;
  margin-top: 20px;
  max-width: 500px !important;
`;
