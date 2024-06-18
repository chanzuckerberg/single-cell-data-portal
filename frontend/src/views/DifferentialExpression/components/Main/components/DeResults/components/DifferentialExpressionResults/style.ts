import styled from "@emotion/styled";
import { Callout, fontBodyS, fontBodyXxs } from "@czi-sds/components";
import { gray100, primary400 } from "src/common/theme";
import { TextField } from "@mui/material";
import Image from "next/image";
import { inputBaseClasses } from "@mui/material/InputBase";
import { inputLabelClasses } from "@mui/material/InputLabel";
import { alertClasses } from "@mui/material/Alert";

const TABLE_WIDTH = "480px";
export const TableWrapper = styled.div`
  width: ${TABLE_WIDTH};
  [class*="StyledCell"] {
    max-width: ${parseFloat(TABLE_WIDTH.replace("px", "")) / 3}px;
    overflow: hidden;
    text-overflow: ellipsis;
    white-space: nowrap;
  }
`;

export const OpenInGE = styled.a`
  color: ${primary400};
  cursor: pointer;
  font-weight: 500;
  ${fontBodyS}
  display: flex;
  align-items: center;
  flex-direction: row;
  column-gap: 2px;
`;
export const StyledTextField = styled(TextField)`
  height: 32px;
  max-width: 140px;
  margin-top: 4px;
  & .${inputBaseClasses.root} {
    padding: 0;
    height: 32px;
  }
  & .${inputLabelClasses.root} {
    margin-top: -8px;
    z-index: 0;
  }
`;

export const TableHeaderWrapper = styled.div`
  display: flex;
  flex-direction: column;
`;

export const EffectSizeHeaderWrapper = styled.div`
  display: flex;
  flex-direction: row;
  column-gap: 8px;
  align-items: center;
  cursor: pointer;
`;

export const CellGroupWrapper = styled.div`
  min-height: 92px;
  height: fit-content;
  width: ${TABLE_WIDTH};
  margin-bottom: 8px;

  padding: 12px;
  background-color: ${gray100};
`;

export const CellGroupTitleWrapper = styled.span`
  display: flex;
  flex-direction: row;
  justify-content: space-between;
  align-items: center;
`;

export const CellGroupTitle = styled.span`
  ${fontBodyS}
  font-weight: 600;
`;

export const EffectSizeIndicator = styled.span`
  ${fontBodyXxs}
  color: #959595;
`;

export const FilterTagsWrapper = styled.div`
  display: flex;
  flex-direction: row;
  flex-wrap: wrap;
  column-gap: 8px;
  margin-top: 8px;
`;

export const StyledCallout = styled(Callout)`
  width: ${TABLE_WIDTH};

  .${alertClasses.icon} {
    margin-top: auto;
    margin-bottom: auto;
  }
`;

export const StyledIconImage = styled(Image)`
  width: 24px;
  height: 24px;

  /* CSS filter generated for #0073FF (primary400) using https://isotropic.co/tool/hex-color-to-css-filter/ */
  filter: invert(37%) sepia(58%) saturate(7350%) hue-rotate(205deg)
    brightness(105%) contrast(105%);
`;

export const StyledTooltipText = styled.div`
  text-align: left;
`;
