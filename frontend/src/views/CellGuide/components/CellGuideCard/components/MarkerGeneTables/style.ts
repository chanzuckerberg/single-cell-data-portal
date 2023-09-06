import {
  CommonThemeProps,
  fontBodyS,
  fontBodyXs,
  fontBodyXxs,
} from "@czi-sds/components";
import styled from "@emotion/styled";
import {
  fontWeightSemibold,
  gray200,
  gray400,
  gray500,
  primary400,
  spacesL,
} from "src/common/theme";
import Link from "../common/Link";

const DIVIDER_WIDTH = 2;

export const StyledLink = styled(Link)`
  min-width: 2.25ch;
  max-width: 3.25ch;
`;

export const TableTitleOuterWrapper = styled.div`
  display: flex;
  justify-content: space-between;
  align-items: center;
  width: 100%;
  flex-wrap: wrap;
`;

export const ReferenceTooltipWrapper = styled.div`
  row-gap: ${spacesL}px;
  display: flex;
  flex-direction: column;
`;
export const NoWrapWrapper = styled.span`
  white-space: nowrap;
`;

export const PublicationLinkWrapper = styled.div`
  display: flex;
  column-gap: 2px;
`;

export const StyledHeadCellContent = styled.div`
  display: flex;
  flex-direction: row;
  justify-content: flex-end;
`;

export const StyledCellNumerical = styled.span`
  display: flex;
  flex-direction: row;
  justify-content: flex-end;
  padding-right: 12px;
`;

export const TableSelectorRow = styled.div`
  display: flex;
  flex-direction: row;
  gap: 16px;
  border-bottom: ${DIVIDER_WIDTH}px solid ${gray200};
`;

export const MarkerStrengthContainer = styled.div`
  ${fontBodyXxs}
`;

interface TableSelectorButtonProps extends CommonThemeProps {
  isActive: boolean;
}
export const TableSelectorButton = styled.button<TableSelectorButtonProps>`
  ${fontBodyS}
  position: relative;
  background: none;
  border: none;
  cursor: pointer;
  padding: 0;
  margin: 0;
  font-weight: ${fontWeightSemibold};
  color: ${(props) => `${props.isActive ? "#000000" : gray500(props)}`};

  &::after {
    content: "";
    position: absolute;
    bottom: -2px;
    left: 0;
    width: 100%;
    height: 2px;
    background-color: ${(props) =>
      `${props.isActive ? primary400(props) : gray200(props)}`};
    transition: background-color 0.3s ease;
  }
`;

export const MarkerGenePagination = styled.div`
  display: flex;
  justify-content: space-between;
  flex-wrap: wrap;
`;

export const MarkerGeneInfo = styled.div`
  ${fontBodyS}
  color: ${gray500};
  align-items: flex-start;
  display: inline-flex;
`;

export const MarkerGeneTooltipText = styled.div`
  ${fontBodyXs}
  font-weight: 500;
`;

export const MarkerGeneTooltipSubtext = styled.div`
  ${fontBodyXxs}
  color: ${gray400};
`;

export const MarkerGeneTableWrapper = styled.div`
  max-width: calc(100vw - ${spacesL}px - ${spacesL}px);
`;
