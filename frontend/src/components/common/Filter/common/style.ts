import styled from "@emotion/styled";
import { Popover } from "@mui/material";
import { CommonThemeProps, getShadows, getSpaces } from "@czi-sds/components";

const shadowM = (props: CommonThemeProps) => getShadows(props)?.m;
const spacesS = (props: CommonThemeProps) => getSpaces(props)?.s;
const spacesXxs = (props: CommonThemeProps) => getSpaces(props)?.xxs;

export const Filter = styled.div`
  display: grid;
  font-feature-settings: normal; /* required; overrides layout.css specification */
  gap: ${spacesS}px;
  margin-bottom: ${spacesS}px;

  &:last-child {
    margin-bottom: 0;
  }
`;

export const FilterPopover = styled(Popover)`
  & .MuiPopover-paper {
    box-shadow: ${shadowM};
    font-feature-settings: normal; /* required; overrides layout.css specification */
  }
`;

export const InfoButtonWrapper = styled.span`
  cursor: pointer;
  margin-left: ${spacesXxs}px;
`;

export const NotificationWrapper = styled.div`
  position: absolute;
  z-index: 999;
  overflow: hidden;
`;
