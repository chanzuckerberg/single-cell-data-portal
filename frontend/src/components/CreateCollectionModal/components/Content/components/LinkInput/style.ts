import { Button, Icon } from "@blueprintjs/core";
import Input from "src/components/common/Form/Input";
import { StyledDiv } from "src/components/common/Form/Input/style";
import { PT_GRID_SIZE_PX } from "src/components/common/theme";
import styled from "styled-components";

export const CollectionLink = styled.div`
  column-gap: 16px;
  display: grid;
  grid-template-columns: 160px 1fr 1fr; /* grid columns for collection link selector, name and url */
  position: relative; /* positions close collection link icon */
`;

export const CloseCollectionLinkIcon = styled(Icon)`
  cursor: pointer;
  padding: 1px;
  position: absolute;
  right: 0;
  top: 0;
`;

/**
 * @deprecated - supersede once filter feature flag is removed (#1718).
 */
export const IconWrapper = styled.div`
  position: relative;
  width: ${3.5 * PT_GRID_SIZE_PX}px;
`;

/**
 * @deprecated - supersede once filter feature flag is removed (#1718).
 */
export const StyledButton = styled(Button)`
  && {
    position: absolute;
    top: 20px;
    right: -${PT_GRID_SIZE_PX}px;
  }
`;

/**
 * @deprecated - supersede once filter feature flag is removed (#1718).
 */
export const StyledLinkTypeButton = styled(Button)`
  height: 34px;
  margin-top: 4px;
  justify-content: space-between;
`;

/**
 * @deprecated - supersede once filter feature flag is removed (#1718).
 */
export const StyledURLInput = styled(Input)`
  /* Blank for LinkWrapper to target */
`;

/**
 * @deprecated - supersede by CollectionLink once filter feature flag is removed (#1718).
 */
export const LinkWrapper = styled.div`
  display: flex;

  & > * :not(:last-child) {
    margin-right: ${2 * PT_GRID_SIZE_PX}px;
  }

  ${StyledDiv} {
    width: 25%;
  }
`;
