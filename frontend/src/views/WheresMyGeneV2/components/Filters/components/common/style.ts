import styled from "@emotion/styled";
import { Dropdown, fontBodyS, fontBodyXxs } from "@czi-sds/components";
import {
  fontWeightMedium,
  fontWeightSemibold,
  gray500,
  spacesXxs,
} from "src/common/theme";

export const Label = styled("div")`
  ${fontBodyS}

  font-weight: ${fontWeightSemibold};
  margin-bottom: ${spacesXxs}px;
`;

export const Wrapper = styled("div")`
  display: flex;
  flex-direction: column;
`;

export const StyledDropdown = styled(Dropdown)`
  width: 100%;
  font-weight: 500;
  .styled-label {
    font-weight: 500;
  }
` as typeof Dropdown;

export const FilterLabel = styled("label")`
  ${fontBodyXxs}

  color: ${gray500};
  font-weight: ${fontWeightMedium};
  margin-bottom: ${spacesXxs}px;
`;
