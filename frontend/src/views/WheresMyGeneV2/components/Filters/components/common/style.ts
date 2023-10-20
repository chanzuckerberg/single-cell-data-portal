import styled from "@emotion/styled";
import { Dropdown, fontBodyS, fontBodyXxxs } from "@czi-sds/components";
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
` as typeof Dropdown;

export const FilterLabel = styled("label")`
  ${fontBodyXxxs}

  color: ${gray500};
  font-weight: ${fontWeightMedium};
  margin-bottom: ${spacesXxs}px;
`;
