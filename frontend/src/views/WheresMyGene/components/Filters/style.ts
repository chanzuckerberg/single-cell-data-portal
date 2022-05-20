import styled from "@emotion/styled";
import { ComplexFilterInputDropdown } from "czifui";

export const StyledComplexFilterInputDropdown = styled(
  ComplexFilterInputDropdown
)`
  padding-left: 0;

  .MuiButton-label {
    margin-left: 0;
    margin-right: 10px;
  }

  &.Mui-disabled {
    border: 0;
  }
`;

export const Wrapper = styled("div")`
  display: flex;
  flex-direction: column;
  gap: 40px;
`;
