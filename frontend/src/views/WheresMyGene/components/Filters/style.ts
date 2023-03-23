import styled from "@emotion/styled";
import { ComplexFilter, ComplexFilterInputDropdown } from "czifui";

export const StyledComplexFilterInputDropdown = styled(
  ComplexFilterInputDropdown
)`
  padding-left: 0;
  /* (thuang): Sex filter is short and doesn't need the default 64px width */
  min-width: 0;

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

export const StyledComplexFilter = styled(ComplexFilter)`
  width: 100%;
  margin-bottom: 16px;
` as typeof ComplexFilter;
