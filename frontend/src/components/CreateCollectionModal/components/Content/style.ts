import { InputGroup } from "@blueprintjs/core";
import { PT_GRID_SIZE_PX } from "src/components/common/theme";
import styled from "styled-components";

export const Form = styled.form`
  width: 580px;
  padding-top: ${PT_GRID_SIZE_PX}px;
`;

export const ContactWrapper = styled.div`
  display: flex;
`;

export const StyledInputGroup = styled(InputGroup)`
  margin-right: ${2 * PT_GRID_SIZE_PX}px;
`;
