import { PT_GRID_SIZE_PX } from "src/components/common/theme";
import styled from "styled-components";
import Input from "../../../common/Form/Input";

export const Form = styled.form`
  width: 580px;
  padding-top: ${PT_GRID_SIZE_PX}px;
`;

export const StyledInput = styled(Input)`
  /* Blank for ContactWrapper to target */
`;

export const ContactWrapper = styled.div`
  display: flex;

  ${StyledInput}:not(:last-child) {
    margin-right: ${2 * PT_GRID_SIZE_PX}px;
  }
`;
