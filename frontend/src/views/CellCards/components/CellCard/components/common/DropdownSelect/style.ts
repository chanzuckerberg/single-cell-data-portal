import styled from "@emotion/styled";
import Select from "@mui/material/Select";

export const StyledSelect = styled(Select)`
  line-height: 1em;
  & > div[role="button"] {
    min-height: 1em;
  }
`;
