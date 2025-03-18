import styled from "@emotion/styled";
import { FormControl as MFormControl } from "@mui/material";
import { fontBodyXxs } from "@czi-sds/components";

import { spacesS, spacesXl, spacesXxl, gray500 } from "src/common/theme";

export const FormControl = styled(MFormControl)`
  display: grid;
  grid-template-columns: 1fr 1fr;
  grid-gap: ${spacesXxl}px !important;
  .data-format-info label{
    display: flex;
  }
  .data-format-checkbox-group{
    display: flex;
    flex-direction: row;
    align-items: center;
  }
  .data-format-checkbox-group label{
    margin-bottom: 0;
  }
  .file-size{
    padding-left: ${spacesS}px;
    color: #6c6c6c;
  }

  .data-format-type{
    margin-bottom: ${spacesS}px;
  }
  .data-type-description{
    ${fontBodyXxs}
    padding-left: ${spacesXl}px;
    color: ${gray500};
    margin-bottom: 0;
  }
`;