import styled from "@emotion/styled";
import { fontBodyXs } from "@czi-sds/components";
import { fontWeightRegular, gray500 } from "src/common/theme";

export const Wrapper = styled.div`
  display: flex;
  padding-bottom: 8px;
  align-items: center;
`;

export const Label = styled.div`
  ${fontBodyXs}

  color: ${gray500};
  font-weight: ${fontWeightRegular};
`;

export const Synonym = styled.div`
  ${fontBodyXs}

  padding: 4px 8px;
  color: black;
  font-weight: ${fontWeightRegular};
`;
