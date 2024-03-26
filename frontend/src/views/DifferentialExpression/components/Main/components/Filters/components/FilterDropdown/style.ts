import { fontBodyXs } from "@czi-sds/components";
import styled from "@emotion/styled";
import { primary400 } from "src/common/theme";

export const Tag = styled.div`
  display: inline-flex;
  align-items: center;
  padding: 4px 8px;
  background-color: ${primary400};
  color: white;
  border-radius: 4px;
  ${fontBodyXs}
  font-weight: 500;
  margin-right: 5px;
  cursor: default;
`;

export const CloseIcon = styled.span`
  display: inline-block;
  margin-left: 8px;
  cursor: pointer;
  &:after {
    content: "×";
  }
`;
