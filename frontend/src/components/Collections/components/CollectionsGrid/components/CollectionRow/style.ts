import styled from "styled-components";
import { DetailsCell } from "./components/common/style";

export const StyledRow = styled.tr`
  align-content: center;
  vertical-align: middle;
`;

// Collection Title Column
export const CollectionTitleText = styled.a`
  font-size: 14px;
  line-height: 18px;
  letter-spacing: -0.1px;
`;

export const RightAlignedDetailsCell = styled(DetailsCell)`
  text-align: right !important;
`;
