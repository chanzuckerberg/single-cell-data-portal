import { GRAY } from "src/components/common/theme";
import styled from "styled-components";

export const StyledRow = styled.tr`
  align-content: center;
  vertical-align: middle;
`;

export const CollectionTitleText = styled.a`
  font-size: 14px;
  line-height: 18px;
  letter-spacing: -0.1px;
`;
export const DOIText = styled.div`
  color: ${GRAY.A};
  font-size: 12px;
`;

export const StyledCell = styled.td`
  vertical-align: middle !important;
`;

export const DetailsCell = styled(StyledCell)`
  color: ${GRAY.A};
`;
