import styled from "@emotion/styled";
import { fontBodyS, fontBodyXs } from "@czi-sds/components";
import {
  Content as CommonContent,
  Header as CommonHeader,
} from "../../common/style";
import { gray600, warning500 } from "src/common/theme";
import {
  DivTableCell,
  DivTableRow,
} from "src/views/WheresMyGeneV2/common/styles";

export const Wrapper = styled.div`
  overflow-y: scroll;
`;

export const Header = styled(CommonHeader)`
  ${fontBodyS}

  padding-bottom: 4px;
  color: ${gray600};
`;

export const Content = styled(CommonContent)`
  width: unset;
  padding: 0 16px;
`;

export const InfoText = styled.div`
  margin-bottom: 16px;
  color: black;
  ${fontBodyXs}
`;

export const StyledLink = styled.a`
  color: #0073ff;
  font-size: 14px;
  font-weight: 400;
`;

export const StyledSpan = styled.span`
  font-size: 14px;
  font-weight: 400;
  color: ${warning500};
`;

export const StyledDivTableRow = styled(DivTableRow)`
  line-height: 20px;
`;

export const StyledDivTableCell = styled(DivTableCell)`
  padding: 4px;
  max-width: 200px;
`;

export const StyledLabel = styled.div`
  width: 50px;
  height: 100%;
  padding-left: 8px;
  padding-right: 8px;
  padding-top: 2px;
  padding-bottom: 2px;
  background: rgba(0, 0, 0, 0.05);
  border-radius: 4px;
  justify-content: center;
  align-items: center;
  display: inline-flex;
  color: black;
  ${fontBodyXs}
`;

export const StyledValue = styled.div`
  display: inline-flex;
  align-items: center;
  color: black;
  padding: 2px;
  width: 100%;
`;
