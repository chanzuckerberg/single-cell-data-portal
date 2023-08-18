import styled from "@emotion/styled";
import { fontBodyS, fontBodyXxxs } from "@czi-sds/components";
import {
  Content as CommonContent,
  Header as CommonHeader,
} from "../../common/style";
import { gray600 } from "src/common/theme";

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

export const ListSubheader = styled.li`
  ${fontBodyXxxs}

  margin-bottom: 4px !important;
  color: ${gray600};
`;

export const InfoText = styled.div`
  margin-bottom: 16px;
`;
