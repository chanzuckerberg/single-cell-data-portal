import styled from "@emotion/styled";
import { fontBodyS, fontBodyXxxs, getColors } from "@czi-sds/components";
import {
  Content as CommonContent,
  Header as CommonHeader,
} from "../../common/style";

export const Wrapper = styled.div`
  overflow-y: scroll;
`;

export const Header = styled(CommonHeader)`
  ${fontBodyS}

  padding-bottom: 4px;

  ${(props) => {
    const colors = getColors(props);

    return `
      color: ${colors?.gray[600]};
    `;
  }}
`;

export const Content = styled(CommonContent)`
  width: unset;
  padding: 0 16px;
`;

export const ListSubheader = styled.li`
  ${fontBodyXxxs}

  margin-bottom: 4px !important;

  ${(props) => {
    const colors = getColors(props);

    return `
      color: ${colors?.gray[600]};
    `;
  }}
`;
