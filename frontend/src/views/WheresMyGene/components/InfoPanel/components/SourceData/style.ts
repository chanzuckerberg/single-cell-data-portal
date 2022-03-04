import styled from "@emotion/styled";
import { fontBodyXxxs, getColors } from "czifui";
import {
  Content as CommonContent,
  Header as CommonHeader,
} from "../../common/style";

export const Wrapper = styled.div`
  overflow-y: scroll;
`;

export const Header = styled(CommonHeader)`
  padding-bottom: 4px;
  border-bottom: 1px solid rgba(16, 22, 26, 0.15);
`;

export const Content = styled(CommonContent)`
  width: unset;
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
