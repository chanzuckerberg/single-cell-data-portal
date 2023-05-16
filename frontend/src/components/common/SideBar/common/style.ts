import styled from "@emotion/styled";
import { CommonThemeProps, getColors } from "czifui";
import SideBar from "src/components/common/SideBar";

const grey300 = (props: CommonThemeProps) => getColors(props)?.gray[300];

export const ListViewSideBar = styled(SideBar)`
  box-shadow: inset -0.5px 0 0 ${grey300};
`;
