import RightSideBar from "src/components/common/RightSideBar";
import styled from "@emotion/styled";

export const StyledRightSideBar = styled(RightSideBar)`
  /**
   * (thuang): This is to ensure the Marker Gene Right Side Bar is on top of the
   * Gene Search Bar when the viewport width is narrow enough that they overlap
   */
  z-index: 3;
`;
