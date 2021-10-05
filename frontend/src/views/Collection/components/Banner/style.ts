import { ORANGE, PT_GRID_SIZE_PX } from "src/components/common/theme";
import styled from "styled-components";

export const Wrapper = styled.div`
  display: flex;
  flex-direction: column;
  justify-content: center;
  align-items: left;
  padding: ${2 * PT_GRID_SIZE_PX}px;

  position: left;
  width: 768px;
  height: 127px;
  left: 200px;
  top: 80px;
  margin: 0 auto ${2 * PT_GRID_SIZE_PX}px auto;
  background: ${ORANGE.F};
`;

export const BannerHeader = styled.div`
  font-size: 16px;
  font-style: normal;
  font-weight: 600;
  line-height: 19px;
  letter-spacing: -0.14px;
  text-align: left;
  padding-bottom: ${PT_GRID_SIZE_PX}px;
  color: ${ORANGE.A};
`;

export const BannerWrapper = styled.div`
  font-size: 14px;
  font-style: normal;
  font-weight: 400;
  line-height: 18px;
  letter-spacing: -0.1px;
  text-align: left;
  color: ${ORANGE.A};
`;
