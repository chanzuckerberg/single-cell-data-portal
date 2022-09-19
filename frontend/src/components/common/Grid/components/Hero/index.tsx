import styled from "@emotion/styled";
import { GRAY, LIGHT_GRAY, PT_TEXT_COLOR } from "src/components/common/theme";

export const GridHero = styled.div`
  align-items: center;
  background-color: ${LIGHT_GRAY.E};
  border-radius: 8px;
  display: flex;
  flex-direction: column;
  padding: 80px;

  h3 {
    color: ${PT_TEXT_COLOR};
    font-size: 22px;
    font-weight: 500;
    letter-spacing: -0.23px;
    line-height: 25px;
    margin-bottom: 8px;
  }

  p {
    color: ${GRAY.A};
    letter-spacing: -0.1px;
    line-height: 18px;
  }

  > *:last-child {
    margin-bottom: 0;
  }
`;
