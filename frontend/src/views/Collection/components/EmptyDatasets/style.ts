import { GRAY } from "src/components/common/theme";
import styled from "styled-components";

export const CenterAlignedDiv = styled.div`
  display: flex;
  flex-direction: column;
  width: 340px;
  justify-content: center;
  align-items: center;
  margin: auto auto;
  height: 229px;
  font-weight: 400;
  line-height: 18px;
  letter-spacing: -0.10000000149011612px;
  text-align: left;
  color: ${GRAY.A};
`;
