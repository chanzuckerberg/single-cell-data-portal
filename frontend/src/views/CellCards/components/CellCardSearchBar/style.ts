import styled from "@emotion/styled";
import { Popper, PopperProps } from "@mui/material";
import { fontBodyXs, fontCapsXxxs, getColors } from "czifui";

interface StyledPopperProps extends PopperProps {
  width: number;
}

export const StyledPopper = styled(Popper)<StyledPopperProps>`
  z-index: 99;
  background-color: white;
  box-shadow: 0px 2px 4px rgba(0, 0, 0, 0.15), 0px 2px 10px rgba(0, 0, 0, 0.15);
  border-radius: 4px;
  max-height: 250px;
  overflow-y: scroll;
  padding: 8px 12px 8px 12px;

  ${(props) => {
    return `
      width: ${props.width}px;
      `;
  }}
`;

export const SectionTitle = styled.div`
  ${fontCapsXxxs}
  font-weight: 600;
  margin-bottom: 10px;

  ${(props) => {
    const colors = getColors(props);
    return `
      color: ${colors?.gray[500]};
      `;
  }}
`;

export const SectionItem = styled.div`
  ${fontBodyXs}
  font-weight: 400;
  margin-bottom: 12px;
  padding-left: 8px;
  cursor: pointer;
`;
