import styled from "@emotion/styled";
import { FormControlLabel } from "@mui/material";

export const Wrapper = styled.div`
  padding: 0 10px;
`;

export const Dot = styled.span`
  border-radius: 50%;
  width: 12px;
  height: 12px;

  ${({ color }: { color: string }) => {
    return `
      background-color: ${color};
    `;
  }}
`;

export const Dots = styled.div`
  display: flex;
  gap: 8px;
  margin-bottom: 10px;
`;

export const ContentWrapper = styled.div`
  display: flex;
  position: relative;
`;

export const StyledFormControlLabel = styled(FormControlLabel)`
  position: absolute;
  left: -120px;
  top: -10px;
`;

export const LabelWrapper = styled.div`
  display: flex;
  justify-content: center;
  align-items: center;
  gap: 5px;
`;

export const FlexDiv = styled.div`
  display: flex;
`;

export const TooltipLink = styled("a")``;
