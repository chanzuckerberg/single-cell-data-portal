import { Icon } from "@blueprintjs/core";
import { LIGHT_GRAY, PRIMARY_BLUE } from "src/components/common/theme";
import styled, { css } from "styled-components";

export const Filter = styled.div`
  display: grid;
  gap: 8px;
  margin-bottom: 12px;

  &:last-child {
    margin-bottom: 0;
  }
`;

export const SelectionIcon = styled(Icon)`
  && {
    align-items: center;
    color: ${PRIMARY_BLUE};
    display: flex;
    height: 18px;
    justify-content: center;
    margin: 0 8px 0 0;

    > svg {
      height: 14px;
      width: 14px;
    }
  }
}
`;

export const scrollbar = css`
  &::-webkit-scrollbar {
    width: 4px;
  }

  &::-webkit-scrollbar-thumb {
    background-clip: content-box;
    background-color: ${LIGHT_GRAY.A};
    border-radius: 4px;
  }
`;
