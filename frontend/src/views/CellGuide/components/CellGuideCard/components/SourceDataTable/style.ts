import styled from "@emotion/styled";
import { CommonThemeProps, Tag, fontBodyXs } from "@czi-sds/components";
import { gray100, spacesL, spacesM, spacesS } from "src/common/theme";

export const StyledTag = styled(Tag)`
  font-weight: 400;
  padding: 2px 8px;
  background-color: rgba(0, 0, 0, 0.05);
  & .MuiChip-label {
    color: #000000;
    ${fontBodyXs}
  }
`;
export const SourceDataTableWrapper = styled.div`
  max-width: calc(100vw - ${spacesL}px - ${spacesL}px);
`;

export const MobileSourceDataTableWrapper = styled.div`
  max-width: calc(100vw - ${spacesL}px - ${spacesL}px);
  display: flex;
  flex-direction: column;
  min-width: 320px;
`;

interface MobileSourceDataTableEntryProps extends CommonThemeProps {
  highlight: boolean;
}

export const MobileSourceDataTableEntry = styled.div<MobileSourceDataTableEntryProps>`
  row-gap: ${spacesS}px;
  padding-top: ${spacesM}px;
  padding-bottom: ${spacesM}px;
  padding-left: ${spacesM}px;
  display: flex;
  flex-direction: column;
  background-color: ${(props) => (props.highlight ? gray100(props) : "white")};
`;

export const MobileSourceDataTableEntryRow = styled.div`
  display: flex;
  flex-direction: row;
  flex-wrap: wrap;
`;
