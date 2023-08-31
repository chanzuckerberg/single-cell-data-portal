import styled from "@emotion/styled";
import { CommonThemeProps, Tag, fontBodyXs } from "@czi-sds/components";
import { spacesM, spacesS, gray300 } from "src/common/theme";

export const StyledTag = styled(Tag)`
  ${fontBodyXs}
  & .MuiChip-label {
    color: #000000;
  }
`;

export const SourceDataTableWrapper = styled.div`
  max-width: calc(100vw - 32px);
`;

export const MobileSourceDataTableWrapper = styled.div`
  max-width: calc(100vw - 32px);
  display: flex;
  flex-direction: column;
  min-width: 320px;
`;

interface MobileSourceDataTableEntryProps extends CommonThemeProps {
  index: number;
}
export const MobileSourceDataTableEntry = styled.div<MobileSourceDataTableEntryProps>`
  row-gap: ${spacesS}px;
  padding-top: ${spacesM}px;
  padding-bottom: ${spacesM}px;
  display: flex;
  flex-direction: column;
  border-top: ${(props) =>
    props.index >= 1 ? `0.5px solid ${gray300(props)}` : "unset"};
`;

export const MobileSourceDataTableEntryRow = styled.div`
  display: flex;
  flex-direction: row;
  flex-wrap: wrap;
`;
