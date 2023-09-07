import styled from "@emotion/styled";
import { CommonThemeProps, fontHeaderXxl, Tag } from "@czi-sds/components";
import { DEFAULT_ONTOLOGY_WIDTH } from "../common/OntologyDagView/common/constants";
import { spacesL, spacesXxl } from "src/common/theme";

export const TOP_PADDING_PX = 32;
export const LEFT_RIGHT_PADDING_PX_XXL = 40;
const TISSUE_CARD_MAX_WIDTH = 1440;

interface TissueViewProps extends CommonThemeProps {
  skinnyMode: boolean;
}

export const TissueCardView = styled.div<TissueViewProps>`
  display: flex;
  flex-direction: row;
  max-width: 100vw;
  ${(props) => {
    const { skinnyMode } = props;
    const space = skinnyMode ? spacesL(props) : spacesXxl(props);
    return `
    padding: ${TOP_PADDING_PX}px ${space}px 0px
      ${space}px;
    `;
  }}
`;

export const Wrapper = styled.div<TissueViewProps>`
  display: flex;
  flex-direction: column;
  align-self: stretch;
  overflow-x: hidden;
  width: 100%;
  ${(props) => {
    const { skinnyMode } = props;
    const maxWidth = skinnyMode
      ? `${DEFAULT_ONTOLOGY_WIDTH}px`
      : `${TISSUE_CARD_MAX_WIDTH}px`;
    return `
    max-width: ${maxWidth};
    `;
  }}
`;

export const TissueCardWrapper = styled.div<TissueViewProps>`
  margin: 0 auto 80px;
  width: ${(props) =>
    props.skinnyMode
      ? `${DEFAULT_ONTOLOGY_WIDTH}px`
      : `${TISSUE_CARD_MAX_WIDTH}px`};
`;

export const TissueCardHeaderInnerWrapper = styled.div`
  display: flex;
  column-gap: 8px;
  align-items: center;
`;

export const DescriptionWrapper = styled.div`
  margin-bottom: 0;
`;

export const TissueCardHeader = styled.div`
  display: flex;
  column-gap: 8px;
  flex-direction: row;
  align-items: center;
  justify-content: space-between;
`;

export const TissueCardName = styled.div`
  ${fontHeaderXxl}
  font-weight: 700;
`;

export const StyledTag = styled(Tag)`
  height: 24px;
  margin: 0;
`;

export const SearchBarWrapper = styled.div`
  margin-bottom: 16px;
  margin-top: 16px;
  width: 240px;
`;
