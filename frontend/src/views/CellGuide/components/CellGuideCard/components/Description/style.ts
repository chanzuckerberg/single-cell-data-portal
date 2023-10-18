import styled from "@emotion/styled";
import {
  Button,
  CommonThemeProps,
  Tag,
  Tooltip,
  fontBodyS,
  fontBodyXs,
  fontBodyXxs,
  fontCapsXxs,
} from "@czi-sds/components";
import {
  gray100,
  gray500,
  gray400,
  fontWeightSemibold,
  spacesXxs,
  spacesS,
} from "src/common/theme";
import Synonyms from "src/components/Synonyms";
import OntologyId from "../../../OntologyId";

interface CellGuideCardDescriptionProps extends CommonThemeProps {
  inSideBar?: boolean;
}
export const CellGuideCardDescription = styled.div<CellGuideCardDescriptionProps>`
  ${fontBodyS}
  font-weight: 400;
  white-space: pre-wrap;
  ${(props) => {
    const { inSideBar } = props;
    const backgroundColor = inSideBar ? "unset" : gray100(props);
    const padding = inSideBar ? "0px" : "12px 16px 12px 16px";
    const borderRadius = inSideBar ? "0px" : "8px";
    return `
      background-color: ${backgroundColor};
      padding: ${padding};
      border-radius: ${borderRadius};
    `;
  }}
`;

export const StyledButton = styled(Button)`
  ${fontBodyS}
  font-weight: 500;
  text-transform: none;
`;

export const Wrapper = styled.div`
  margin-top: 8px;
`;

interface DescriptionWrapperProps extends CommonThemeProps {
  maxHeight?: number;
  inSideBar?: boolean;
}
export const DescriptionWrapper = styled.div<DescriptionWrapperProps>`
  ${(props) => `max-height: ${props.maxHeight}px;`}
  overflow: hidden;
  position: relative;
  ${(props) =>
    props.maxHeight &&
    !props.inSideBar &&
    `
    &:after {
      content: "";
      text-align: right;
      position: absolute;
      bottom: 0;
      right: 0;
      width: 100%;
      height: 40px;
      background: linear-gradient(
        to bottom,
        rgba(248, 248, 248, 0),
        rgba(248, 248, 248, 1) 100%
      );
    }
  `}
`;
export const Source = styled.div`
  ${fontBodyS}
  margin-top: 16px;
  display: flex;
  justify-content: space-between;
  gap: 20px;
  color: ${gray500};
`;

export const SourceLink = styled.div`
  ${fontBodyS}
  white-space: nowrap;
  align-self: end;
`;

export const DescriptionHeader = styled.div`
  ${fontCapsXxs}
  font-weight: ${fontWeightSemibold};
  color: ${gray500};
  margin-bottom: 8px;
`;

export const StyledTooltip = styled(Tooltip)`
  margin-right: 16px;
`;

export const ChatGptTooltipText = styled.div`
  ${fontBodyXs}
  font-weight: 500;
`;

export const ChatGptTooltipSubtext = styled.div`
  ${fontBodyXxs}
  color: ${gray400};
`;

export const StyledSynonyms = styled(Synonyms)`
  margin-top: ${spacesXxs}px;
`;

export const StyledOntologyId = styled(OntologyId)`
  margin-top: ${spacesXxs}px;
`;

export const FlexContainer = styled.div`
  display: flex;
  flex-direction: column;
`;

export const StyledTag = styled(Tag)`
  border-radius: 4px;
  padding: 2px 8px 2px 8px;
  gap: 6px;
  margin: auto 0;
  .MuiChip-label {
    ${fontBodyXs}
    font-weight: 500 !important;
    letter-spacing: -0.04px !important;
  }
`;

export const ValidatedWrapper = styled.div`
  margin-top: 16px;
`;

export const ValidatedInlineWrapper = styled.span`
  ${fontBodyS}
  display: inline;
  color: ${gray500};
`;

export const ReferencesWrapper = styled.div`
  display: flex;
  margin-top: ${spacesXxs}px;
  padding-bottom: ${spacesS}px;
  gap: ${spacesXxs}px;
  height: 36px;
  align-items: center;
`;
