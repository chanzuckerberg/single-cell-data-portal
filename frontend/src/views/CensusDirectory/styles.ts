import {
  Button,
  fontBodyS,
  fontCapsXxxs,
  fontHeaderL,
  fontHeaderXl,
  fontHeaderXxl,
} from "@czi-sds/components";
import styled from "@emotion/styled";
import { RadioGroup } from "@mui/material";
import {
  fontWeightBold,
  fontWeightMedium,
  fontWeightRegular,
  fontWeightSemibold,
  gray400,
  spacesDefault,
  spacesL,
  spacesXl,
  spacesXxs,
  textSecondary,
} from "src/common/theme";

export const Content = styled.div`
  box-sizing: content-box;
  padding: 0px 40px;
  display: flex;
  flex-direction: column;
  margin: 80px auto;
  max-width: 1200px;
`;

export const Header = styled.h1`
  ${fontHeaderXxl}
  margin-bottom: ${spacesDefault}px;
  font-weight: ${fontWeightBold};
`;

export const Paragraph = styled.p`
  ${fontBodyS}
  font-weight: ${fontWeightRegular};
  margin-bottom: 0;
`;

export const DirectoryDescription = styled(Paragraph)`
  margin-bottom: 80px;
`;

export const TierContainer = styled.div`
  margin-bottom: 120px;
`;

export const TierTitle = styled.h3`
  ${fontHeaderXl}
  margin-bottom: ${spacesDefault}px;
  font-weight: ${fontWeightSemibold};
`;

export const TierDescription = styled.p`
  ${fontBodyS}
  color: ${textSecondary};
  font-weight: ${fontWeightRegular};
  margin-bottom: 0;
`;

export const ProjectTitle = styled.h4`
  ${fontHeaderL}
  font-weight: ${fontWeightSemibold};
  margin-bottom: ${spacesDefault}px;
`;

export const ProjectSubmitter = styled.h4`
  ${fontBodyS}
  font-weight: ${fontWeightSemibold};
  margin-bottom: ${spacesDefault}px;
`;

export const ProjectDesctiption = styled(Paragraph)`
  max-width: 85ch;
`;

export const ProjectContainer = styled.div`
  display: flex;
  flex-direction: row;
  justify-content: space-between;
  margin-top: 20px;
`;
export const ProjectButtons = styled.div`
  display: flex;
  flex-direction: row;
  gap: ${spacesDefault}px;
`;
export const ProjectDetails = styled.div`
  display: flex;
  flex-direction: column;
`;
export const DetailsContainer = styled.div`
  display: flex;
  flex-direction: row;
  gap: ${spacesXl}px;
  margin-top: ${spacesL}px;
`;

export const StyledButton = styled(Button)`
  font-weight: ${fontWeightMedium};
  min-width: 80px;
`;

export const ItemContainer = styled.div`
  display: flex;
  flex-direction: column;
  gap: ${spacesXxs}px;
`;

export const ItemLabel = styled.div`
  ${fontCapsXxxs}
  font-weight: ${fontWeightSemibold};
  font-feature-settings:
    "clig" off,
    "liga" off;
  color: ${gray400};
`;

export const StyledRadioGroup = styled(RadioGroup)`
  padding: 8px;
`;
