import styled from "@emotion/styled";
import { fontBodyXs, fontBodyS, Callout } from "@czi-sds/components";
import {
  fontWeightRegular,
  gray500,
  primary400,
  warning100,
  warning400,
} from "src/common/theme";

export const GeneSummary = styled.div`
  ${fontBodyXs}

  padding: 8px 0;
  font-weight: 400;
  color: black;
`;

export const Label = styled.div`
  ${fontBodyXs}

  color: ${gray500};
  font-weight: ${fontWeightRegular};
`;

export const Link = styled.a`
  ${fontBodyS}
  font-weight: 400;
  color: ${primary400};
`;

export const GeneName = styled.div`
  ${fontBodyS}
  font-weight: 500;
  color: black;
  padding-top: 4px;
`;

export const StyledCallout = styled(Callout)`
  ${fontBodyXs}
  align-items: center;
  width: 100%;
  margin-top: 8px;
  margin-bottom: 0;
`;

export const WarningBanner = styled.div`
  ${fontBodyXs}

  padding: 8px;
  display: flex;
  align-items: center;
  margin-top: 8px;

  background-color: ${warning100};
  svg {
    fill: ${warning400};
  }
`;

export const OutLinksWrapper = styled.div`
  display: flex;
  justify-content: space-between;
`;
