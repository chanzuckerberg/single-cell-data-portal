import styled from "@emotion/styled";
import { fontBodyXs } from "@czi-sds/components";
import { textSecondary } from "src/common/theme";
import { OLD_GRAY } from "src/components/common/theme";

export const Caption = styled.div`
  ${fontBodyXs}
  color: ${textSecondary};

  p {
    margin: 0;

    & + p {
      margin-top: 20px;
    }
  }

  b {
    font-weight: 500;
  }
`;

/**
 * @deprecated by Caption styles once "DOWNLOAD_UX" feature flag is removed.
 */
export const DownloadCaption = styled.div`
  color: ${OLD_GRAY.TIPS};
  font-size: 12px;
  margin-top: -3px;
  width: 490px;
`;
