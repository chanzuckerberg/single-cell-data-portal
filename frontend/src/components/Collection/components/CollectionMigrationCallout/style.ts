import { Callout, Classes } from "@blueprintjs/core";
import styled from "@emotion/styled";
import { ORANGE } from "src/components/common/theme";

export const CollectionMigrationCallout = styled(Callout)`
  background-color: ${ORANGE.F};
  border-radius: unset; /* overrides bp callout border radius rule */
  color: ${ORANGE.A};
  letter-spacing: -0.1px;
  line-height: 18px;
  margin-bottom: 24px;
  padding: 16px;
  width: 600px;

  /* Title */
  .${Classes.HEADING} {
    color: inherit;
    letter-spacing: -0.144px;
    line-height: 19px;
  }

  *:last-child {
    margin-bottom: 0;
  }
`;
