import { Callout } from "@blueprintjs/core";
import { BLUE, PT_GRID_SIZE_PX } from "src/components/common/theme";
import styled from "styled-components";

export const CollectionRevisionCallout = styled(Callout)`
  color: ${BLUE.A};
  margin-bottom: ${PT_GRID_SIZE_PX * 2}px;
  padding: ${PT_GRID_SIZE_PX}px ${PT_GRID_SIZE_PX * 1.5}px;
  width: fit-content; /* required; overrides specificity of bp3 callout width rule */
`;
