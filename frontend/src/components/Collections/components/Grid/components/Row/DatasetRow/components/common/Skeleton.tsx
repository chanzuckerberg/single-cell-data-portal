import { Classes } from "@blueprintjs/core";
import React from "react";
import { PT_GRID_SIZE_PX } from "src/components/common/theme";
import styled from "styled-components";

export const Skeleton = styled(() => (
  <div className={Classes.SKELETON}>PLACEHOLDER_TEXT</div>
))`
  padding: 0 !important;
  padding-left: ${PT_GRID_SIZE_PX * 2}px !important;
`;
