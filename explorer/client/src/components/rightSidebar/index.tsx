/* Core dependencies */
import React, { CSSProperties } from "react";
import { connect } from "react-redux";

/* App dependencies */
import { RootState } from "../../reducers";
import GeneExpression from "../geneExpression";
import * as globals from "../../globals";

/* Styles */
export const STYLE_RIGHT_SIDEBAR: CSSProperties = {
  /* x y blur spread color */
  borderLeft: `1px solid ${globals.lightGrey}`,
  display: "flex",
  flexDirection: "column",
  height: "inherit",
  overflowY: "inherit",
  padding: globals.leftSidebarSectionPadding,
  position: "relative",
  width: "inherit",
};

const RightSidebar = () => (
  <div style={STYLE_RIGHT_SIDEBAR}>
    <GeneExpression />
  </div>
);

export default connect((state: RootState) => ({
  scatterplotXXaccessor: state.controls.scatterplotXXaccessor,
  scatterplotYYaccessor: state.controls.scatterplotYYaccessor,
}))(RightSidebar);
