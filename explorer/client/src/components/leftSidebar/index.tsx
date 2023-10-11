/* Core dependencies */
import React, { CSSProperties } from "react";
import { connect } from "react-redux";

/* App dependencies */
import Categorical from "../categorical";
import * as globals from "../../globals";
import DynamicScatterplot from "../scatterplot/scatterplot";
import GeneInfo from "../geneExpression/geneInfo/geneInfo";
import Continuous from "../continuous/continuous";
import { RootState } from "../../reducers";

/* Styles */
export const STYLE_LEFT_SIDEBAR: CSSProperties = {
  /* x y blur spread color */
  borderRight: `1px solid ${globals.lightGrey}`,
  display: "flex",
  flexDirection: "column",
  height: "100%",
};

interface Props {
  scatterplotXXaccessor: string;
  scatterplotYYaccessor: string;
  geneIsOpen: boolean;
}

const LeftSideBar = (props: Props) => {
  const { scatterplotXXaccessor, scatterplotYYaccessor, geneIsOpen } = props;
  return (
    <div style={STYLE_LEFT_SIDEBAR}>
      <div
        style={{
          height: "100%",
          width: globals.leftSidebarWidth,
          overflowY: "auto",
        }}
      >
        <Categorical />
        <Continuous />
      </div>
      {scatterplotXXaccessor && scatterplotYYaccessor ? (
        <DynamicScatterplot />
      ) : null}
      {geneIsOpen ? (
        <GeneInfo
          geneSummary=""
          geneName=""
          gene=""
          geneUrl=""
          geneSynonyms={[]}
          showWarningBanner
        />
      ) : null}
    </div>
  );
};

export default connect((state: RootState) => ({
  scatterplotXXaccessor: state.controls.scatterplotXXaccessor,
  scatterplotYYaccessor: state.controls.scatterplotYYaccessor,
  geneIsOpen: state.controls.geneIsOpen,
}))(LeftSideBar);
