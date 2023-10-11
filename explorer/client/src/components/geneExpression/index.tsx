import React from "react";
import { connect } from "react-redux";
import { Button, H4, Icon } from "@blueprintjs/core";
import { IconNames } from "@blueprintjs/icons";

import GeneSet from "./geneSet";
import QuickGene from "./quickGene";
import CreateGenesetDialogue from "./menus/createGenesetDialogue";
import { track } from "../../analytics";
import { EVENTS } from "../../analytics/events";
import { Dataframe, DataframeValue } from "../../util/dataframe";

// eslint-disable-next-line @typescript-eslint/no-explicit-any --- FIXME: disabled temporarily on migrate to TS.
type State = any;

// @ts-expect-error ts-migrate(1238) FIXME: Unable to resolve signature of class decorator whe... Remove this comment to see the full error message
@connect((state: RootState) => ({
  genesets: state.genesets.genesets,
  annoMatrix: state.annoMatrix,
}))
// eslint-disable-next-line @typescript-eslint/ban-types --- FIXME: disabled temporarily on migrate to TS.
class GeneExpression extends React.Component<{}, State> {
  // eslint-disable-next-line @typescript-eslint/ban-types --- FIXME: disabled temporarily on migrate to TS.
  constructor(props: {}) {
    super(props);
    this.state = {
      geneSetsExpanded: true,
      geneIds: null,
      geneNames: null,
    };
  }

  async componentDidMount(): Promise<void> {
    // @ts-expect-error ts-migrate(2339) FIXME: Property 'annoMatrix' does not exist on type 'Readon... Remove this comment to see the full error message
    const { annoMatrix } = this.props;
    const { schema } = annoMatrix;
    const varIndex = schema.annotations.var.index;

    let dfIds;
    const df: Dataframe = await annoMatrix.fetch("var", varIndex);
    this.setState({
      geneNames: df.col(varIndex).asArray() as DataframeValue[],
    });

    const geneIdCol = "feature_id";
    // if feature id column is available in var
    if (annoMatrix.getMatrixColumns("var").includes(geneIdCol)) {
      dfIds = await annoMatrix.fetch("var", geneIdCol);
      this.setState({
        geneIds: dfIds.col(geneIdCol).asArray() as DataframeValue[],
      });
    }
  }

  // eslint-disable-next-line @typescript-eslint/explicit-module-boundary-types --- FIXME: disabled temporarily on migrate to TS.
  renderGeneSets = () => {
    const sets = [];
    // @ts-expect-error ts-migrate(2339) FIXME: Property 'genesets' does not exist on type 'Readon... Remove this comment to see the full error message
    const { genesets } = this.props;
    const { geneIds, geneNames } = this.state;

    for (const [name, geneset] of genesets) {
      const genesetIds = [];
      const genesetNames = [];

      // find ensembl IDs for each gene in the geneset
      for (const gene of geneset.genes) {
        let geneId;
        if (geneIds) {
          geneId = geneIds[geneNames.indexOf(gene[0])];
        } else {
          geneId = "";
        }
        genesetIds.push(geneId);
        genesetNames.push(gene[0]);
      }

      sets.push(
        <GeneSet
          key={name}
          // @ts-expect-error ts-migrate(2322) FIXME: Type '{ key: any; setGenes: any; setName: any; gen... Remove this comment to see the full error message
          setGenes={geneset.genes}
          setName={name}
          genesetDescription={geneset.genesetDescription}
          geneIds={genesetIds}
          geneNames={genesetNames}
        />,
      );
    }
    return sets;
  };

  // eslint-disable-next-line @typescript-eslint/explicit-module-boundary-types --- FIXME: disabled temporarily on migrate to TS.
  handleActivateCreateGenesetMode = () => {
    track(EVENTS.EXPLORER_OPEN_CREATE_GENESET_DIALOG_BUTTON_CLICKED);

    // @ts-expect-error ts-migrate(2339) FIXME: Property 'dispatch' does not exist on type 'Readon... Remove this comment to see the full error message
    const { dispatch } = this.props;
    const { geneSetsExpanded } = this.state;
    dispatch({ type: "geneset: activate add new geneset mode" });
    if (!geneSetsExpanded) {
      // eslint-disable-next-line @typescript-eslint/no-explicit-any --- FIXME: disabled temporarily on migrate to TS.
      this.setState((state: any) => ({ ...state, geneSetsExpanded: true }));
    }
  };

  // eslint-disable-next-line @typescript-eslint/explicit-module-boundary-types --- FIXME: disabled temporarily on migrate to TS.
  handleExpandGeneSets = () => {
    // eslint-disable-next-line @typescript-eslint/no-explicit-any --- FIXME: disabled temporarily on migrate to TS.
    this.setState((state: any) => {
      if (!state.geneSetsExpanded) {
        track(EVENTS.EXPLORER_GENESET_HEADING_EXPAND_BUTTON_CLICKED);
      }

      return {
        ...state,
        geneSetsExpanded: !state.geneSetsExpanded,
      };
    });
  };

  // eslint-disable-next-line @typescript-eslint/explicit-module-boundary-types --- FIXME: disabled temporarily on migrate to TS.
  render() {
    const { geneSetsExpanded } = this.state;
    return (
      <div>
        <QuickGene />
        <div>
          <div
            style={{
              display: "flex",
              flexDirection: "row",
              justifyContent: "space-between",
            }}
          >
            <H4
              role="menuitem"
              // @ts-expect-error ts-migrate(2322) FIXME: Type 'string' is not assignable to type 'number | ... Remove this comment to see the full error message
              tabIndex="0"
              data-testclass="geneset-heading-expand"
              onKeyPress={this.handleExpandGeneSets}
              style={{
                cursor: "pointer",
              }}
              onClick={this.handleExpandGeneSets}
            >
              Gene Sets{" "}
              {geneSetsExpanded ? (
                <Icon icon={IconNames.CHEVRON_DOWN} />
              ) : (
                <Icon icon={IconNames.CHEVRON_RIGHT} />
              )}
            </H4>

            <div style={{ marginBottom: 10, position: "relative", top: -2 }}>
              <Button
                data-testid="open-create-geneset-dialog"
                onClick={this.handleActivateCreateGenesetMode}
                intent="primary"
              >
                Create new
              </Button>
            </div>
          </div>
          <CreateGenesetDialogue />
        </div>

        {geneSetsExpanded && <div>{this.renderGeneSets()}</div>}
      </div>
    );
  }
}

export default GeneExpression;
