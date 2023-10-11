import React from "react";

import { Button, Icon } from "@blueprintjs/core";
import { connect } from "react-redux";
import { Icon as InfoCircle, IconButton } from "czifui";
import Truncate from "../util/truncate";
import HistogramBrush from "../brushableHistogram";
import { RootState } from "../../reducers";

import actions from "../../actions";

import { track } from "../../analytics";
import { EVENTS } from "../../analytics/events";
import { DataframeValue } from "../../util/dataframe";

const MINI_HISTOGRAM_WIDTH = 110;

type State = RootState;

interface Props {
  gene: string;
  // eslint-disable-next-line @typescript-eslint/no-explicit-any -- FIXME
  quickGene: any;
  // eslint-disable-next-line @typescript-eslint/no-explicit-any -- FIXME
  removeGene: any;
  geneId: DataframeValue;
  isGeneExpressionComplete: boolean;
  onGeneExpressionComplete: () => void;
}

// @ts-expect-error ts-migrate(1238) FIXME: Unable to resolve signature of class decorator whe... Remove this comment to see the full error message
@connect((state: RootState, ownProps: Props) => {
  const { gene } = ownProps;

  return {
    isColorAccessor:
      state.colors.colorAccessor === gene &&
      state.colors.colorMode !== "color by categorical metadata",
    isScatterplotXXaccessor: state.controls.scatterplotXXaccessor === gene,
    isScatterplotYYaccessor: state.controls.scatterplotYYaccessor === gene,
    isGeneInfo: state.controls.gene === gene && state.controls.geneIsOpen,
  };
})
class Gene extends React.Component<Props, State> {
  constructor(props: Props) {
    super(props);
    this.state = {
      geneIsExpanded: false,
    };
  }

  onColorChangeClick = (): void => {
    // @ts-expect-error ts-migrate(2339) FIXME: Property 'dispatch' does not exist on type 'Readon... Remove this comment to see the full error message
    const { dispatch, gene } = this.props;
    track(EVENTS.EXPLORER_COLORBY_GENE_BUTTON_CLICKED);
    dispatch(actions.requestSingleGeneExpressionCountsForColoringPOST(gene));
  };

  handleGeneExpandClick = (): void => {
    track(EVENTS.EXPLORER_MAXIMIZE_GENE_BUTTON_CLICKED);

    const { geneIsExpanded } = this.state;
    this.setState({ geneIsExpanded: !geneIsExpanded });
  };

  handleSetGeneAsScatterplotX = (): void => {
    track(EVENTS.EXPLORER_PLOT_X_BUTTON_CLICKED);
    // @ts-expect-error ts-migrate(2339) FIXME: Property 'dispatch' does not exist on type 'Readon... Remove this comment to see the full error message
    const { dispatch, gene } = this.props;
    dispatch({
      type: "set scatterplot x",
      data: gene,
    });
  };

  handleSetGeneAsScatterplotY = (): void => {
    track(EVENTS.EXPLORER_PLOT_Y_BUTTON_CLICKED);
    // @ts-expect-error ts-migrate(2339) FIXME: Property 'dispatch' does not exist on type 'Readon... Remove this comment to see the full error message
    const { dispatch, gene } = this.props;
    dispatch({
      type: "set scatterplot y",
      data: gene,
    });
  };

  handleDeleteGeneFromSet = (): void => {
    // @ts-expect-error ts-migrate(2339) FIXME: Property 'dispatch' does not exist on type 'Readon... Remove this comment to see the full error message
    const { dispatch, gene, geneset } = this.props;
    dispatch(actions.genesetDeleteGenes(geneset, [gene]));
  };

  handleDisplayGeneInfo = async (): Promise<void> => {
    // @ts-expect-error ts-migrate(2339) FIXME: Property 'dispatch' does not exist on type 'Readon... Remove this comment to see the full error message
    const { dispatch, gene, geneId, isGeneInfo } = this.props;
    track(EVENTS.EXPLORER_GENE_INFO_BUTTON_CLICKED, {
      gene,
    });

    if (isGeneInfo) {
      dispatch({
        type: "clear gene info",
      });
      return;
    }

    dispatch({
      type: "load gene info",
      gene,
    });

    const info = await actions.fetchGeneInfo(geneId, gene);
    if (!info) {
      dispatch({
        type: "open gene info",
        gene,
        url: "",
        name: "",
        synonyms: [],
        summary: "",
        infoError: "fetch gene info failed",
      });
      return;
    }
    dispatch({
      type: "open gene info",
      gene,
      url: info.ncbi_url,
      name: info.name,
      synonyms: info.synonyms,
      summary: info.summary,
      infoError: null,
      showWarningBanner: info.show_warning_banner,
    });
  };

  render(): JSX.Element {
    const {
      gene,
      // @ts-expect-error ts-migrate(2339) FIXME: Property 'geneDescription' does not exist on type ... Remove this comment to see the full error message
      geneDescription,
      // @ts-expect-error ts-migrate(2339) FIXME: Property 'isColorAccessor' does not exist on type ... Remove this comment to see the full error message
      isColorAccessor,
      // @ts-expect-error ts-migrate(2339) FIXME: Property 'isScatterplotXXaccessor' does not exist on type ... Remove this comment to see the full error message
      isScatterplotXXaccessor,
      // @ts-expect-error ts-migrate(2339) FIXME: Property 'isScatterplotYYaccessor' does not exist on type ... Remove this comment to see the full error message
      isScatterplotYYaccessor,
      // @ts-expect-error ts-migrate(2339) FIXME: Property 'isGeneInfo' does not exist on type ... Remove this comment to see the full error message
      isGeneInfo,
      quickGene,
      removeGene,
      onGeneExpressionComplete,
      isGeneExpressionComplete,
    } = this.props;
    const { geneIsExpanded } = this.state;
    const geneSymbolWidth = 60 + (geneIsExpanded ? MINI_HISTOGRAM_WIDTH : 0);

    return (
      <div>
        <div
          style={{
            marginLeft: 5,
            marginRight: 0,
            marginTop: 2,
            display: "flex",
            justifyContent: "space-between",
            alignItems: "center",
          }}
        >
          <div
            role="menuitem"
            // @ts-expect-error ts-migrate(2322) FIXME: Type 'string' is not assignable to type 'number | ... Remove this comment to see the full error message
            tabIndex="0"
            data-testclass="gene-expand"
            data-testid={`${gene}:gene-expand`}
            onKeyPress={() => {}}
            style={{
              cursor: "pointer",
              display: "flex",
              justifyContent: "space-between",
              width: "100%",
            }}
          >
            <div>
              <Truncate
                tooltipAddendum={geneDescription && `: ${geneDescription}`}
              >
                <span
                  style={{
                    width: geneSymbolWidth,
                    display: "inline-block",
                  }}
                  data-testid={`${gene}:gene-label`}
                >
                  {gene}
                </span>
              </Truncate>
            </div>
            <div style={{ display: "inline-block", marginLeft: "0" }}>
              <Button
                small
                minimal
                intent={isGeneInfo ? "primary" : "none"}
                data-testid={`get-info-${gene}`}
                active={isGeneInfo}
                onClick={this.handleDisplayGeneInfo}
                disabled={!isGeneExpressionComplete}
              >
                <IconButton
                  disabled={!isGeneExpressionComplete}
                  sdsSize="small"
                >
                  <div style={{ filter: "grayscale(100%)" }}>
                    <InfoCircle
                      sdsIcon="infoCircle"
                      sdsSize="s"
                      sdsType="iconButton"
                    />
                  </div>
                </IconButton>
              </Button>
            </div>
            {!geneIsExpanded ? (
              <div style={{ width: MINI_HISTOGRAM_WIDTH }}>
                <HistogramBrush
                  // @ts-expect-error ts-migrate(2769) FIXME: No overload matches this call.
                  isUserDefined
                  field={gene}
                  mini
                  width={MINI_HISTOGRAM_WIDTH}
                  onGeneExpressionComplete={onGeneExpressionComplete}
                />
              </div>
            ) : null}
          </div>
          <div style={{ flexShrink: 0, marginLeft: 2 }}>
            <Button
              minimal
              small
              data-testid={`delete-from-geneset:${gene}`}
              onClick={() => {
                track(EVENTS.EXPLORER_DELETE_FROM_GENESET_BUTTON_CLICKED);

                if (quickGene) {
                  removeGene(gene)();
                } else {
                  this.handleDeleteGeneFromSet();
                }
              }}
              intent="none"
              style={{ fontWeight: 700, marginRight: 2 }}
              icon={<Icon icon="trash" iconSize={10} />}
            />
            <Button
              minimal
              small
              data-testid={`plot-x-${gene}`}
              onClick={this.handleSetGeneAsScatterplotX}
              active={isScatterplotXXaccessor}
              intent={isScatterplotXXaccessor ? "primary" : "none"}
              style={{ fontWeight: 700, marginRight: 2 }}
            >
              x
            </Button>
            <Button
              minimal
              small
              data-testid={`plot-y-${gene}`}
              onClick={this.handleSetGeneAsScatterplotY}
              active={isScatterplotYYaccessor}
              intent={isScatterplotYYaccessor ? "primary" : "none"}
              style={{ fontWeight: 700, marginRight: 2 }}
            >
              y
            </Button>
            <Button
              minimal
              small
              data-testclass="maximize"
              data-testid={`maximize-${gene}`}
              onClick={this.handleGeneExpandClick}
              active={geneIsExpanded}
              intent="none"
              icon={<Icon icon="maximize" iconSize={10} />}
              style={{ marginRight: 2 }}
            />
            <Button
              minimal
              small
              data-testclass="colorby"
              data-testid={`colorby-${gene}`}
              onClick={this.onColorChangeClick}
              active={isColorAccessor}
              intent={isColorAccessor ? "primary" : "none"}
              icon={<Icon icon="tint" iconSize={12} />}
            />
          </div>
        </div>
        {geneIsExpanded && (
          <HistogramBrush
            // @ts-expect-error ts-migrate(2769) FIXME: No overload matches this call.
            isUserDefined
            field={gene}
            onGeneExpressionComplete={() => {}}
          />
        )}
      </div>
    );
  }
}

export default Gene;
