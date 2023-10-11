import React from "react";
import { FaChevronRight, FaChevronDown } from "react-icons/fa";
import Gene from "./gene";
import Truncate from "../util/truncate";
import * as globals from "../../globals";
import GenesetMenus from "./menus/genesetMenus";
import EditGenesetNameDialogue from "./menus/editGenesetNameDialogue";
import HistogramBrush from "../brushableHistogram";

import { diffexpPopNamePrefix1, diffexpPopNamePrefix2 } from "../../globals";
import { track } from "../../analytics";
import { EVENTS } from "../../analytics/events";

// eslint-disable-next-line @typescript-eslint/no-explicit-any --- FIXME: disabled temporarily on migrate to TS.
type State = any;

// eslint-disable-next-line @typescript-eslint/ban-types --- FIXME: disabled temporarily on migrate to TS.
class GeneSet extends React.Component<{}, State> {
  isGeneExpressionLoadComplete = false;

  // eslint-disable-next-line @typescript-eslint/ban-types --- FIXME: disabled temporarily on migrate to TS.
  constructor(props: {}) {
    super(props);
    this.state = {
      isOpen: false,
      genesetLoadCount: 0,
    };
  }

  // eslint-disable-next-line @typescript-eslint/explicit-module-boundary-types --- FIXME: disabled temporarily on migrate to TS.
  onGenesetMenuClick = () => {
    const { isOpen } = this.state;

    if (!isOpen) {
      track(EVENTS.EXPLORER_GENESET_EXPAND_BUTTON_CLICKED);
    }

    this.setState({ isOpen: !isOpen });
  };

  onGeneExpressionComplete = (): void => {
    const { genesetLoadCount } = this.state;
    this.setState({ genesetLoadCount: genesetLoadCount + 1 });
  };

  // eslint-disable-next-line @typescript-eslint/explicit-module-boundary-types, @typescript-eslint/no-explicit-any -- - FIXME: disabled temporarily on migrate to TS.
  componentDidUpdate = () => {
    // @ts-expect-error ts-migrate(2339) FIXME: Property 'setGenes' does not exist on type 'Readonl... Remove this comment to see the full error message
    const { setGenes } = this.props;
    const { genesetLoadCount } = this.state;
    if (genesetLoadCount + 1 >= setGenes.size) {
      this.isGeneExpressionLoadComplete = true;
    }
  };

  // eslint-disable-next-line @typescript-eslint/explicit-module-boundary-types --- FIXME: disabled temporarily on migrate to TS.
  renderGenes() {
    // @ts-expect-error ts-migrate(2339) FIXME: Property 'setName' does not exist on type 'Readonl... Remove this comment to see the full error message
    const { setName, setGenes, geneIds, geneNames } = this.props;
    const setGenesNames = [...setGenes.keys()];
    return (
      <div data-testclass="gene-set-genes">
        {setGenesNames.map((gene) => {
          const { geneDescription } = setGenes.get(gene);
          const geneId = geneIds[geneNames.indexOf(gene)];
          return (
            <Gene
              key={gene}
              gene={gene}
              // @ts-expect-error ts-migrate(2322) FIXME: Type '{ key: any; gene: any; geneDescription: any;... Remove this comment to see the full error message
              geneDescription={geneDescription}
              geneset={setName}
              geneId={geneId}
              onGeneExpressionComplete={this.onGeneExpressionComplete}
              isGeneExpressionComplete={this.isGeneExpressionLoadComplete}
            />
          );
        })}
      </div>
    );
  }

  // eslint-disable-next-line @typescript-eslint/explicit-module-boundary-types --- FIXME: disabled temporarily on migrate to TS.
  render() {
    // @ts-expect-error ts-migrate(2339) FIXME: Property 'setName' does not exist on type 'Readonl... Remove this comment to see the full error message
    const { setName, genesetDescription, setGenes } = this.props;
    const { isOpen } = this.state;
    const genesetNameLengthVisible = 150; /* this magic number determines how much of a long geneset name we see */
    const genesetIsEmpty = setGenes.size === 0;
    let testClass = "geneset-expand";

    if (setName.includes(diffexpPopNamePrefix1))
      testClass = "pop-1-geneset-expand";
    else if (setName.includes(diffexpPopNamePrefix2))
      testClass = "pop-2-geneset-expand";

    return (
      <div data-testclass="geneset" style={{ marginBottom: 3 }}>
        <div
          style={{
            display: "flex",
            justifyContent: "space-between",
            alignItems: "baseline",
          }}
        >
          <span
            role="menuitem"
            // @ts-expect-error ts-migrate(2322) FIXME: Type 'string' is not assignable to type 'number | ... Remove this comment to see the full error message
            tabIndex="0"
            data-testclass={testClass}
            data-testid={`${setName}:geneset-expand`}
            onKeyPress={
              /* TODO(colinmegill): #2101: click handler on span */ () => {}
            }
            style={{
              cursor: "pointer",
              userSelect: "none",
            }}
            onClick={this.onGenesetMenuClick}
          >
            <Truncate
              isGenesetDescription
              tooltipAddendum={
                genesetDescription ? `: ${genesetDescription}` : ""
              }
            >
              <span
                style={{
                  maxWidth: globals.leftSidebarWidth - genesetNameLengthVisible,
                }}
                data-testid={`${setName}:geneset-name`}
              >
                {setName}
              </span>
            </Truncate>
            {isOpen ? (
              <FaChevronDown
                data-testclass="geneset-expand-is-expanded"
                style={{ fontSize: 10, marginLeft: 5 }}
              />
            ) : (
              <FaChevronRight
                data-testclass="geneset-expand-is-not-expanded"
                style={{ fontSize: 10, marginLeft: 5 }}
              />
            )}
          </span>
          <div>
            {/* @ts-expect-error ts-migrate(2322) FIXME: Type '{ isOpen: any; genesetsEditable: true; genes... Remove this comment to see the full error message */}
            <GenesetMenus isOpen={isOpen} genesetsEditable geneset={setName} />
          </div>
        </div>

        <div style={{ marginLeft: 15, marginTop: 5, marginRight: 0 }}>
          {isOpen && genesetIsEmpty && (
            <p style={{ fontStyle: "italic", color: "lightgrey" }}>
              No genes to display
            </p>
          )}
        </div>
        {isOpen && !genesetIsEmpty && (
          <HistogramBrush
            // @ts-expect-error ts-migrate(2769) FIXME: No overload matches this call.
            isGeneSetSummary
            field={setName}
            setGenes={setGenes}
            onGeneExpressionComplete={() => {}}
          />
        )}
        {isOpen && !genesetIsEmpty && this.renderGenes()}
        <EditGenesetNameDialogue
          // @ts-expect-error ts-migrate(2322) FIXME: Type '{ parentGeneset: any; parentGenesetDescripti... Remove this comment to see the full error message
          parentGeneset={setName}
          parentGenesetDescription={genesetDescription}
        />
      </div>
    );
  }
}

export default GeneSet;
