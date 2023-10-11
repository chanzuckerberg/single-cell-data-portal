import { H4, Icon, MenuItem } from "@blueprintjs/core";
import { IconNames } from "@blueprintjs/icons";
import React, { useState, useEffect, useRef, useMemo } from "react";
import fuzzysort from "fuzzysort";
import { Suggest } from "@blueprintjs/select";
import { useSelector, useDispatch } from "react-redux";

import Gene from "./gene";

import { postUserErrorToast } from "../framework/toasters";
import actions from "../../actions";
import { Dataframe, DataframeValue } from "../../util/dataframe";
import { track } from "../../analytics";
import { EVENTS } from "../../analytics/events";

// eslint-disable-next-line @typescript-eslint/no-explicit-any --- FIXME: disabled temporarily on migrate to TS.
const usePrevious = (value: any) => {
  const ref = useRef();
  useEffect(() => {
    ref.current = value;
  });
  return ref.current;
};

// eslint-disable-next-line @typescript-eslint/explicit-module-boundary-types --- FIXME: disabled temporarily on migrate to TS.
function QuickGene() {
  const dispatch = useDispatch();

  const [isExpanded, setIsExpanded] = useState(true);
  const [geneNames, setGeneNames] = useState([] as DataframeValue[]);
  const [geneIds, setGeneIds] = useState([] as DataframeValue[]);
  const [, setStatus] = useState("pending");

  const { annoMatrix, userDefinedGenes, userDefinedGenesLoading } = useSelector(
    (state) => ({
      // eslint-disable-next-line @typescript-eslint/no-explicit-any --- FIXME: disabled temporarily on migrate to TS.
      annoMatrix: (state as any).annoMatrix,
      // eslint-disable-next-line @typescript-eslint/no-explicit-any --- FIXME: disabled temporarily on migrate to TS.
      userDefinedGenes: (state as any).quickGenes.userDefinedGenes,
      // eslint-disable-next-line @typescript-eslint/no-explicit-any --- FIXME: disabled temporarily on migrate to TS.
      userDefinedGenesLoading: (state as any).quickGenes
        .userDefinedGenesLoading,
    }),
  );

  const prevProps = usePrevious({ annoMatrix });

  useEffect(() => {
    if (!annoMatrix) return;

    fetch();

    async function fetch() {
      // eslint-disable-next-line @typescript-eslint/no-explicit-any --- FIXME: disabled temporarily on commit
      if (annoMatrix !== (prevProps as any)?.annoMatrix) {
        const { schema } = annoMatrix;
        const varIndex = schema.annotations.var.index;

        setStatus("pending");
        try {
          const df: Dataframe = await annoMatrix.fetch("var", varIndex);
          let dfIds: Dataframe;
          const geneIdCol = "feature_id";

          // if feature id column is available in var
          if (annoMatrix.getMatrixColumns("var").includes(geneIdCol)) {
            dfIds = await annoMatrix.fetch("var", geneIdCol);
            setGeneIds(dfIds.col("feature_id").asArray() as DataframeValue[]);
          }

          setStatus("success");
          setGeneNames(df.col(varIndex).asArray() as DataframeValue[]);
        } catch (error) {
          setStatus("error");
          throw error;
        }
      }
    }
  }, [annoMatrix, prevProps]);

  const handleExpand = () => setIsExpanded(!isExpanded);

  const renderGene = (
    // eslint-disable-next-line @typescript-eslint/no-explicit-any --- FIXME: disabled temporarily on migrate to TS.
    fuzzySortResult: any,
    // eslint-disable-next-line @typescript-eslint/no-explicit-any --- FIXME: disabled temporarily on migrate to TS.
    { handleClick, modifiers }: any,
  ) => {
    if (!modifiers.matchesPredicate) {
      return null;
    }
    /* the fuzzysort wraps the object with other properties, like a score */
    const geneName = fuzzySortResult.target;

    return (
      <MenuItem
        active={modifiers.active}
        disabled={modifiers.disabled}
        data-testid={`suggest-menu-item-${geneName}`}
        key={geneName}
        // eslint-disable-next-line @typescript-eslint/no-explicit-any --- FIXME: disabled temporarily on migrate to TS.
        onClick={(g: any /* this fires when user clicks a menu item */) => {
          track(EVENTS.EXPLORER_SUGGEST_MENU_ITEM_CLICKED);
          handleClick(g);
        }}
        text={geneName}
      />
    );
  };

  // eslint-disable-next-line @typescript-eslint/no-explicit-any --- FIXME: disabled temporarily on migrate to TS.
  const handleClick = (g: any) => {
    if (!g) return;
    const gene = g.target;
    if (userDefinedGenes.indexOf(gene) !== -1) {
      postUserErrorToast("That gene already exists");
    } else if (geneNames.indexOf(gene) === undefined) {
      postUserErrorToast("That doesn't appear to be a valid gene name.");
    } else {
      dispatch({ type: "single user defined gene start" });
      dispatch(actions.requestUserDefinedGene(gene));
      dispatch({ type: "single user defined gene complete" });
    }
  };

  // eslint-disable-next-line @typescript-eslint/no-explicit-any --- FIXME: disabled temporarily on migrate to TS.
  const filterGenes = (query: any, genes: any) =>
    /* fires on load, once, and then for each character typed into the input */
    fuzzysort.go(query, genes, {
      limit: 5,
      threshold: -10000, // don't return bad results
    });

  const QuickGenes = useMemo((): JSX.Element => {
    // eslint-disable-next-line @typescript-eslint/no-explicit-any --- FIXME: disabled temporarily on commit
    const removeGene = (gene: any) => () => {
      dispatch({ type: "clear user defined gene", gene });
    };

    // eslint-disable-next-line @typescript-eslint/no-explicit-any --- FIXME: disabled temporarily on migrate to TS.
    return userDefinedGenes.map((gene: any) => {
      let geneId = geneIds[geneNames.indexOf(gene)];
      if (!geneId) {
        geneId = "";
      }

      return (
        <>
          <Gene
            key={`quick=${gene}`}
            gene={gene}
            removeGene={removeGene}
            quickGene
            geneId={geneId}
            isGeneExpressionComplete
            onGeneExpressionComplete={() => {}}
          />
        </>
      );
    });
  }, [userDefinedGenes, geneNames, geneIds, dispatch]);

  return (
    <div style={{ width: "100%", marginBottom: "16px" }}>
      <H4
        role="menuitem"
        // @ts-expect-error ts-migrate(2322) FIXME: Type 'string' is not assignable to type 'number | ... Remove this comment to see the full error message
        tabIndex="0"
        data-testclass="quickgene-heading-expand"
        onKeyPress={handleExpand}
        style={{
          cursor: "pointer",
        }}
        onClick={handleExpand}
      >
        Genes{" "}
        {isExpanded ? (
          <Icon icon={IconNames.CHEVRON_DOWN} />
        ) : (
          <Icon icon={IconNames.CHEVRON_RIGHT} />
        )}
      </H4>
      {isExpanded && (
        <>
          <div style={{ marginBottom: "8px" }}>
            <Suggest
              resetOnSelect
              closeOnSelect
              resetOnClose
              itemDisabled={userDefinedGenesLoading ? () => true : () => false}
              noResults={<MenuItem disabled text="No matching genes." />}
              onItemSelect={(g) => {
                /* this happens on 'enter' */
                handleClick(g);
              }}
              initialContent={<MenuItem disabled text="Enter a gene…" />}
              inputProps={{
                // @ts-expect-error ts-migrate(2322) FIXME: Type '{ "data-testid": string; placeholder: string... Remove this comment to see the full error message
                "data-testid": "gene-search",
                placeholder: "Quick Gene Search",
                leftIcon: IconNames.SEARCH,
                fill: true,
              }}
              inputValueRenderer={() => ""}
              // @ts-expect-error ts-migrate(2322) FIXME: Type '(query: any, genes: any) => Fuzzysort.Result... Remove this comment to see the full error message
              itemListPredicate={filterGenes}
              itemRenderer={renderGene}
              items={geneNames || ["No genes"]}
              popoverProps={{ minimal: true }}
              fill
            />
          </div>
          {QuickGenes}
        </>
      )}
    </div>
  );
}

export default QuickGene;
