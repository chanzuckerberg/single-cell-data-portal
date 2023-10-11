/**
 * Gene set state. Geneset UI state is in a different reducer.
 *
 * geneset reducer state is a Map object, where:
 *  key: the geneset name, a string.
 *  val: the geneset defined as an object ("geneset object")
 *
 * A geneset object is:
 * {
 *    genesetName: <string>  # same as the map key
 *    genesetDescription: <string>
 *    genes: Map<<string>, {
 *      geneSymbol: <string>,  # same as the map key
 *      geneDescription: <string>
 *    }>
 * }
 *
 * Geneset and genes Map order is significant, and will be preserved across
 * CRUD operations on either.
 *
 * This reducer does light error checking, but not as much as the backend
 * routes. Do not rely on it to enforce geneset integrity - eg, no duplicate
 * genes in a geneset.
 */
import { AnyAction } from "redux";
import { diffexpPopNamePrefix1, diffexpPopNamePrefix2 } from "../globals";

interface Gene {
  geneSymbol: string;
  geneDescription: string;
}

interface Geneset {
  genesetName: string;
  genesetDescription: string;
  genes: Map<Gene["geneSymbol"], Gene>;
}

export type Genesets = Map<Geneset["genesetName"], Geneset>;

export interface State {
  initialized: boolean;
  lastTid?: number;
  genesets: Genesets;
}

const GeneSets = (
  state: State = {
    initialized: false,
    lastTid: undefined,
    genesets: new Map(),
  },
  action: AnyAction,
): State => {
  switch (action.type) {
    /**
     * Creates a new & empty geneset with the given name and description.
     * {
     *    type: "geneset: create",
     *    genesetName: string, // gene set name
     *    genesetDescription: string, // geneset description
     * }
     *
     */
    case "geneset: create": {
      const { genesetName, genesetDescription } = action;
      if (
        typeof genesetName !== "string" ||
        !genesetName ||
        genesetDescription === undefined
      )
        throw new Error("geneset: create -- name or description unspecified.");
      if (state.genesets.has(genesetName))
        throw new Error("geneset: create -- name already defined.");

      const genesets = new Map([
        [
          genesetName,
          {
            genesetName,
            genesetDescription,
            genes: new Map(),
          },
        ],
        ...state.genesets,
      ]); // clone and add new geneset to beginning

      return {
        ...state,
        genesets,
      };
    }

    /**
     * Deletes the named geneset, if it exists. Throws if it does not.
     * {
     *    type: "geneset: delete",
     *    genesetName: string
     * }
     */
    case "geneset: delete": {
      const { genesetName } = action;
      if (!state.genesets.has(genesetName))
        throw new Error("geneset: delete -- geneset name does not exist.");

      const genesets = new Map(state.genesets); // clone
      genesets.delete(genesetName);
      return {
        ...state,
        genesets,
      };
    }

    /**
     * Update the named geneset with a new name and description.  Preserves the existing
     * order of the geneset, even when the genesetName changes.
     * {
     *    type: "geneset: update",
     *    genesetName: string, current name of geneset to be updated
     *    update: {
     *        genesetName: string, new name
     *        genesetDescription: string, new description
     *    }
     * }
     *
     * For example, if you want to update JUST the description:
     *   dispatch({
     *     action: "geneset: update",
     *     genesetName: "foo",
     *     update: { genesetName: "foo", genesetDescription: "a new description"}
     *   })
     */
    case "geneset: update": {
      const { genesetName, update } = action;

      if (
        typeof genesetName !== "string" ||
        !genesetName ||
        !state.genesets.has(genesetName)
      )
        throw new Error(
          "geneset: update -- geneset name unspecified or does not exist.",
        );

      /* now that we've confirmed the gene set exists, check for duplicates */
      const genesetNameIsDuplicate = state.genesets.has(update.genesetName);
      const descriptionIsDuplicate =
        state.genesets.get(update.genesetName)?.genesetDescription ===
        update.genesetDescription;

      if (genesetNameIsDuplicate && descriptionIsDuplicate)
        throw new Error(
          "geneset: update -- update specified existing name and description.",
        );

      const prevGs = state.genesets.get(genesetName);
      const newGs = {
        ...update,
        genes: prevGs?.genes,
      }; // clone

      // clone the map, preserving current insert order, but mapping name->newName.
      const genesets = new Map();
      for (const [name, gs] of state.genesets) {
        if (name === genesetName) genesets.set(newGs.genesetName, newGs);
        else genesets.set(name, gs);
      }

      return {
        ...state,
        genesets,
      };
    }

    /**
     * Adds genes to the geneset.  They are appended to the END of the geneset, in the
     * order provided. Duplicates or genes already in the geneset, will be ignored.
     * {
     *    type: "geneset: add genes"
     *    genesetName: <string>, // gene set name
     *    genes: Array<{
     *      geneSymbol: <string>,
     *      geneDescription: <string>
     *    }>
     * }
     *
     * Example:
     *   dispatch({
     *     type: "add genes",
     *     genesetName: "foo",
     *     genes: [ { geneSymbol: "FOXP", geneDescription: "test" }]
     *   });
     */
    case "geneset: add genes": {
      const { genesetName, genes } = action;

      const genesets = cloneOneGenesetInGenesets(state.genesets, genesetName);
      const geneset = genesets.get(genesetName) as Geneset;

      // add
      const newGenes = geneset.genes;
      for (const gene of genes) {
        const { geneSymbol } = gene;
        const geneDescription = gene?.geneDescription ?? "";
        // ignore genes already present
        if (!newGenes.has(geneSymbol))
          newGenes.set(geneSymbol, {
            geneSymbol,
            geneDescription,
          });
      }

      return {
        ...state,
        genesets,
      };
    }

    /**
     * Delete genes from the named geneset. Will throw if the genesetName does
     * not exist.  Will ignore geneSymbols that do not exist.
     * {
     *    type: "geneset: delete genes",
     *    genesetName: <string>, // the geneset from which to delete genes
     *    geneSymbols: [<string>, ...], // the gene symbols to delete.
     * }
     *
     * Example:
     *  dispatch({
     *    type: "geneset: delete genes",
     *    genesetName: "a geneset name",
     *    geneSymbols: ["F5"]
     *  })
     */
    case "geneset: delete genes": {
      const { genesetName, geneSymbols } = action;

      const genesets = cloneOneGenesetInGenesets(state.genesets, genesetName);
      const geneset = genesets.get(genesetName) as Geneset;

      // delete
      const { genes } = geneset;
      for (const geneSymbol of geneSymbols) {
        genes.delete(geneSymbol);
      }

      return {
        ...state,
        genesets,
      };
    }

    /**
     * Set/update the description of the gene.  NOTE that this does not allow the name
     * of the gene to change - only "geneset: add" and "geneset: delete" can change
     * the genes in a geneset.  Use this to update a gene description AFTER you add it
     * to the geneset.
     * {
     *    type: "geneset: set gene description",
     *    genesetName: <string>, // the geneset to update
     *    update: {
     *      geneSymbol: <string>, // the gene to update, MUST exist already in the geneset
     *      geneDescription: <string>
     *    }
     * }
     *
     * Example:
     *  dispatch({
     *      type: "geneset: set gene description",
     *      genesetName: "my fav geneset",
     *      update: {
     *        geneSymbol: "F5",
     *        geneDescription: "tada, moar description"
     *      }
     *  })
     */
    case "geneset: set gene description": {
      const { genesetName, update } = action;

      const genesets = cloneOneGenesetInGenesets(state.genesets, genesetName);
      const geneset = genesets.get(genesetName) as Geneset;

      const { geneSymbol, geneDescription } = update;
      const gene = geneset.genes.get(geneSymbol);

      if (!gene) {
        throw new Error("geneset: set gene description -- no such gene");
      }

      geneset.genes.set(geneSymbol, {
        geneSymbol,
        geneDescription,
      });

      return {
        ...state,
        genesets,
      };
    }

    case "request differential expression success": {
      const { data } = action;

      const dateString = new Date().toLocaleString();

      const genesetNames = {
        positive: `${diffexpPopNamePrefix1} (${dateString})`,
        negative: `${diffexpPopNamePrefix2} (${dateString})`,
      };

      const diffExpGeneSets = [];
      for (const polarity of Object.keys(genesetNames)) {
        const genes = new Map(
          // eslint-disable-next-line @typescript-eslint/no-explicit-any --- FIXME: disabled temporarily on migrate to TS.
          data[polarity].map((diffExpGene: any) => [
            diffExpGene[0],
            {
              geneSymbol: diffExpGene[0],
            },
          ]),
        );
        diffExpGeneSets.push([
          // @ts-expect-error ts-migrate(7053) FIXME: Element implicitly has an 'any' type because expre... Remove this comment to see the full error message
          genesetNames[polarity],
          {
            // @ts-expect-error ts-migrate(7053) FIXME: Element implicitly has an 'any' type because expre... Remove this comment to see the full error message
            genesetName: genesetNames[polarity],
            genesetDescription: "",
            genes,
          },
        ]);
      }

      // @ts-expect-error ts-migrate(2769) FIXME: No overload matches this call.
      const genesets = new Map([
        ...diffExpGeneSets,
        ...state.genesets,
      ]) as Genesets; // clone

      return {
        ...state,
        genesets,
      };
    }

    default:
      return state;
  }
};

export default GeneSets;

function cloneOneGenesetInGenesets(genesets: Genesets, genesetName: string) {
  const oldGeneset = genesets.get(genesetName);

  if (!oldGeneset) {
    throw new Error(
      "geneset: set gene description -- geneset name does not exist.",
    );
  }

  const newGenesets = new Map(genesets);
  const newGeneset = {
    ...oldGeneset,
    genes: new Map(oldGeneset.genes),
  };
  newGenesets.set(genesetName, newGeneset);

  return newGenesets;
}
