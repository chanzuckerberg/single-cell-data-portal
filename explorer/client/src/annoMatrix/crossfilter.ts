/*
Row crossfilter proxy for an AnnoMatrix.  This wraps Crossfilter,
providing a number of services, and ensuring that the crossfilter and
AnnoMatrix stay in sync:
	- on-demand index creation as data is loaded
	- transparently mapping between queries and crossfilter index names.
  - for mutation of the matrix by user annotations, maintain synchronization
    between Crossfilter and AnnoMatrix.
*/
import Crossfilter, {
  CrossfilterSelector,
  CrossfilterDimensionParameters,
} from "../util/typedCrossfilter";
import { _getColumnSchema } from "./schema";
import { Field } from "../common/types/schema";
import AnnoMatrix from "./annoMatrix";
import { Dataframe, LabelType } from "../util/dataframe";
import { Query } from "./query";
import { TypedArray } from "../common/types/arraytypes";
import { LabelArray } from "../util/dataframe/types";
import * as globals from "../globals";

function _dimensionNameFromDf(field: Field, df: Dataframe): string {
  const colNames = df.colIndex.labels();
  return _dimensionName(field, colNames);
}

function _dimensionName(
  field: Field,
  colNames: LabelType | LabelArray,
): string {
  if (!Array.isArray(colNames)) return `${field}/${colNames}`;
  return `${field}/${colNames.join(":")}`;
}

export default class AnnoMatrixObsCrossfilter {
  annoMatrix: AnnoMatrix;

  obsCrossfilter: Crossfilter<Dataframe>;

  constructor(
    annoMatrix: AnnoMatrix,
    _obsCrossfilter: Crossfilter<Dataframe> | null = null,
  ) {
    this.annoMatrix = annoMatrix;
    this.obsCrossfilter =
      _obsCrossfilter || new Crossfilter<Dataframe>(annoMatrix._cache.obs);
    this.obsCrossfilter = this.obsCrossfilter.setData(annoMatrix._cache.obs);
  }

  size(): number {
    return this.obsCrossfilter.size();
  }

  /**
   * Drop the crossfilter dimension. Do not change the annoMatrix. Useful when we
   * want to stop tracking the selection state, but aren't sure we want to blow the
   * annomatrix cache.
   */
  dropDimension(field: Field, query: Query): AnnoMatrixObsCrossfilter {
    const { annoMatrix } = this;
    let { obsCrossfilter } = this;
    const keys = annoMatrix
      .getCacheKeys(field, query)
      // @ts-expect-error ts-migrate --- suppressing TS defect (https://github.com/microsoft/TypeScript/issues/44373).
      // Compiler is complaining that expression is not callable on array union types. Remove suppression once fixed.
      .filter((k?: string | number) => k !== undefined);
    const dimName = _dimensionName(field, keys as string[]);
    if (obsCrossfilter.hasDimension(dimName)) {
      obsCrossfilter = obsCrossfilter.delDimension(dimName);
    }
    return new AnnoMatrixObsCrossfilter(annoMatrix, obsCrossfilter);
  }

  /**
  Selection state - API is identical to ImmutableTypedCrossfilter, as these
  are just wrappers to lazy create indices.
  **/

  async select(
    field: Field,
    query: Query,
    spec: CrossfilterSelector,
  ): Promise<AnnoMatrixObsCrossfilter> {
    const { annoMatrix } = this;
    let { obsCrossfilter } = this;

    if (!annoMatrix?._cache?.[field]) {
      throw new Error("Unknown field name");
    }
    if (field === "var") {
      throw new Error("unable to obsSelect upon the var dimension");
    }

    // set the correct number of bins for lossy compression of float arrays based on field
    let nBins;
    switch (field) {
      case "emb": {
        nBins = globals.numBinsEmb;
        break;
      }
      case "obs":
      case "X": {
        nBins = globals.numBinsObsX;
        break;
      }
      default: {
        nBins = null;
        break;
      }
    }

    // grab the data, so we can grab the index.
    const df = await annoMatrix.fetch(field, query, nBins);
    if (!df) {
      throw new Error("Dataframe cannot be `undefined`");
    }
    const dimName = _dimensionNameFromDf(field, df);
    if (!obsCrossfilter.hasDimension(dimName)) {
      // lazy index generation - add dimension when first used
      obsCrossfilter = this._addObsCrossfilterDimension(
        annoMatrix,
        obsCrossfilter,
        field,
        df,
      );
    }

    // select
    obsCrossfilter = obsCrossfilter.select(dimName, spec);
    return new AnnoMatrixObsCrossfilter(annoMatrix, obsCrossfilter);
  }

  selectAll(): AnnoMatrixObsCrossfilter {
    /*
		Select all on any dimension in this field.
		*/
    const { annoMatrix } = this;
    const currentDims = this.obsCrossfilter.dimensionNames();
    const obsCrossfilter = currentDims.reduce(
      (xfltr, dim) => xfltr.select(dim, { mode: "all" }),
      this.obsCrossfilter,
    );
    return new AnnoMatrixObsCrossfilter(annoMatrix, obsCrossfilter);
  }

  countSelected(): number {
    /* if no data yet indexed in the crossfilter, just say everything is selected */
    if (this.obsCrossfilter.size() === 0) return this.annoMatrix.nObs;
    return this.obsCrossfilter.countSelected();
  }

  allSelectedMask(): Uint8Array {
    /* if no data yet indexed in the crossfilter, just say everything is selected */
    if (
      this.obsCrossfilter.size() === 0 ||
      this.obsCrossfilter.dimensionNames().length === 0
    ) {
      /* fake the mask */
      return new Uint8Array(this.annoMatrix.nObs).fill(1);
    }
    return this.obsCrossfilter.allSelectedMask();
  }

  allSelectedLabels(): LabelArray {
    /* if no data yet indexed in the crossfilter, just say everything is selected */
    if (
      this.obsCrossfilter.size() === 0 ||
      this.obsCrossfilter.dimensionNames().length === 0
    ) {
      return this.annoMatrix.rowIndex.labels();
    }

    const mask = this.obsCrossfilter.allSelectedMask();
    const index = this.annoMatrix.rowIndex.isubsetMask(mask);
    return index.labels();
  }

  fillByIsSelected<A extends TypedArray>(
    array: A,
    selectedValue: A[0],
    deselectedValue: A[0],
  ): A {
    /* if no data yet indexed in the crossfilter, just say everything is selected */
    if (
      this.obsCrossfilter.size() === 0 ||
      this.obsCrossfilter.dimensionNames().length === 0
    ) {
      array.fill(selectedValue);
      return array;
    }
    return this.obsCrossfilter.fillByIsSelected(
      array,
      selectedValue,
      deselectedValue,
    );
  }

  /**
   ** Private below
   **/

  _addObsCrossfilterDimension(
    annoMatrix: AnnoMatrix,
    obsCrossfilter: Crossfilter<Dataframe>,
    field: Field,
    df: Dataframe,
  ): Crossfilter<Dataframe> {
    if (field === "var") return obsCrossfilter;
    const dimName = _dimensionNameFromDf(field, df);
    const dimParams = this._getObsDimensionParams(field, df);
    if (!dimName || !dimParams) return obsCrossfilter;
    obsCrossfilter = obsCrossfilter.setData(annoMatrix._cache.obs);
    obsCrossfilter = obsCrossfilter.addDimension(dimName, dimParams);
    return obsCrossfilter;
  }

  _getColumnBaseType(field: Field, col: LabelType): string {
    /* Look up the primitive type for this field/col */
    const colSchema = _getColumnSchema(this.annoMatrix.schema, field, col);
    return colSchema.type;
  }

  _getObsDimensionParams(
    field: Field,
    df: Dataframe,
  ): CrossfilterDimensionParameters | undefined {
    /* return the crossfilter dimension type and params for this field/dataframe */

    if (field === Field.emb) {
      /* assumed to be 2D and float, as that is all the schema supports */
      return {
        type: "spatial",
        X: df.icol(0).asArray() as TypedArray,
        Y: df.icol(1).asArray() as TypedArray,
      };
    }

    /* assumed to be 1D */
    const col = df.icol(0);
    const colName: LabelType = df.colIndex.getLabel(0) as LabelType;
    const type = this._getColumnBaseType(field, colName);
    if (type === "string" || type === "categorical" || type === "boolean") {
      return { type: "enum", value: col.asArray() };
    }
    if (type === "int32") {
      return {
        type: "scalar",
        value: col.asArray(),
        ValueArrayCtor: Int32Array,
      };
    }
    if (type === "float32") {
      return {
        type: "scalar",
        value: col.asArray(),
        ValueArrayCtor: Float32Array,
      };
    }
    // Currently not supporting boolean and categorical types.
    console.error(
      `Warning - unknown metadata schema (${type}) for field ${field} ${colName}.`,
    );
    // skip it - we don't know what to do with this type

    return undefined;
  }
}
