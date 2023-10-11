/* eslint-disable max-classes-per-file -- Classes are interrelated*/

/*
Views on the annomatrix.  all API here is defined in viewCreators.js and annoMatrix.js.
*/
import clip from "../util/clip";
import AnnoMatrix from "./annoMatrix";
import { _whereCacheCreate, WhereCache } from "./whereCache";
import { _isContinuousType, _getColumnSchema } from "./schema";
import { Dataframe, DataframeValueArray, LabelType } from "../util/dataframe";
import { Query } from "./query";
import { ArraySchema, Field } from "../common/types/schema";
import { LabelIndexBase } from "../util/dataframe/labelIndex";
import { isTypedArray } from "../common/types/arraytypes";

type MapFn = (
  field: Field,
  colLabel: LabelType,
  colSchema: ArraySchema,
  colData: DataframeValueArray,
  df: Dataframe,
) => DataframeValueArray;

abstract class AnnoMatrixView extends AnnoMatrix {
  constructor(viewOf: AnnoMatrix, rowIndex: LabelIndexBase | null = null) {
    const nObs = rowIndex ? rowIndex.size() : viewOf.nObs;
    super(viewOf.schema, nObs, viewOf.nVar, rowIndex || viewOf.rowIndex);
    this.viewOf = viewOf;
    this.isView = true;
  }
}

class AnnoMatrixMapView extends AnnoMatrixView {
  mapFn: MapFn;

  /*
    A view which knows how to transform its data.
    */
  constructor(viewOf: AnnoMatrix, mapFn: MapFn) {
    super(viewOf);
    this.mapFn = mapFn;
  }

  async _doLoad(
    field: Field,
    query: Query,
    nBins: number | null = null,
  ): Promise<[WhereCache | null, Dataframe]> {
    const df = await this.viewOf._fetch(field, query, nBins);
    const dfMapped = df.mapColumns(
      (colData: DataframeValueArray, _colIdx: number, colLabel: LabelType) => {
        const colSchema = _getColumnSchema(this.schema, field, colLabel);
        return this.mapFn(field, colLabel, colSchema, colData, df);
      },
    );
    const whereCacheUpdate = _whereCacheCreate(
      field,
      query,
      dfMapped.colIndex.labels(),
    );
    return [whereCacheUpdate, dfMapped];
  }
}

export class AnnoMatrixClipView extends AnnoMatrixMapView {
  clipRange: [number, number];

  isClipped: boolean;

  /*
    A view which is a clipped transformation of its parent
    */
  constructor(viewOf: AnnoMatrix, qmin: number, qmax: number) {
    super(
      viewOf,
      (
        field: Field,
        colLabel: LabelType,
        colSchema: ArraySchema,
        colData: DataframeValueArray,
        df: Dataframe,
      ) => _clipAnnoMatrix(field, colLabel, colSchema, colData, df, qmin, qmax),
    );
    this.isClipped = true;
    this.clipRange = [qmin, qmax];
    Object.seal(this);
  }
}

export class AnnoMatrixRowSubsetView extends AnnoMatrixView {
  /*
    A view which is a subset of total rows.
    */
  constructor(viewOf: AnnoMatrix, rowIndex: LabelIndexBase) {
    super(viewOf, rowIndex);
    Object.seal(this);
  }

  async _doLoad(
    field: Field,
    query: Query,
    nBins: number | null = null,
  ): Promise<[WhereCache | null, Dataframe]> {
    const df = await this.viewOf._fetch(field, query, nBins);

    // don't try to row-subset the var dimension.
    if (field === Field.var) {
      return [null, df];
    }
    const dfSubset = df.subset(null, null, this.rowIndex);
    const whereCacheUpdate = _whereCacheCreate(
      field,
      query,
      dfSubset.colIndex.labels(),
    );
    return [whereCacheUpdate, dfSubset];
  }
}

/*
Utility functions below
*/

function _clipAnnoMatrix(
  field: Field,
  colLabel: LabelType,
  colSchema: ArraySchema,
  colData: DataframeValueArray,
  df: Dataframe,
  qmin: number,
  qmax: number,
): DataframeValueArray {
  /* only clip obs and var scalar columns */
  if (field !== Field.obs && field !== Field.X) return colData;
  if (!_isContinuousType(colSchema) || !isTypedArray(colData)) return colData;
  if (qmin < 0) qmin = 0;
  if (qmax > 1) qmax = 1;
  if (qmin === 0 && qmax === 1) return colData;

  const quantiles = df.col(colLabel).summarizeContinuous().percentiles;
  const lower = quantiles[100 * qmin];
  const upper = quantiles[100 * qmax];
  const clippedData = clip(colData.slice(), lower, upper, Number.NaN);
  return clippedData;
}

/* eslint-enable max-classes-per-file -- enable*/
