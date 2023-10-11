import {
  TypedArray,
  SortableArray,
  NumberArrayConstructor,
  IntTypedArray,
  UnsignedIntTypedArray,
} from "../../common/types/arraytypes";

/**
 * Cross-filter shared types.
 */
export type NonFloatSortableArray =
  | Array<number | string | boolean>
  | IntTypedArray
  | UnsignedIntTypedArray;
export type IndexArray = Uint32Array;

/**
 * Crossfilter selection - see Crossfilter.select()
 */
type SelectableValue = string | number;
export enum SelectionMode {
  All = "all",
  None = "none",
  Exact = "exact",
  Range = "range",
  WithinRect = "within-rect",
  WithinPolygon = "within-polygon",
}

export interface SelectAll {
  mode: SelectionMode.All | "all";
}

export interface SelectNone {
  mode: SelectionMode.None | "none";
}

export interface SelectExact {
  mode: SelectionMode.Exact | "exact";
  values: SelectableValue | SelectableValue[];
}

export interface SelectRange {
  mode: SelectionMode.Range | "range";
  inclusive?: boolean;
  lo: SelectableValue;
  hi: SelectableValue;
}

export interface SelectWithinRect {
  mode: SelectionMode.WithinRect | "within-rect";
  minX: number;
  minY: number;
  maxX: number;
  maxY: number;
}

export interface SelectWithinPolygon {
  mode: SelectionMode.WithinPolygon | "within-polygon";
  polygon: [number, number][];
}

export type CrossfilterSelector =
  | SelectAll
  | SelectNone
  | SelectExact
  | SelectRange
  | SelectWithinRect
  | SelectWithinPolygon;

/**
 * Add dimension to crossfilter - see Crossfilter.addDimension()
 */
export interface ScalarDimensionParameters {
  type: "scalar";
  ValueArrayCtor?: NumberArrayConstructor;
  value: SortableArray;
}

export interface EnumDimensionParameters {
  type: "enum";
  value: SortableArray;
}

export interface SpatialDimensionParameters {
  type: "spatial";
  X: TypedArray;
  Y: TypedArray;
}

export type CrossfilterDimensionParameters =
  | ScalarDimensionParameters
  | EnumDimensionParameters
  | SpatialDimensionParameters;
