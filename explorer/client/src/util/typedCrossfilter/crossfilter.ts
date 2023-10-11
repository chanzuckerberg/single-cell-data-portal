// eslint-disable-next-line max-classes-per-file -- classes are interrelated
import {
  AnyArray,
  isAnyArray,
  TypedArray,
  SortableArray,
} from "../../common/types/arraytypes";
import {
  IndexArray,
  SelectionMode,
  CrossfilterSelector,
  SelectExact,
  SelectRange,
  SelectWithinRect,
  SelectWithinPolygon,
  ScalarDimensionParameters,
  EnumDimensionParameters,
  SpatialDimensionParameters,
  CrossfilterDimensionParameters,
} from "./types";
import { Dataframe } from "../dataframe";
import PositiveIntervals, {
  Interval,
  IntervalArray,
} from "./positiveIntervals";
import BitArray from "./bitArray";
import {
  sortArray,
  lowerBound,
  binarySearch,
  lowerBoundIndirect,
  upperBoundIndirect,
} from "./sort";
import { makeSortIndex, NotImplementedError } from "./util";

/**
 * The result of a dimension selection
 * @internal
 */
type DimensionSelectionResult = {
  ranges: IntervalArray;
  index: IndexArray | null;
};

/**
 * Table of dimensions added to the crossfilter
 * @internal
 */
interface Dimension {
  id: number | undefined; // ID is often undefined when selection cache is clear
  dim: _ImmutableBaseDimension;
  name: string;
  selection: DimensionSelectionResult;
}
interface Dimensions {
  [name: string]: Dimension;
}

/**
 * Selection performance cache
 * @internal
 */
interface SelectionCache {
  bitArray: BitArray;
  allSelectedMask?: Uint8Array;
  countSelected?: number;
}
type PartialSelectionCache = Partial<SelectionCache>;

/**
 * The type of our "data" (aka records)
 */
export type CrossfilterData = AnyArray | Dataframe;
export function isCrossfilterData(tbd: unknown): tbd is CrossfilterData {
  return isAnyArray(tbd) || tbd instanceof Dataframe;
}

export default class ImmutableTypedCrossfilter<
  RecordType extends CrossfilterData = CrossfilterData,
> {
  public data: RecordType;

  private dimensions: Dimensions;

  private selectionCache: PartialSelectionCache;

  constructor(
    data: RecordType,
    dimensions: Dimensions = {},
    selectionCache: PartialSelectionCache = {},
  ) {
    /*
    Typically, parameter 'data' is one of:
      - Array of objects/records
      - Dataframe (util/dataframe)
    Other parameters are only used internally.

    Object field description:
    - data: reference to the array of records in the crossfilter
    - selectionCache: object which may contains a bit array indicating
      the flatted selection state for all dimensions, plus other state
      summarizing the selection.  This is lazily created and is effectively
      a performance cache.  Methods which return a new crossfilter,
      such as select(), addDimension() and delDimension(), will pass
      the cache forward to the new object, as the typical "immutable API"
      usage pattern is to retain the new crossfilter and discard the old.
      If the cache object is empty, it will be rebuilt.
    - dimensions: contains each dimension and its current state:
        - id: bit offset in the cached bit array
        - dim: the dimension object
        - name: the dimension name
        - selection: the dimension's current selection
    */
    if (!isCrossfilterData(data))
      throw new TypeError("Unsupported crossfilter data type.");
    this.data = data;
    this.selectionCache = selectionCache; /* { BitArray, ... }*/
    this.dimensions = dimensions; /* name: { id, dim, name, selection } */
    Object.preventExtensions(this);
  }

  size(): number {
    return this.data.length;
  }

  all(): RecordType {
    return this.data;
  }

  setData(data: RecordType): ImmutableTypedCrossfilter<RecordType> {
    if (this.data === data) return this;

    if (!isCrossfilterData(data))
      throw new TypeError("Unsupported crossfilter data type.");

    // please leave, WIP
    // console.log("...crossfilter set data, will drop cache");
    return new ImmutableTypedCrossfilter<RecordType>(data, this.dimensions);
  }

  dimensionNames(): string[] {
    /* return array of all dimensions (by name) */
    return Object.keys(this.dimensions);
  }

  hasDimension(name: string): boolean {
    return !!this.dimensions[name];
  }

  addDimension(
    name: string,
    params: CrossfilterDimensionParameters,
  ): ImmutableTypedCrossfilter<RecordType> {
    /*
    Add a new dimension to this crossfilter, of type DimensionType.
    Remainder of parameters are dimension-type-specific.
    */
    const { data } = this;
    const { bitArray } = this.selectionCache;

    if (this.dimensions[name] !== undefined) {
      throw new Error(`Adding duplicate dimension name ${name}`);
    }

    this._clearSelectionCache();

    let id: number | undefined;
    if (bitArray) {
      id = bitArray.allocDimension();
      bitArray.selectAll(id);
    }

    let dim: _ImmutableBaseDimension;
    switch (params.type) {
      case "scalar":
        dim = new ImmutableScalarDimension(name, this.data.length, params);
        break;
      case "enum":
        dim = new ImmutableEnumDimension(name, this.data.length, params);
        break;
      case "spatial":
        dim = new ImmutableSpatialDimension(name, this.data.length, params);
        break;
      default:
        throw new TypeError("Unknown dimension type.");
    }

    Object.freeze(dim);
    const dimensions: Dimensions = {
      ...this.dimensions,
      [name]: {
        id,
        dim,
        name,
        selection: dim.select({ mode: SelectionMode.All }),
      },
    };

    return new ImmutableTypedCrossfilter<RecordType>(data, dimensions, {
      bitArray,
    });
  }

  delDimension(name: string): ImmutableTypedCrossfilter<RecordType> {
    const { data } = this;
    const { bitArray } = this.selectionCache;
    const dimensions = { ...this.dimensions };
    if (dimensions[name] === undefined) {
      throw new ReferenceError(`Unable to delete unknown dimension ${name}`);
    }

    const { id } = dimensions[name];
    delete dimensions[name];
    this._clearSelectionCache();
    if (bitArray) {
      if (id === undefined)
        throw new Error("Internal selection cache inconsistent.");
      bitArray.freeDimension(id);
    }

    return new ImmutableTypedCrossfilter<RecordType>(data, dimensions, {
      bitArray,
    });
  }

  renameDimension(
    oldName: string,
    newName: string,
  ): ImmutableTypedCrossfilter<RecordType> {
    const { [oldName]: dim, ...dimensions } = this.dimensions;
    const { data, selectionCache } = this;
    const newDimensions = {
      ...dimensions,
      [newName]: {
        ...dim,
        name: newName,
        dim: dim.dim.rename(newName),
      },
    };
    return new ImmutableTypedCrossfilter<RecordType>(
      data,
      newDimensions,
      selectionCache,
    );
  }

  select(
    name: string,
    spec: CrossfilterSelector,
  ): ImmutableTypedCrossfilter<RecordType> {
    /*
    select on named dimension, as indicated by `spec`.   Spec is an object
    specifying the selection, and must contain at least a `mode` field.

    Examples:
      select("foo", {mode: "all"});
      select("bar", {mode: "none"});
      select("mumble", {mode: "exact", values: "blue"});
      select("mumble", {mode: "exact", values: ["red", "green", "blue"]});
      select("blort", {mode: "range", lo: 0, hi: 999.99});
    */
    const { data, selectionCache } = this;
    this.selectionCache = {};
    const dimensions = { ...this.dimensions };
    const { dim, id, selection: oldSelection } = dimensions[name];
    const newSelection = dim.select(spec);
    newSelection.ranges = PositiveIntervals.canonicalize(newSelection.ranges);
    dimensions[name] = { id, dim, name, selection: newSelection };
    const newSelectionCache = ImmutableTypedCrossfilter._dimSelnHasUpdated(
      selectionCache,
      id,
      newSelection,
      oldSelection,
    );
    return new ImmutableTypedCrossfilter<RecordType>(
      data,
      dimensions,
      newSelectionCache,
    );
  }

  private static _dimSelnHasUpdated(
    selectionCache: PartialSelectionCache,
    id: number | undefined,
    newSeln: DimensionSelectionResult,
    oldSeln: DimensionSelectionResult,
  ): PartialSelectionCache {
    /*
    Selection has updated from oldSeln to newSeln.   Update the
    bit array if it exists.  If not, we will lazy create it when
    needed.
    */
    if (!selectionCache || !selectionCache.bitArray) return {};

    const { bitArray } = selectionCache;
    if (id === undefined)
      throw new Error("Internal selection cache inconsistent.");

    /*
      if both new and old selection use the same index, we can
      perform an incremental update. If the index changed, we have
      to do a suboptimal full deselect/select.
      */
    let adds: IntervalArray;
    let dels: IntervalArray;
    if (newSeln.index === oldSeln.index) {
      adds = PositiveIntervals.difference(newSeln.ranges, oldSeln.ranges);
      dels = PositiveIntervals.difference(oldSeln.ranges, newSeln.ranges);
    } else {
      // please leave, WIP
      // console.log("suboptimal selection update - index changed");
      adds = newSeln.ranges;
      dels = oldSeln.ranges;
    }

    /*
      allow dimensions to return selected ranges in either dimension sort
      order (indirect via index), or in original record order.

      If sort index exists in the dimension, assume sort ordered ranges.
      */
    const oldIndex = oldSeln.index;
    if (oldIndex) {
      dels.forEach((interval) =>
        bitArray.deselectIndirectFromRange(id, oldIndex, interval),
      );
    } else {
      dels.forEach((interval) => bitArray.deselectFromRange(id, interval));
    }

    const newIndex = newSeln.index;
    if (newIndex) {
      adds.forEach((interval) =>
        bitArray.selectIndirectFromRange(id, newIndex, interval),
      );
    } else {
      adds.forEach((interval) => bitArray.selectFromRange(id, interval));
    }

    return { bitArray };
  }

  private _getSelectionCache(): SelectionCache {
    if (!this.selectionCache) this.selectionCache = {};

    if (!this.selectionCache.bitArray) {
      // console.log("...rebuilding crossfilter cache...");
      const bitArray = new BitArray(this.data.length);
      Object.keys(this.dimensions).forEach((name) => {
        const { selection } = this.dimensions[name];
        const id = bitArray.allocDimension();
        this.dimensions[name].id = id;
        const { ranges, index } = selection;
        ranges.forEach((range) => {
          if (index) {
            bitArray.selectIndirectFromRange(id, index, range);
          } else {
            bitArray.selectFromRange(id, range);
          }
        });
      });
      this.selectionCache.bitArray = bitArray;
    }
    // Guaranteed bitarray is now set.
    return this.selectionCache as SelectionCache;
  }

  private _clearSelectionCache(): PartialSelectionCache {
    this.selectionCache = {};
    return this.selectionCache;
  }

  private _setSelectionCache(
    vals: PartialSelectionCache = {},
  ): PartialSelectionCache {
    Object.assign(this.selectionCache, vals);
    return this.selectionCache;
  }

  allSelected(): CrossfilterData {
    /*
    return array of all records currently selected by all dimensions
    */
    const selectionCache = this._getSelectionCache();
    const { bitArray } = selectionCache;
    const { data } = this;
    if (Array.isArray(data)) {
      const res: unknown[] = [];
      for (let i = 0, len = data.length; i < len; i += 1) {
        if (bitArray.isSelected(i)) {
          res.push(data[i]);
        }
      }
      return res as SortableArray;
    }
    if (data instanceof Dataframe) {
      /* else, Dataframe-like */
      return data.isubsetMask(this.allSelectedMask(), null);
    }
    /* or error! */
    throw new TypeError("Unsupported data type!");
  }

  allSelectedMask(): Uint8Array {
    /*
    return Uint8Array containing selection state (truthy/falsey) for each record.
    */
    const selectionCache = this._getSelectionCache();
    let { allSelectedMask } = selectionCache;

    if (allSelectedMask !== undefined) return allSelectedMask;

    allSelectedMask = selectionCache.bitArray.fillBySelection(
      new Uint8Array(this.data.length),
      1,
      0,
    );
    this._setSelectionCache({ allSelectedMask });
    return allSelectedMask;
  }

  countSelected(): number {
    /*
    return number of records selected on all dimensions
    */
    const selectionCache = this._getSelectionCache();
    let { countSelected } = selectionCache;

    if (countSelected !== undefined) return countSelected;

    countSelected = selectionCache.bitArray.selectionCount();
    this._setSelectionCache({ countSelected });
    return countSelected;
  }

  isElementSelected(i: number): boolean {
    /*
    return truthy/falsey if this record is selected on all dimensions
    */
    const selectionCache = this._getSelectionCache();
    return selectionCache.bitArray.isSelected(i);
  }

  fillByIsSelected<A extends TypedArray>(
    array: A,
    selectedValue: A[0],
    deselectedValue: A[0],
  ): A {
    /*
    fill array with one of two values, based upon selection state.
    */
    const selectionCache = this._getSelectionCache();
    return selectionCache.bitArray.fillBySelection(
      array,
      selectedValue,
      deselectedValue,
    );
  }
}

/*
Base dimension object.

A Dimension is an index, accessed via a select() method.  The protocol
for a dimension:
  - constructor - first param is name, remainder is whatever params are
    required to initialize the dimension.
  - select - one and only param is the selection specifier.  Returns an
    array of record IDs.
  - name - the dimension name/label.
*/
class _ImmutableBaseDimension {
  name: string;

  constructor(name: string) {
    this.name = name;
  }

  clone() {
    return Object.assign(Object.create(Object.getPrototypeOf(this)), this);
  }

  rename(name: string) {
    const d = this.clone();
    d.name = name;
    return d;
  }

  select(spec: CrossfilterSelector): DimensionSelectionResult {
    const { mode } = spec;
    if (mode === undefined) {
      throw new Error("select spec does not contain 'mode'");
    }
    throw new Error(
      `select mode ${mode} not implemented by dimension ${this.name}`,
    );
  }
}

class ImmutableScalarDimension extends _ImmutableBaseDimension {
  protected index: IndexArray;

  protected value: SortableArray;

  // Two constructor modes - caller can provide a pre-created value array,
  // a map function which will create it, or another array which
  // will be copied into the user-specified value array.

  constructor(name: string, length: number, params: ScalarDimensionParameters) {
    super(name);

    const { value, ValueArrayCtor } = params;

    let array;
    if (ValueArrayCtor !== undefined) {
      if (value instanceof ValueArrayCtor) {
        // user has provided the final typed array - just use it
        if (value.length !== length) {
          throw new RangeError(
            "ScalarDimension values length must equal crossfilter data record count",
          );
        }
        array = value;
      } else if (isAnyArray(value)) {
        // Create value array from user-provided array ctor. Typically used
        // only by enumerated dimensions.  May be overridden in subclass.
        array = new ValueArrayCtor(length);
        for (let i = 0; i < length; i += 1) {
          array[i] = value[i];
        }
      }
    }
    if (!array) {
      throw new NotImplementedError(
        "dimension value must be function or value array type",
      );
    }

    this.value = array;

    // create sort index
    this.index = makeSortIndex(array);
  }

  select(spec: CrossfilterSelector): DimensionSelectionResult {
    const { index } = this;
    switch (spec.mode) {
      case SelectionMode.All:
        return { ranges: [[0, this.value.length]], index };
      case SelectionMode.None:
        return { ranges: [], index };
      case SelectionMode.Exact:
        return this.selectExact(spec);
      case SelectionMode.Range:
        return this.selectRange(spec);
      default:
        return super.select(spec);
    }
  }

  selectExact(spec: SelectExact): DimensionSelectionResult {
    const { value, index } = this;
    let { values } = spec;
    if (!Array.isArray(values)) {
      values = [values];
    }
    const ranges: IntervalArray = [];
    for (let v = 0, len = values.length; v < len; v += 1) {
      const r: Interval = [
        lowerBoundIndirect(value, index, values[v], 0, value.length),
        upperBoundIndirect(value, index, values[v], 0, value.length),
      ];
      if (r[0] <= r[1]) {
        ranges.push(r);
      }
    }
    return { ranges, index };
  }

  selectRange(spec: SelectRange): DimensionSelectionResult {
    const { value, index } = this;
    /*
    if !inclusive: [lo, hi) else [lo, hi]
    */
    const { lo, hi, inclusive } = spec;
    const ranges: IntervalArray = [];
    const r: Interval = [
      lowerBoundIndirect(value, index, lo, 0, value.length),
      inclusive
        ? upperBoundIndirect(value, index, hi, 0, value.length)
        : lowerBoundIndirect(value, index, hi, 0, value.length),
    ];
    if (r[0] < r[1]) ranges.push(r);
    return { ranges, index };
  }
}

class ImmutableEnumDimension extends ImmutableScalarDimension {
  private enumIndex: SortableArray;

  constructor(name: string, length: number, params: EnumDimensionParameters) {
    const { value } = params;
    const s = new Set(value);
    const enumIndex = sortArray(Array.from(s));
    const indexArr = new Uint32Array(length);
    for (let i = 0; i < length; i += 1) {
      const v = value[i];
      const e = lowerBound(enumIndex, v, 0, enumIndex.length);
      indexArr[i] = e;
    }

    super(name, length, {
      type: "scalar",
      value: indexArr,
      ValueArrayCtor: Uint32Array,
    });
    this.enumIndex = enumIndex;
  }

  selectExact(spec: SelectExact): DimensionSelectionResult {
    const { enumIndex } = this;
    let { values } = spec;
    if (!Array.isArray(values)) {
      values = [values];
    }
    return super.selectExact({
      mode: spec.mode,
      values: values.map((v) =>
        binarySearch(enumIndex, v, 0, enumIndex.length),
      ),
    });
  }

  selectRange(spec: SelectRange): DimensionSelectionResult {
    throw new Error(
      `select mode ${spec.mode} not implemented by dimension ${this.name}`,
    );
  }
}

class ImmutableSpatialDimension extends _ImmutableBaseDimension {
  private X: TypedArray;

  private Xindex: IndexArray;

  private Y: TypedArray;

  private Yindex: IndexArray;

  constructor(
    name: string,
    length: number,
    params: SpatialDimensionParameters,
  ) {
    super(name);

    const { X, Y } = params;
    if (X.length !== Y.length || X.length !== length) {
      throw new RangeError(
        "SpatialDimension values must have same dimensionality as crossfilter",
      );
    }
    this.X = X;
    this.Y = Y;

    this.Xindex = makeSortIndex(X);
    this.Yindex = makeSortIndex(Y);
  }

  select(spec: CrossfilterSelector): DimensionSelectionResult {
    switch (spec.mode) {
      case SelectionMode.All:
        return { ranges: [[0, this.X.length]], index: null };
      case SelectionMode.None:
        return { ranges: [], index: null };
      case SelectionMode.WithinRect:
        return this.selectWithinRect(spec);
      case SelectionMode.WithinPolygon:
        return this.selectWithinPolygon(spec);
      default:
        return super.select(spec);
    }
  }

  selectWithinRect(spec: SelectWithinRect): DimensionSelectionResult {
    /*
      { mode: "within-rect", minX: 1, minY: 0, maxX: 3, maxY: 9 }
    */
    const { minX, minY, maxX, maxY } = spec;
    const { X, Y } = this;
    const ranges: IntervalArray = [];
    let start = -1;
    for (let i = 0, l = X.length; i < l; i += 1) {
      const x = X[i];
      const y = Y[i];
      const inside = minX <= x && x < maxX && minY <= y && y < maxY;
      if (inside && start === -1) start = i;
      if (!inside && start !== -1) {
        ranges.push([start, i]);
        start = -1;
      }
    }
    if (start !== -1) ranges.push([start, X.length]);
    return { ranges, index: null };
  }

  /*
  Relatively brute force filter by polygon.

  Currently uses d3.polygonContains() to test for polygon inclusion, which itself
  uses a ray casting (crossing number) algorithm.  There are a series of optimizations
  to make this faster:
    * first sliced by X or Y, using an index on the axis
    * then the polygon bounding box is used for trivial rejection
    * then the polygon test is applied
  */

  selectWithinPolygon(spec: SelectWithinPolygon): DimensionSelectionResult {
    /*
      { mode: "within-polygon", polygon: [ [x0, y0], ... ] }
    */
    const { polygon } = spec;
    const [minX, minY, maxX, maxY] = polygonBoundingBox(polygon);
    const { X, Y, Xindex, Yindex } = this;
    const { length } = X;
    let slice;
    let index;
    if (maxY - minY > maxX - minX) {
      slice = [
        lowerBoundIndirect(X, Xindex, minX, 0, length),
        lowerBoundIndirect(X, Xindex, maxX, 0, length),
      ];
      index = Xindex;
    } else {
      slice = [
        lowerBoundIndirect(Y, Yindex, minY, 0, length),
        lowerBoundIndirect(Y, Yindex, maxY, 0, length),
      ];
      index = Yindex;
    }

    const ranges: IntervalArray = [];
    let start = -1;
    for (let i = slice[0], e = slice[1]; i < e; i += 1) {
      const rid = index[i];
      const x = X[rid];
      const y = Y[rid];
      const inside =
        minX <= x &&
        x < maxX &&
        minY <= y &&
        y < maxY &&
        withinPolygon(polygon, x, y);

      if (inside && start === -1) start = i;
      if (!inside && start !== -1) {
        ranges.push([start, i]);
        start = -1;
      }
    }
    if (start !== -1) ranges.push([start, slice[1]]);
    return { ranges, index };
  }
}

/* return bounding box of the polygon */
function polygonBoundingBox(polygon: [number, number][]) {
  let minX = Number.MAX_VALUE;
  let minY = Number.MAX_VALUE;
  let maxX = Number.MIN_VALUE;
  let maxY = Number.MIN_VALUE;
  for (let i = 0, l = polygon.length; i < l; i += 1) {
    const point = polygon[i];
    const [x, y] = point;
    if (x < minX) minX = x;
    if (y < minY) minY = y;
    if (x > maxX) maxX = x;
    if (y > maxY) maxY = y;
  }
  return [minX, minY, maxX, maxY];
}

/**
 *  withinPolygon determines if a point is within a polygon
 *  Code adapted from https://github.com/d3/d3-polygon/blob/master/src/contains.js
 *  @param {array} polygon - is an array of point arrays of format [[x1, y1], [x2, y2], ...]
 *  @param {float} x - point x coordinate
 *  @param {float} y - point y coordinate
 *  @type {boolean}
 */
function withinPolygon(polygon: [number, number][], x: number, y: number) {
  const n = polygon.length;
  let p = polygon[n - 1];
  let x0 = p[0];
  let y0 = p[1];
  let x1;
  let y1;
  let inside = false;

  for (let i = 0; i < n; i += 1) {
    p = polygon[i];
    x1 = p[0];
    y1 = p[1];

    if (y1 > y !== y0 > y && x < ((x0 - x1) * (y - y1)) / (y0 - y1) + x1)
      inside = !inside;
    x0 = x1;
    y0 = y1;
  }
  return inside;
}
