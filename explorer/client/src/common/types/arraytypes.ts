/**
 * Utility type and interface definitions.
 */

import { DictEncoded16Array } from "../../util/stateManager/dict-encoded16_array";
import { DictEncoded8Array } from "../../util/stateManager/dict-encoded8_array";
import { DictEncoded32Array } from "../../util/stateManager/dict-encoded32_array";

/**
 * Arrays of numeric values.
 */
export type TypedArray =
  | Int8Array
  | Uint8Array
  | Uint8ClampedArray
  | Int16Array
  | Uint16Array
  | Int32Array
  | Uint32Array
  | Float32Array
  | Float64Array;

export type UnsignedIntTypedArray =
  | Uint8Array
  | Uint8ClampedArray
  | Uint16Array
  | Uint32Array;
export type FloatTypedArray = Float32Array | Float64Array;
export type IntTypedArray = Int8Array | Int16Array | Int32Array;

export type TypedArrayConstructor =
  | Int8ArrayConstructor
  | Uint8ArrayConstructor
  | Int16ArrayConstructor
  | Uint16ArrayConstructor
  | Int32ArrayConstructor
  | Uint32ArrayConstructor
  | Float32ArrayConstructor
  | Float64ArrayConstructor;

/** this general type is the basis of all data stored in dataframe/crossfilter */
export type SortableArray = Array<number | string | boolean> | TypedArray;

/** used if you accept array-ish, but don't care about the contents of the object */
export type AnyArray = Array<unknown> | TypedArray;

export interface GenericArrayConstructor<T> {
  new (
    ...args: ConstructorParameters<
      typeof Int8Array &
        typeof Uint8Array &
        typeof Uint8ClampedArray &
        typeof Int16Array &
        typeof Uint16Array &
        typeof Int32Array &
        typeof Uint32Array &
        typeof Float32Array &
        typeof Float64Array &
        typeof DictEncoded8Array &
        typeof DictEncoded16Array &
        typeof DictEncoded32Array &
        typeof Array &
        T
    >
  ): T;
}

export type NumberArray = Array<number> | TypedArray;
export type NumberArrayConstructor = GenericArrayConstructor<
  Array<number> | TypedArray
>;

export type Int8 = Int8Array[0];
export type Uint8 = Uint8Array[0];
export type Int16 = Int16Array[0];
export type Uint16 = Uint16Array[0];
export type Int32 = Int32Array[0];
export type Uint32 = Uint32Array[0];
export type Float32 = Float32Array[0];
export type Float64 = Float64Array[0];

/**
 * Test if the parameter is a TypedArray.
 * @param tbd - value to be tested
 * @returns true if `tbd` is a TypedArray, false if not.
 */
export function isTypedArray(tbd: unknown): tbd is TypedArray {
  return (
    ArrayBuffer.isView(tbd) &&
    Object.prototype.toString.call(tbd) !== "[object DataView]"
  );
}

/**
 * Test if the paramter is a float TypedArray
 * @param tbd - value to be tested
 * @returns - true if `tbd` is a float typed array.
 */
export function isFloatTypedArray(tbd: unknown): tbd is FloatTypedArray {
  return tbd instanceof Float32Array || tbd instanceof Float64Array;
}

/**
 * Test if the paramter is an unsigned int TypedArray
 * @param tbd - value to be tested
 * @returns - true if `tbd` is an unsigned int typed array.
 */
export function isUnsignedIntTypedArray(
  tbd: unknown,
): tbd is UnsignedIntTypedArray {
  return (
    tbd instanceof Uint8Array ||
    tbd instanceof Uint8ClampedArray ||
    tbd instanceof Uint16Array ||
    tbd instanceof Uint32Array
  );
}

/**
 * Test if the parameter is a categorical typed array
 * @param c - value to be tested
 * @returns - true if `c` is a dictionary encoded array
 */
export const isDictEncodedTypedArray = (c: unknown): c is DictEncodedArray =>
  c instanceof DictEncoded8Array ||
  c instanceof DictEncoded16Array ||
  c instanceof DictEncoded32Array;

export type DictEncodedArray =
  | DictEncoded8Array
  | DictEncoded16Array
  | DictEncoded32Array;

/**
 * Test if the paramter is an int TypedArray
 * @param tbd - value to be tested
 * @returns - true if `tbd` is an int typed array.
 */
export function isIntTypedArray(tbd: unknown): tbd is IntTypedArray {
  return (
    tbd instanceof Int8Array ||
    tbd instanceof Int16Array ||
    tbd instanceof Int32Array
  );
}

/**
 * Test if the parameter is a TypedArray or Array
 * @param tbd - value to be tested
 * @returns - true if `tbd` is a TypedArray or Array
 */
export function isAnyArray(tbd: unknown): tbd is AnyArray {
  return Array.isArray(tbd) || isTypedArray(tbd);
}
