import {
  isTypedArray,
  isFloatTypedArray,
  FloatTypedArray,
  SortableArray,
} from "../../common/types/arraytypes";

import { NonFloatSortableArray, IndexArray } from "./types";

/* eslint-disable no-bitwise -- code relies on bitwise ops */

/*
 ** fast sort and search, with separate code paths for floats (NaN ordering),
 ** indirect and direct search/sort.
 */

/*
Comparators for float sort.   -Infinity < finite < Infinity < NaN
*/
function lt(a: number, b: number) {
  if (Number.isNaN(b)) return !Number.isNaN(a);
  return a < b;
}

function gt(a: number, b: number) {
  if (Number.isNaN(a)) return !Number.isNaN(b);
  return a > b;
}

/*
insertion sort, used for small arrays (controlled by SMALL_ARRAY constant)
*/
const SMALL_ARRAY = 32;
function insertionsort<T extends NonFloatSortableArray>(
  a: T,
  lo: number,
  hi: number,
): T {
  for (let i = lo + 1; i < hi + 1; i += 1) {
    const x = a[i];
    let j;
    for (j = i; j > lo && a[j - 1] > x; j -= 1) {
      a[j] = a[j - 1];
    }
    a[j] = x;
  }
  return a;
}

function insertionsortFloats<T extends FloatTypedArray>(
  a: T,
  lo: number,
  hi: number,
): T {
  for (let i = lo + 1; i < hi + 1; i += 1) {
    const x = a[i];
    let j;
    for (j = i; j > lo && gt(a[j - 1], x); j -= 1) {
      a[j] = a[j - 1];
    }
    a[j] = x;
  }
  return a;
}

function insertionsortIndirect(
  a: IndexArray,
  s: NonFloatSortableArray,
  lo: number,
  hi: number,
): IndexArray {
  for (let i = lo + 1; i < hi + 1; i += 1) {
    const x = a[i];
    const t = s[x];
    let j;
    for (j = i; j > lo && s[a[j - 1]] > t; j -= 1) {
      a[j] = a[j - 1];
    }
    a[j] = x;
  }
  return a;
}

function insertionsortFloatsIndirect(
  a: IndexArray,
  s: FloatTypedArray,
  lo: number,
  hi: number,
): IndexArray {
  for (let i = lo + 1; i < hi + 1; i += 1) {
    const x = a[i];
    const t = s[x];
    let j;
    for (j = i; j > lo && gt(s[a[j - 1]], t); j -= 1) {
      a[j] = a[j - 1];
    }
    a[j] = x;
  }
  return a;
}

/*
Quicksort - used for larger arrays
*/
function quicksort<T extends NonFloatSortableArray>(
  a: T,
  lo: number,
  hi: number,
): T {
  if (hi - lo < SMALL_ARRAY) {
    return insertionsort(a, lo, hi);
  }
  if (lo < hi) {
    // partition
    const mid = Math.floor((lo + hi) / 2);
    const p = a[mid];
    let i = lo - 1;
    let j = hi + 1;
    while (i < j) {
      do {
        i += 1;
      } while (a[i] < p);
      do {
        j -= 1;
      } while (a[j] > p);
      if (i < j) {
        const tmp = a[i];
        a[i] = a[j];
        a[j] = tmp;
      }
    }
    // sort
    quicksort(a, lo, j);
    quicksort(a, j + 1, hi);
  }
  return a;
}

function quicksortFloats<T extends FloatTypedArray>(
  a: T,
  lo: number,
  hi: number,
): T {
  if (hi - lo < SMALL_ARRAY) {
    return insertionsortFloats(a, lo, hi);
  }
  if (lo < hi) {
    // partition
    const mid = Math.floor((lo + hi) / 2);
    const p = a[mid];
    let i = lo - 1;
    let j = hi + 1;
    while (i < j) {
      do {
        i += 1;
      } while (lt(a[i], p));
      do {
        j -= 1;
      } while (gt(a[j], p));
      if (i < j) {
        const tmp = a[i];
        a[i] = a[j];
        a[j] = tmp;
      }
    }
    // sort
    quicksortFloats(a, lo, j);
    quicksortFloats(a, j + 1, hi);
  }
  return a;
}

function quicksortIndirect(
  a: IndexArray,
  s: NonFloatSortableArray,
  lo: number,
  hi: number,
): IndexArray {
  if (hi - lo < SMALL_ARRAY) {
    return insertionsortIndirect(a, s, lo, hi);
  }
  if (lo < hi) {
    // partition
    const mid = Math.floor((lo + hi) / 2);
    const p = a[mid];
    const t = s[p];
    let i = lo - 1;
    let j = hi + 1;
    while (i < j) {
      do {
        i += 1;
      } while (s[a[i]] < t);
      do {
        j -= 1;
      } while (s[a[j]] > t);
      if (i < j) {
        const tmp = a[i];
        a[i] = a[j];
        a[j] = tmp;
      }
    }
    // sort
    quicksortIndirect(a, s, lo, j);
    quicksortIndirect(a, s, j + 1, hi);
  }
  return a;
}

function quicksortFloatsIndirect(
  a: IndexArray,
  s: FloatTypedArray,
  lo: number,
  hi: number,
): IndexArray {
  if (hi - lo < SMALL_ARRAY) {
    return insertionsortFloatsIndirect(a, s, lo, hi);
  }
  if (lo < hi) {
    // partition
    const mid = Math.floor((lo + hi) / 2);
    const p = a[mid];
    const t = s[p];
    let i = lo - 1;
    let j = hi + 1;
    while (i < j) {
      do {
        i += 1;
      } while (lt(s[a[i]], t));
      do {
        j -= 1;
      } while (gt(s[a[j]], t));
      if (i < j) {
        const tmp = a[i];
        a[i] = a[j];
        a[j] = tmp;
      }
    }
    // sort
    quicksortFloatsIndirect(a, s, lo, j);
    quicksortFloatsIndirect(a, s, j + 1, hi);
  }
  return a;
}

/*
Convenience wrappers, handling optimization paths and default
handlers for NaN comparisons.  Sorts in place.
*/
export function sortArray<T extends SortableArray>(arr: T): T {
  if (Array.isArray(arr)) {
    return quicksort(arr, 0, arr.length - 1);
  }
  if (isTypedArray(arr)) {
    if (isFloatTypedArray(arr)) {
      return quicksortFloats(arr, 0, arr.length - 1);
    }
    return quicksort(arr, 0, arr.length - 1);
  }
  /* else unsupported */
  throw new Error("sortArray received unsupported object type");
}

export function sortIndex(
  index: IndexArray,
  source: SortableArray,
): IndexArray {
  if (isFloatTypedArray(source))
    return quicksortFloatsIndirect(index, source, 0, index.length - 1);
  return quicksortIndirect(index, source, 0, index.length - 1);
}

// Search for `value` in the sorted array `arr`, in the range [first, last).
// Return the first (left most) index where arr[index] >= value.
//
// In other words, return array index I where:
//  arr[i] < value for all tarr[lo:I]
//  arr[i] >= value for all tarr[I:last]
//
// The same semantics/behavior as:
//    C++: lower_bound()
//    Python: bisect.bisect_left()
//
function lowerBoundNonFloat<T extends NonFloatSortableArray>(
  valueArray: T,
  value: T[0],
  first: number,
  last: number,
): number {
  let lfirst = first;
  let llast = last;
  // this is just a binary search
  while (lfirst < llast) {
    const middle = (lfirst + llast) >>> 1;
    if (valueArray[middle] < value) {
      lfirst = middle + 1;
    } else {
      llast = middle;
    }
  }
  return lfirst;
}

// lowerBound, but with NaN handling
//
// If the underlying array is a Float32Array or Float64Array, will enforce
// the ordering -Infinity < finite < Infinity < NaN.
//
function lowerBoundFloat(
  valueArray: FloatTypedArray,
  value: number,
  first: number,
  last: number,
): number {
  let lfirst = first;
  let llast = last;
  // this is just a binary search
  while (lfirst < llast) {
    const middle = (lfirst + llast) >>> 1;
    if (lt(valueArray[middle], value)) {
      lfirst = middle + 1;
    } else {
      llast = middle;
    }
  }
  return lfirst;
}

export function lowerBound<T extends SortableArray>(
  valueArray: T,
  value: T[0],
  first: number,
  last: number,
): number {
  if (isFloatTypedArray(valueArray)) {
    // TSC is not able to infer that T[0] is a number, so we revert to a cast
    return lowerBoundFloat(valueArray, value as number, first, last);
  }
  return lowerBoundNonFloat(valueArray, value, first, last);
}

// Inlined performance optimization - used to indirect through a sort map.
//
function lowerBoundNonFloatIndirect<T extends NonFloatSortableArray>(
  valueArray: T,
  indexArray: IndexArray,
  value: T[0],
  first: number,
  last: number,
): number {
  let lfirst = first;
  let llast = last;
  // this is just a binary search
  while (lfirst < llast) {
    const middle = (lfirst + llast) >>> 1;
    if (valueArray[indexArray[middle]] < value) {
      lfirst = middle + 1;
    } else {
      llast = middle;
    }
  }
  return lfirst;
}

function lowerBoundFloatIndirect(
  valueArray: FloatTypedArray,
  indexArray: IndexArray,
  value: number,
  first: number,
  last: number,
): number {
  let lfirst = first;
  let llast = last;
  // this is just a binary search
  while (lfirst < llast) {
    const middle = (lfirst + llast) >>> 1;
    if (lt(valueArray[indexArray[middle]], value)) {
      lfirst = middle + 1;
    } else {
      llast = middle;
    }
  }
  return lfirst;
}

export function lowerBoundIndirect<T extends SortableArray>(
  valueArray: T,
  indexArray: IndexArray,
  value: T[0],
  first: number,
  last: number,
): number {
  if (isFloatTypedArray(valueArray)) {
    return lowerBoundFloatIndirect(
      valueArray,
      indexArray,
      value as number, // TSC is not able to infer that T[0] is a number, so we revert to a cast
      first,
      last,
    );
  }
  return lowerBoundNonFloatIndirect(valueArray, indexArray, value, first, last);
}

// Search for `value in the sorted array `arr`, in the range [first, last).
// Return the first value where arr[index] > value.
//
// In other words, return array index I, where:
//  arr[i] <= value for all tarr[lo:I]
//  arr[i] > value for all tarr[I:last]
//
// The same semantics/behavior as:
//    C++: upper_bound()
//    Python: bisect.bisect_right()
//
function upperBoundNonFloat<T extends NonFloatSortableArray>(
  valueArray: T,
  value: T[0],
  first: number,
  last: number,
): number {
  let lfirst = first;
  let llast = last;
  // this is just a binary search
  while (lfirst < llast) {
    const middle = (lfirst + llast) >>> 1;
    if (valueArray[middle] > value) {
      llast = middle;
    } else {
      lfirst = middle + 1;
    }
  }
  return lfirst;
}

function upperBoundFloat(
  valueArray: FloatTypedArray,
  value: number,
  first: number,
  last: number,
): number {
  let lfirst = first;
  let llast = last;
  // this is just a binary search
  while (lfirst < llast) {
    const middle = (lfirst + llast) >>> 1;
    if (gt(valueArray[middle], value)) {
      llast = middle;
    } else {
      lfirst = middle + 1;
    }
  }
  return lfirst;
}

export function upperBound<T extends SortableArray>(
  valueArray: T,
  value: T[0],
  first: number,
  last: number,
): number {
  if (isFloatTypedArray(valueArray)) {
    // TSC is not able to infer that T[0] is a number, so we revert to a cast
    return upperBoundFloat(valueArray, value as number, first, last);
  }
  return upperBoundNonFloat(valueArray, value, first, last);
}

// Inline performance optimization
//
function upperBoundNonFloatIndirect<T extends NonFloatSortableArray>(
  valueArray: T,
  indexArray: IndexArray,
  value: T[0],
  first: number,
  last: number,
): number {
  let lfirst = first;
  let llast = last;
  // this is just a binary search
  while (lfirst < llast) {
    const middle = (lfirst + llast) >>> 1;
    if (valueArray[indexArray[middle]] > value) {
      llast = middle;
    } else {
      lfirst = middle + 1;
    }
  }
  return lfirst;
}

function upperBoundFloatIndirect(
  valueArray: FloatTypedArray,
  indexArray: IndexArray,
  value: number,
  first: number,
  last: number,
): number {
  let lfirst = first;
  let llast = last;
  // this is just a binary search
  while (lfirst < llast) {
    const middle = (lfirst + llast) >>> 1;
    if (gt(valueArray[indexArray[middle]], value)) {
      llast = middle;
    } else {
      lfirst = middle + 1;
    }
  }
  return lfirst;
}

export function upperBoundIndirect<T extends SortableArray>(
  valueArray: T,
  indexArray: IndexArray,
  value: T[0],
  first: number,
  last: number,
): number {
  if (isFloatTypedArray(valueArray)) {
    // TSC is not able to infer that T[0] is a number, so we revert to a cast
    return upperBoundFloatIndirect(
      valueArray,
      indexArray,
      value as number,
      first,
      last,
    );
  }
  return upperBoundNonFloatIndirect(valueArray, indexArray, value, first, last);
}

// Search for `value` in the sorted array `arr`, in the range [first, last).
// Return the first index where arr[index] == value, OR if value not present,
// return `last`
//
// The same semantics/behavior as:
//    C++: binary_search()
//
export function binarySearch<T extends SortableArray>(
  valueArray: T,
  value: T[0],
  first: number,
  last: number,
): number {
  const index = lowerBound(valueArray, value, first, last);
  if (index !== last && value === valueArray[index]) return index;
  return last;
}
/* eslint-enable no-bitwise -- enable */
