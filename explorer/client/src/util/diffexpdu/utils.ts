import pako from "pako";

export function compressBuffer(
  buf: Uint8Array,
  options: Record<string, unknown> = {},
): Uint8Array {
  return pako.deflateRaw(buf, {
    level: 3,
    windowBits: -15,
    memLevel: 9,
    ...options,
  });
}

export function decompressBuffer(buf: Uint8Array): Uint8Array {
  return pako.inflateRaw(buf);
}

/**
 * Pairwise iterator, eg, [...pairwise([0, 1, 2])] -> [[0,1], [1,2], [2,3]].
 */
export function* pairwise<T>(iterable: Iterable<T>): Generator<Array<T>> {
  const iterator = iterable[Symbol.iterator]();
  let current = iterator.next();
  let next = iterator.next();
  while (!next.done) {
    yield [current.value, next.value];
    current = next;
    next = iterator.next();
  }
}

/**
 * Groupby iterator, grouping by return value of key(). Behavior only defined if
 * iterator is pre-sorted.
 *
 * Example:  [...groupby([0, 0, 1, 2, 2, 3])] -> [[0, 0], [1], [2, 2], [3]]
 */
export function* groupby<T>(
  iterable: Iterable<T>,
  key: (a: T) => any = (a) => a, // eslint-disable-line @typescript-eslint/no-explicit-any -- we actually WANT any...
): Generator<Array<T>> {
  const iterator = iterable[Symbol.iterator]();
  let current = iterator.next();
  let acc = [];
  while (!current.done) {
    acc.push(current.value);
    const next = iterator.next();
    if (next.done || key(current.value) !== key(next.value)) {
      yield acc;
      acc = [];
    }
    current = next;
  }
}

/**
 * Concatenate TypedArray
 */
type TypedArray =
  | Uint8Array
  | Uint16Array
  | Uint32Array
  | Int8Array
  | Int16Array
  | Int32Array
  | Float32Array
  | Float64Array;
type Constructor<T extends unknown = unknown> = new (...args: unknown[]) => T;
export function concat<SrcType extends TypedArray, DstType extends TypedArray>(
  Ctor: Constructor<DstType>,
  tArrays: SrcType[],
): DstType {
  let totalLength = 0;
  for (const arr of tArrays) {
    totalLength += arr.length;
  }
  const dst = new Ctor(totalLength);
  let offset = 0;
  for (const arr of tArrays) {
    dst.set(arr, offset);
    offset += arr.length;
  }
  return dst;
}

export function nonzeroBits(mask: number): number[] {
  const ids = [];
  for (let b = 0; b < 8; b += 1) {
    // eslint-disable-next-line no-bitwise -- a error in a bit op function
    if (mask & (1 << b)) {
      ids.push(b);
    }
  }
  return ids;
}
