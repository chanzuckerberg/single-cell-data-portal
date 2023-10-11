/* eslint-disable no-bitwise -- relies on bitwise ops */

/*
The format is documented in `dev_docs/diffexpdu.md`. It is implemented in both Python
and Typescript.  Implementations must maintain compatibility.
*/

import { PduBuffer } from "./pdubuffer";
import {
  compressBuffer,
  decompressBuffer,
  nonzeroBits,
  pairwise,
  groupby,
  concat,
} from "./utils";

interface ListPartition {
  listId: number;
  key: number;
  startIdx: number;
  endIdx: number;
  arr: Uint32Array;
}

const MAGIC_NUMBER = 0xce;

enum BlockType {
  BitArray = 0,
  Uint16Array = 1,
  Uint16ArrayInverted = 2,
}

interface BlockDescription {
  blockType: BlockType;
  listIdMask: number;
  nElem: number;
  key: number;
  nBytes: number;
}

export function packPostingsLists(
  pdu: PduBuffer,
  postingsLists: Array<Uint32Array>,
  sorted = false,
): void {
  if (postingsLists.length === 0) {
    throw new Error("must specify one or more postings lists");
  }
  if (!sorted) {
    throw new Error("all postings lists must be sorted");
  }

  const allPartitions = partitionLists(postingsLists);
  allPartitions.sort((a, b) => a.key - b.key);

  const encodedBlocks = [];
  for (const partition of allPartitions) {
    encodedBlocks.push(encodePartition(partition));
  }

  // encode list header
  const blockCount = encodedBlocks.length;
  packPostingsListsHeader(pdu, 1, blockCount);

  // encode each block description
  for (const [header] of encodedBlocks) {
    packBlockDescription(pdu, header);
  }

  // encode each block
  for (const [, block] of encodedBlocks) {
    packBlock(pdu, block);
  }
}

function packPostingsListsHeader(
  pdu: PduBuffer,
  nLists: number,
  nBlocks: number,
): void {
  pdu.ensureFreeCapacity(4);
  pdu.addUint8(MAGIC_NUMBER);
  pdu.addUint8(nLists - 1);
  pdu.addUint16(nBlocks - 1);
}

function packBlockDescription(
  pdu: PduBuffer,
  description: BlockDescription,
): void {
  pdu.ensureFreeCapacity(8);
  pdu.addUint8(description.blockType);
  pdu.addUint8(description.listIdMask);
  pdu.addUint16(description.nElem - 1);
  pdu.addUint16(description.key);
  pdu.addUint16(description.nBytes - 1);
}

function packBlock(pdu: PduBuffer, block: Uint8Array): void {
  pdu.ensureFreeCapacity(block.byteLength);
  pdu.add(block);
}

/**
 * Given a partition, choose a block type which is likely near-optimal for
 * encoding.  Heuristics live here.
 *
 * A description of the heuristics and their rationale is documented
 * in the format description (`dev_docs/diffexpdu.md`)
 *
 * @param partition
 * @returns the type of block
 */
function chooseBlockType(partition: ListPartition): BlockType {
  const nElem = partition.endIdx - partition.startIdx;
  const interval =
    (partition.arr[partition.endIdx - 1] & 0xffff) +
    1 -
    (partition.arr[partition.startIdx] & 0xffff);
  const density = nElem / interval;

  const minBitArrayThreshold = 2 ** 11;
  const maxBitArrayThreshold = 2 ** 16 - minBitArrayThreshold;
  const minDensity = 0.125;
  const maxDensity = 0.875;
  if (
    minBitArrayThreshold < nElem &&
    nElem < maxBitArrayThreshold &&
    minDensity < density &&
    density < maxDensity
  )
    return BlockType.BitArray;

  return BlockType.Uint16Array;
}

function encodePartition(
  partition: ListPartition,
): [BlockDescription, Uint8Array] {
  const arr = partition.arr.subarray(partition.startIdx, partition.endIdx);
  const blockType = chooseBlockType(partition);
  if (blockType === BlockType.BitArray) {
    return encodeBitarrayBlock(partition, arr);
  }
  return encodeUint16Block(partition, arr);
}

function encodeBitarrayBlock(
  partition: ListPartition,
  list: Uint32Array,
): [BlockDescription, Uint8Array] {
  const buf = new PduBuffer(2 ** 13);
  buf.advance(2 ** 13);
  const buf32 = buf.asUint32Array();
  for (let i = 0; i < list.length; i += 1) {
    const val = list[i] & 0xffff;
    buf32[(val / 32) >> 0] |= 1 << (val & 0x1f);
  }
  const compressedBuf = compressBuffer(buf.asUint8Array());
  const header = {
    blockType: BlockType.BitArray,
    listIdMask: 1 << partition.listId,
    nElem: list.length,
    key: list[0] >> 16,
    nBytes: compressedBuf.length,
  };
  return [header, compressedBuf];
}

function encodeUint16Block(
  partition: ListPartition,
  list: Uint32Array,
): [BlockDescription, Uint8Array] {
  const nElem = list.length;
  const key = list[0] >> 16;
  let blockType: BlockType;
  let buf: Uint8Array;

  const nElemInverted = intervalInvertedLength(list);

  if (nElemInverted + 2 >= nElem) {
    // Encode as-is.
    blockType = BlockType.Uint16Array;
    buf = byteShuffle16(delta16(list));
  } else {
    // Encode the inverted values as it will be shorter.
    // inverted array is prefaced by [intervalStart, intervalEnd]
    blockType = BlockType.Uint16ArrayInverted;
    buf = new Uint8Array(2 * nElemInverted + 4);
    const dv = new DataView(buf.buffer);
    dv.setUint16(0, list[0] & 0xffff, true); // interval start
    dv.setUint16(2, list[list.length - 1] & 0xffff, true); // interval end INCLUSIVE
    byteShuffle16(delta16(intervalInvert(list)), buf.subarray(4));
  }

  const compressedBuf = compressBuffer(buf);
  const header = {
    blockType,
    listIdMask: 1 << partition.listId,
    nElem,
    key,
    nBytes: compressedBuf.length,
  };
  return [header, compressedBuf];
}

export function unpackPostingsLists(buf: Uint8Array): Array<Uint32Array> {
  const dv = new DataView(buf.buffer, buf.byteOffset, buf.byteLength);
  const [startOffset, , , blockDescriptions] = unpackPostingsListsHeader(dv);

  let offset = startOffset;
  const subLists: [number, Uint32Array][] = [];
  for (const blockDescription of blockDescriptions) {
    const subList = unpackBlock(blockDescription, buf.subarray(offset));
    subLists.push(...subList);
    offset += blockDescription.nBytes;
  }

  // now group by listId and concatenate
  subLists.sort((a, b) => a[0] - b[0]);
  const concatLists = [];
  for (const listGroup of groupby(subLists, (sl) => sl[0])) {
    concatLists.push(
      concat(
        Uint32Array,
        listGroup.map((lg) => lg[1]),
      ),
    );
  }
  return concatLists;
}

function unpackPostingsListsHeader(
  dv: DataView,
): [number, number, number, Array<BlockDescription>] {
  const magic = dv.getUint8(0);
  if (magic !== MAGIC_NUMBER)
    throw new Error("Malformed postings list header - incorrect magic number");
  const nLists = dv.getUint8(1) + 1;
  const nBlocks = dv.getUint16(2, true) + 1;
  let offset = 4;

  const blockDescriptions = [];
  for (let i = 0; i < nBlocks; i += 1) {
    const blockType = dv.getUint8(offset + 0);
    const listIdMask = dv.getUint8(offset + 1);
    const nElem = dv.getUint16(offset + 2, true) + 1;
    const key = dv.getUint16(offset + 4, true);
    const nBytes = dv.getUint16(offset + 6, true) + 1;
    blockDescriptions.push({
      blockType,
      listIdMask,
      nElem,
      key,
      nBytes,
    });
    offset += 8;
  }
  return [offset, nLists, nBlocks, blockDescriptions];
}

function unpackBlock(
  blockDescription: BlockDescription,
  buf: Uint8Array,
): Array<[number, Uint32Array]> {
  const { blockType } = blockDescription;
  const decoded: [number, Uint32Array][] = []; // array of [listId, list]
  const listIds = nonzeroBits(blockDescription.listIdMask);
  // Current spec only supports one list per block.
  if (listIds.length !== 1) throw new Error("Unexpected multi-list block.");

  buf = decompressBuffer(buf);

  switch (blockType) {
    case BlockType.BitArray: {
      decoded.push(decodeBitarrayBlock(blockDescription, buf, listIds));
      break;
    }
    case BlockType.Uint16Array: {
      decoded.push(decodeUint16Block(blockDescription, buf, listIds));
      break;
    }
    case BlockType.Uint16ArrayInverted: {
      decoded.push(decodeUint16InvertedBlock(blockDescription, buf, listIds));
      break;
    }
    default:
      throw new Error("Unknown block type");
  }

  // Add key to each block result
  const key = blockDescription.key << 16;
  const { nElem } = blockDescription;
  for (const [, list] of decoded) {
    if (list.length !== nElem)
      throw new Error("Block corrupted - nElem mismatch");
    for (let i = 0; i < nElem; i += 1) {
      list[i] += key;
    }
  }

  return decoded;
}

function decodeBitarrayBlock(
  blockDescription: BlockDescription,
  buf: Uint8Array,
  listIds: number[],
): [number, Uint32Array] {
  const { nElem } = blockDescription;
  const buf32 = new Uint32Array(buf.buffer, buf.byteOffset);
  const list = new Uint32Array(blockDescription.nElem);
  for (let i = 0, elem = 0; i < 2 ** 16 && elem < nElem; i += 1) {
    if (buf32[(i / 32) >> 0] & (1 << (i & 0x1f))) {
      list[elem] = i;
      elem += 1;
    }
  }
  return [listIds[0], list];
}

function decodeUint16Block(
  _: BlockDescription,
  buf: Uint8Array,
  listIds: number[],
): [number, Uint32Array] {
  return [listIds[0], unDelta16(unByteShuffle16(buf))];
}

function decodeUint16InvertedBlock(
  _: BlockDescription,
  buf: Uint8Array,
  listIds: number[],
): [number, Uint32Array] {
  const dv = new DataView(buf.buffer, buf.byteOffset, buf.byteLength);
  const interval: [number, number] = [
    dv.getUint16(0, true),
    dv.getUint16(2, true) + 1,
  ];
  const list = intervalInvert(
    unDelta16(unByteShuffle16(buf.subarray(4))),
    interval,
  );
  return [listIds[0], list];
}

function partitionLists(plists: Array<Uint32Array>): Array<ListPartition> {
  const partitions = [];
  for (let listId = 0; listId < plists.length; listId += 1) {
    const plist = plists[listId];
    const boundaries = findPartitionBoundaries(plist);
    for (const [startIdx, endIdx] of pairwise(boundaries)) {
      partitions.push({
        listId,
        key: plist[startIdx] >> 16,
        startIdx,
        endIdx,
        arr: plist,
      });
    }
  }
  return partitions;
}

function findPartitionBoundaries(arr: Uint32Array): Array<number> {
  const partitions: number[] = [];
  if (arr.length === 0) return partitions;

  let currentPartition = -1;
  for (let idx = 0; idx < arr.length; idx += 1) {
    const val = arr[idx];
    const partition = val >> 16;
    if (partition !== currentPartition) {
      partitions.push(idx);
      currentPartition = val >> 16;
    }
  }
  partitions.push(arr.length); // terminal
  return partitions;
}

/**
 * Compute delta between each element, for lowest 16 bits ONLY.
 *
 * @param src - Uint32Array, sorted, with common top 16 bits.
 * @returns Uint16Array, with difference by position
 */
function delta16(src: Uint32Array): Uint16Array {
  const nElem = src.length;
  const dst = new Uint16Array(nElem);
  let prev = src[0] & 0xffff0000;
  let i: number;

  // unroll loop by blocks of eight
  for (i = 0; i < (nElem & ~7); i += 8) {
    dst[i + 7] = src[i + 7] - src[i + 6];
    dst[i + 6] = src[i + 6] - src[i + 5];
    dst[i + 5] = src[i + 5] - src[i + 4];
    dst[i + 4] = src[i + 4] - src[i + 3];
    dst[i + 3] = src[i + 3] - src[i + 2];
    dst[i + 2] = src[i + 2] - src[i + 1];
    dst[i + 1] = src[i + 1] - src[i + 0];
    dst[i + 0] = src[i + 0] - prev;
    prev = src[i + 7];
  }
  for (; i < nElem; i += 1) {
    dst[i] = src[i] - prev;
    prev = src[i];
  }
  return dst;
}

function unDelta16(src: Uint16Array): Uint32Array {
  const nElem = src.length;
  const dst = new Uint32Array(nElem);
  let prev = 0;
  for (let i = 0; i < nElem; i += 1) {
    dst[i] = prev + src[i];
    prev = dst[i];
  }
  return dst;
}

function byteShuffle16(src: Uint16Array, dst?: Uint8Array): Uint8Array {
  const nElem = src.length;
  const srcU8 = new Uint8Array(src.buffer);
  if (dst === undefined) dst = new Uint8Array(2 * nElem);

  // unroll loop by blocks of four
  let i = 0;
  for (; i < (nElem & ~3); i += 4) {
    dst[i + 0] = srcU8[(i + 0) * 2];
    dst[i + 1] = srcU8[(i + 1) * 2];
    dst[i + 2] = srcU8[(i + 2) * 2];
    dst[i + 3] = srcU8[(i + 3) * 2];

    dst[nElem + (i + 0)] = srcU8[(i + 0) * 2 + 1];
    dst[nElem + (i + 1)] = srcU8[(i + 1) * 2 + 1];
    dst[nElem + (i + 2)] = srcU8[(i + 2) * 2 + 1];
    dst[nElem + (i + 3)] = srcU8[(i + 3) * 2 + 1];
  }
  for (; i < nElem; i += 1) {
    dst[i] = srcU8[i * 2];
    dst[nElem + i] = srcU8[i * 2 + 1];
  }

  return dst;
}

function unByteShuffle16(src: Uint8Array, dst?: Uint8Array): Uint16Array {
  const typeSize = 2; // Uint16
  const nElem = src.length / typeSize;
  if (dst === undefined) dst = new Uint8Array(src.length);

  for (let i = 0; i < nElem; i += 1) {
    for (let j = 0; j < typeSize; j += 1) {
      dst[i * typeSize + j] = src[j * nElem + i];
    }
  }
  return new Uint16Array(dst.buffer);
}

/**
 * For the given increasing monotone `arr`, return the length of
 * the array required to store the array if it was inverted across
 * the given interval.
 *
 * Examples:
 *    [0, 4] -> 3
 *    [0, 1, 3] -> 1
 *
 * @param arr - Uint32Array, containing an increasing monotone.
 * @param interval - the interval across which to invert the values. If not specified, defaults to the range of values in `arr`.
 */
function intervalInvertedLength(
  arr: Uint32Array,
  interval?: [number, number],
): number {
  if (interval === undefined) {
    interval = [arr[0], arr[arr.length - 1] + 1];
  }
  return interval[1] - interval[0] - arr.length;
}

/**
 * Invert the array values over the given interval. Array must be sorted monotone.
 * If interval not specified, will use [arr[0], arr[-1]+1]
 * @param src - source array
 * @param interval - [start, end)
 * @returns inverted array
 */
function intervalInvert(
  src: Uint32Array,
  interval?: [number, number],
): Uint32Array {
  if (interval === undefined) {
    interval = [src[0], src[src.length - 1] + 1];
  }
  const dst = new Uint32Array(intervalInvertedLength(src, interval));
  for (let i = interval[0], srcIdx = 0, dstIdx = 0; i < interval[1]; i += 1) {
    if (srcIdx < src.length && src[srcIdx] === i) {
      srcIdx += 1;
    } else {
      dst[dstIdx] = i;
      dstIdx += 1;
    }
  }
  return dst;
}

/* eslint-enable no-bitwise -- relies on bitwise ops */
