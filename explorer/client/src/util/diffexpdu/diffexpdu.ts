import { PduBuffer } from "./pdubuffer";
import { packPostingsLists, unpackPostingsLists } from "./postingslist";

/*
format is documented in dev_docs/diffexpdu.md
*/
const MAGIC_NUMBER = 0xde;

/**
 * The diffex mode. Currently, only TopN is implemented.
 */
export enum DiffExMode {
  TopN = 0,
  VarFilter = 1,
}

/**
 * The diffex arguments. IMPORTANT:
 *  - set1 and set2 must be monotonically increasing (sorted, unique)
 *  - set1 and set2 must be disjoint (no elements in common)
 */
export interface DiffExArguments {
  mode: DiffExMode;
  params: {
    N: number;
  };
  set1: Uint32Array;
  set2: Uint32Array;
}

/**
 * Encode (pack) differential expression arguments into a (smaller) binary buffer.
 *
 * @param {DiffExArguments} args - the differential expression arguments
 * @returns {Uint8Array} the binary encoded buffer
 */
export function packDiffExPdu(args: DiffExArguments): Uint8Array {
  if (args.mode !== DiffExMode.TopN)
    throw new Error("Modes other than TopN are unsupported.");
  if (args.set1.length === 0 || args.set2.length === 0)
    throw new Error("Cell sets must be nonzero length");

  const pdu = new PduBuffer(1024);

  // diffex header
  packDiffExHeader(pdu, args.mode, args.params.N);

  // cell set 1 and 2
  packPostingsLists(pdu, [args.set1, args.set2], true);

  return pdu.asUint8Array();
}

/**
 * Unpack a buffer containing binary encoded differential expression arguments.
 *
 * @param {Uint8Array} buf - the buffer containing the encoded arguments.
 * @returns the decoded arguments
 */
export function unpackDiffExPdu(buf: Uint8Array): DiffExArguments {
  const dv = new DataView(buf.buffer, buf.byteOffset, buf.byteLength);
  const { mode, N } = unpackDiffExHeader(dv);
  const [set1, set2] = unpackPostingsLists(buf.subarray(4));
  return {
    mode,
    params: {
      N,
    },
    set1,
    set2,
  };
}

function packDiffExHeader(pdu: PduBuffer, mode: number, N: number): void {
  pdu.addUint8(MAGIC_NUMBER);
  pdu.addUint8(mode);
  pdu.addUint16(N);
}

function unpackDiffExHeader(dv: DataView): { mode: number; N: number } {
  const magic = dv.getUint8(0);
  if (magic !== MAGIC_NUMBER)
    throw new Error("Malformed diffexpdu - incorrect magic number");
  const mode = dv.getUint8(1);
  const N = dv.getUint16(2, true);
  return { mode, N };
}
