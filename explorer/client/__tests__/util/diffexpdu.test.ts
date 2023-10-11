import {
  packDiffExPdu,
  unpackDiffExPdu,
  DiffExMode,
  DiffExArguments,
} from "../../src/util/diffexpdu";

import { range } from "../../src/util/range";

describe("diffexpdu", () => {
  test("roundtrip diffex pdu", () => {
    // single PDU round-trip
    const deArgs: DiffExArguments = {
      mode: DiffExMode.TopN,
      params: { N: 50 },
      set1: new Uint32Array(range(0, 100, 2)),
      set2: new Uint32Array(range(1, 100, 2)),
    };
    const buf = packDiffExPdu(deArgs);
    const decoded = unpackDiffExPdu(buf);
    expect(decoded).toMatchObject(deArgs);
  });

  test("vary density", () => {
    // Test wide variety of densities
    const densityChoices = [
      0.0001, 0.01, 0.1, 0.2, 0.4, 0.6, 0.8, 0.99, 0.9999,
    ];

    const nobs = 4 * 2 ** 16;
    for (const set1Density of densityChoices) {
      const nElem1 = Math.floor(nobs * set1Density);

      for (const set2Density of densityChoices) {
        const nElem2 = Math.floor(nobs * set2Density);
        const arr = shuffledIntArray(0, nobs, 1);
        const set1 = arr.subarray(0, nElem1);
        const set2 = arr.subarray(nElem1, nElem1 + nElem2);
        set1.sort();
        set2.sort();

        const deArgs: DiffExArguments = {
          mode: DiffExMode.TopN,
          params: { N: Math.floor(Math.random() * 100) },
          set1,
          set2,
        };

        const buf = packDiffExPdu(deArgs);
        const decoded = unpackDiffExPdu(buf);
        expect(decoded).toMatchObject(deArgs);
      }
    }
  });
});

// Fisher-Yates shuffled array.  Returns a shuffled array
// from start to stop, by step.
//
function shuffledIntArray(
  start: number,
  stop: number,
  step: number,
): Uint32Array {
  // Create and fill array with values
  const len = Math.max(Math.ceil((stop - start) / step), 0);
  const arr = new Uint32Array(len);
  for (let idx = 0, val = start; idx < len; idx += 1, val += step) {
    arr[idx] = val;
  }
  // Shuffle
  let i = len;
  while (i) {
    const j = Math.floor(Math.random() * i);
    i -= 1;
    const temp = arr[i];
    arr[i] = arr[j];
    arr[j] = temp;
  }
  return arr;
}
