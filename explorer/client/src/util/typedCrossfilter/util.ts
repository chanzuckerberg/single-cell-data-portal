import { sortIndex } from "./sort";
import { rangeFill as fillRange } from "../range";
import { SortableArray } from "../../common/types/arraytypes";
import { IndexArray } from "./types";

/*
    Utility functions, private to this module.
*/

export function makeSortIndex(src: SortableArray): IndexArray {
  const index = fillRange(new Uint32Array(src.length));
  sortIndex(index, src);
  return index;
}

export class NotImplementedError extends Error {
  constructor(msg: string) {
    super(msg);

    // Maintains proper stack trace for where our error was thrown (only available on V8)
    if (Error.captureStackTrace) {
      Error.captureStackTrace(this, NotImplementedError);
    }
  }
}
