/**
 * This module exports a subclass of Int8Array for
 * categorical data types. The subclass contains the
 * `codeMapping` dictionary which maps codes to values,
 * and a helper function `vat` that returns the value at
 * a particular index.
 */
import { CodeMapping } from "./code_mapping_interfaces";

export class DictEncoded8Array extends Int8Array {
  codeMapping: CodeMapping;

  constructor(props: ArrayBufferLike) {
    super(props);
    this.codeMapping = {};
  }

  setCodeMapping = (codeMap: CodeMapping): void => {
    this.codeMapping = codeMap;
  };

  vat = (index: number): string => this.codeMapping?.[this[index]] ?? "NaN";
}
