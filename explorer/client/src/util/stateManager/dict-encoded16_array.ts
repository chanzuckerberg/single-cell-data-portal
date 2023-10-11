/**
 * This module exports a subclass of Int16Array for
 * categorical data types. The subclass contains the
 * `codeMapping` dictionary which maps codes to values,
 * and a helper function `vat` that returns the value at
 * a particular index.
 */
import { CodeMapping } from "./code_mapping_interfaces";

export class DictEncoded16Array extends Int16Array {
  codeMapping: CodeMapping;

  constructor(props: ArrayBufferLike) {
    super(props);
    this.codeMapping = {};
  }

  setCodeMapping = (codeMap: CodeMapping) => {
    this.codeMapping = codeMap;
  };

  vat = (index: number): string => this.codeMapping?.[this[index]] ?? "NaN";
}
