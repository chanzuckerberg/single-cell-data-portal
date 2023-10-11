import { flatbuffers } from "flatbuffers";
import { NetEncoding, NetEncodingDict } from "./fbs_data_types";
import {
  TypedArray,
  isDictEncodedTypedArray,
  DictEncodedArray,
  isFloatTypedArray,
} from "../../common/types/arraytypes";
import {
  Dataframe,
  IdentityInt32Index,
  DenseInt32Index,
  KeyIndex,
  DataframeValueArray,
} from "../dataframe";

import { DictEncoded8Array } from "./dict-encoded8_array";
import { DictEncoded16Array } from "./dict-encoded16_array";
import { DictEncoded32Array } from "./dict-encoded32_array";
import { TypedFBArray } from "./net-encoding/typed-f-b-array";
import { Matrix } from "./net-encoding/matrix";
import { Column } from "./net-encoding/column";

const utf8Decoder = new TextDecoder("utf-8");
const DEFAULT_NUM_BINS = 5000;

/**
 * Matrix flatbuffer decoding support. See fbs/matrix.fbs
 */

/**
 * Decode TypedFBArray
 */

function decodeDictArray(
  uType: TypedFBArray,
  uValF: Column["u"],
): DictEncodedArray {
  const TypeClass =
    NetEncodingDict[TypedFBArray[uType] as keyof typeof NetEncodingDict];
  const arr = uValF(new TypeClass());
  const codesArray = arr.codesArray();
  let codesToValues = arr.dictArray();
  codesToValues = JSON.parse(utf8Decoder.decode(codesToValues));

  let data: DictEncodedArray;
  if (uType === TypedFBArray.DictEncoded8FBArray) {
    data = new DictEncoded8Array(codesArray);
    data.setCodeMapping(codesToValues);
  } else if (uType === TypedFBArray.DictEncoded16FBArray) {
    data = new DictEncoded16Array(codesArray);
    data.setCodeMapping(codesToValues);
  } else {
    data = new DictEncoded32Array(codesArray);
    data.setCodeMapping(codesToValues);
  }
  return data;
}

function decodeNumericArray(
  uType: TypedFBArray,
  uValF: Column["u"],
  inplace = false,
): TypedArray {
  const TypeClass =
    NetEncoding[TypedFBArray[uType] as keyof typeof NetEncoding];
  const arr = uValF(new TypeClass());
  let dataArray = arr.dataArray();
  if (!inplace) {
    /* force copy to release underlying FBS buffer */
    dataArray = new dataArray.constructor(dataArray);
  }
  return dataArray;
}

function decodeIntCodedArray(
  uType: TypedFBArray,
  uValF: Column["u"],
): TypedArray {
  const TypeClass =
    NetEncoding[TypedFBArray[uType] as keyof typeof NetEncoding];
  const arr = uValF(new TypeClass());
  const codesArray = arr.codesArray();
  const dataArray = new Float32Array(codesArray.length);
  const maxValue = arr.max();
  const minValue = arr.min();
  const numBins = arr.nbins();
  for (let i = 0; i < dataArray.length; i += 1) {
    dataArray[i] = (codesArray[i] / numBins) * (maxValue - minValue) + minValue;
  }
  return dataArray;
}

function decodeJSONArray(
  uType: TypedFBArray,
  uValF: Column["u"],
): DataframeValueArray {
  const TypeClass =
    NetEncoding[TypedFBArray[uType] as keyof typeof NetEncoding];
  const arr = uValF(new TypeClass());
  const dataArray = arr.dataArray();
  const json = utf8Decoder.decode(dataArray);
  return JSON.parse(json);
}

function decodeTypedArray(
  uType: TypedFBArray,
  uValF: Column["u"],
  inplace = false,
): TypedArray | DataframeValueArray | DictEncodedArray | null {
  if (uType === TypedFBArray.NONE) {
    return null;
  }
  switch (uType) {
    case TypedFBArray.JSONEncodedFBArray: {
      return decodeJSONArray(uType, uValF);
    }
    case TypedFBArray.DictEncoded8FBArray:
    case TypedFBArray.DictEncoded16FBArray:
    case TypedFBArray.DictEncoded32FBArray: {
      return decodeDictArray(uType, uValF);
    }
    case TypedFBArray.Int16EncodedXFBArray: {
      return decodeIntCodedArray(uType, uValF);
    }
    default: {
      return decodeNumericArray(uType, uValF, inplace);
    }
  }
}

/**
 * Parameter: Uint8Array or ArrayBuffer containing raw flatbuffer Matrix
 * Returns: object containing decoded Matrix:
 * {
 *   nRows: num,
 *   nCols: num,
 *   columns: [
 *     each column, which will be a TypedArray or Array
 *   ],
 *   colIdx: [] | null,
 * }
 */
export function decodeMatrixFBS(arrayBuffer: any, inplace = false) {
  const bb = new flatbuffers.ByteBuffer(new Uint8Array(arrayBuffer));
  const matrix = Matrix.getRootAsMatrix(bb);

  const nRows = matrix.nRows();
  const nCols = matrix.nCols();

  /* decode columns */
  const columnsLength = matrix.columnsLength();
  const columns = Array(columnsLength).fill(null);
  for (let c = 0; c < columnsLength; c += 1) {
    const col = matrix.columns(c);
    if (!col) {
      columns[c] = null;
    } else {
      const arr = decodeTypedArray(col.uType(), col.u.bind(col), inplace);
      columns[c] = arr;
    }
  }

  /* decode col_idx */
  const colIdx = decodeTypedArray(
    matrix.colIndexType(),
    matrix.colIndex.bind(matrix),
    inplace,
  );
  return {
    nRows,
    nCols,
    columns,
    colIdx,
    rowIdx: null,
  };
}

function encodeDictArray(
  builder: flatbuffers.Builder,
  uType: TypedFBArray,
  uData: DictEncodedArray,
): any {
  const [dCodes, dDict] = _getCodesAndDictVectorOffset(builder, uType, uData);
  builder.startObject(2);
  builder.addFieldOffset(0, dCodes, 0);
  builder.addFieldOffset(0, dDict, 0);
  return builder.endObject();
}

function _getCodesAndDictVectorOffset(
  builder: flatbuffers.Builder,
  uType: TypedFBArray,
  uData: DictEncodedArray,
): [number, number] {
  const uTypeName = TypedFBArray[uType];
  const ArrayType = NetEncodingDict[uTypeName as keyof typeof NetEncodingDict];

  const json = JSON.stringify(uData.codeMapping);
  const utf8Encoder = new TextEncoder();
  const jsonUTF8 = utf8Encoder.encode(json);

  let dCodes;
  if (ArrayType === NetEncodingDict.DictEncoded8FBArray) {
    dCodes = NetEncodingDict.DictEncoded8FBArray.createCodesVector(
      builder,
      uData as DictEncoded8Array,
    );
  } else if (ArrayType === NetEncodingDict.DictEncoded16FBArray) {
    dCodes = NetEncodingDict.DictEncoded16FBArray.createCodesVector(
      builder,
      uData as DictEncoded16Array,
    );
  } else if (ArrayType === NetEncodingDict.DictEncoded32FBArray) {
    dCodes = NetEncodingDict.DictEncoded32FBArray.createCodesVector(
      builder,
      uData as DictEncoded32Array,
    );
  } else {
    throw new Error(`unsupported dictionary-encoded array type ${uTypeName}`);
  }

  const dDict = ArrayType.createDictVector(builder, jsonUTF8);
  return [dCodes, dDict];
}

function _getNumericDataVectorOffset(
  builder: flatbuffers.Builder,
  uType: TypedFBArray,
  uData: TypedArray,
): number {
  const uTypeName = TypedFBArray[uType];
  const ArrayType = NetEncoding[uTypeName as keyof typeof NetEncoding];
  let dArray;
  if (ArrayType === NetEncoding.Int32FBArray) {
    dArray = NetEncoding.Int32FBArray.createDataVector(
      builder,
      uData as Int32Array,
    );
  } else if (ArrayType === NetEncoding.Float32FBArray) {
    dArray = NetEncoding.Float32FBArray.createDataVector(
      builder,
      uData as Float32Array,
    );
  } else if (ArrayType === NetEncoding.Float64FBArray) {
    dArray = NetEncoding.Float64FBArray.createDataVector(
      builder,
      uData as Float64Array,
    );
  } else if (ArrayType === NetEncoding.Uint32FBArray) {
    dArray = NetEncoding.Uint32FBArray.createDataVector(
      builder,
      uData as Uint32Array,
    );
  } else {
    throw new Error(`unsupported numeric array type ${uTypeName}`);
  }
  return dArray;
}

function encodeNumericArray(
  builder: flatbuffers.Builder,
  uType: TypedFBArray,
  uData: TypedArray,
): number {
  const dArray = _getNumericDataVectorOffset(builder, uType, uData);
  builder.startObject(1);
  builder.addFieldOffset(0, dArray, 0);
  return builder.endObject();
}

function encodeIntCodedArray(
  builder: flatbuffers.Builder,
  uData: TypedArray,
): number {
  const codesArray = new Int16Array(uData.length);
  const maxValue = Math.max(...codesArray);
  const minValue = Math.min(...codesArray);
  for (let i = 0; i < codesArray.length; i += 1) {
    codesArray[i] = Math.floor(
      ((uData[i] - minValue) / (maxValue - minValue)) * DEFAULT_NUM_BINS,
    );
  }
  const cArray = NetEncoding.Int16EncodedXFBArray.createCodesVector(
    builder,
    codesArray,
  );
  builder.startObject(4);
  builder.addFieldOffset(0, cArray, 0);
  NetEncoding.Int16EncodedXFBArray.addMax(builder, maxValue);
  NetEncoding.Int16EncodedXFBArray.addMin(builder, minValue);
  NetEncoding.Int16EncodedXFBArray.addNbins(builder, DEFAULT_NUM_BINS);
  return builder.endObject();
}

function encodeJSONArray(
  builder: flatbuffers.Builder,
  uData: Array<string> | Array<boolean>,
): any {
  const json = JSON.stringify(uData);
  const utf8Encoder = new TextEncoder();
  const jsonUTF8 = utf8Encoder.encode(json);
  const dArray = NetEncoding.JSONEncodedFBArray.createDataVector(
    builder,
    jsonUTF8,
  );
  builder.startObject(1);
  builder.addFieldOffset(0, dArray, 0);
  return builder.endObject();
}

function encodeTypedArray(
  builder: flatbuffers.Builder,
  uType: any,
  uData: any,
) {
  switch (uType) {
    case TypedFBArray.JSONEncodedFBArray: {
      return encodeJSONArray(builder, uData);
    }
    case TypedFBArray.DictEncoded8FBArray:
    case TypedFBArray.DictEncoded16FBArray:
    case TypedFBArray.DictEncoded32FBArray: {
      return encodeDictArray(builder, uType, uData);
    }
    case TypedFBArray.Int16EncodedXFBArray: {
      return encodeIntCodedArray(builder, uData);
    }
    default: {
      return encodeNumericArray(builder, uType, uData);
    }
  }
}

/**
 * Encode the dataframe as an FBS Matrix
 */

export function encodeMatrixFBS(
  df: Dataframe,
  encodeSparse = false,
): Uint8Array {
  /* row indexing not supported currently */
  if (!(df.rowIndex instanceof IdentityInt32Index)) {
    throw new Error("FBS does not support row index encoding at this time");
  }

  const shape = df.dims;
  const utf8Encoder = new TextEncoder();
  const builder = new flatbuffers.Builder(1024);

  let encColIndex;
  let encColIndexUType;
  let encColumns;

  if (shape[0] > 0 && shape[1] > 0) {
    const columns = df.columns().map((col) => col.asArray());

    const cols = columns.map((carr) => {
      let { name } = carr.constructor;
      if (encodeSparse) {
        name = `Sparse${name}`;
      }
      name = name.replace("Array", "FBArray");
      const uType = TypedFBArray[name as keyof typeof TypedFBArray];
      const tarr = encodeTypedArray(builder, uType, carr);
      Column.startColumn(builder);
      Column.addUType(builder, uType);
      Column.addU(builder, tarr);
      return Column.endColumn(builder);
    });

    encColumns = Matrix.createColumnsVector(builder, cols);

    if (df.colIndex && shape[1] > 0) {
      const colIndexType = df.colIndex.constructor;
      if (colIndexType === IdentityInt32Index) {
        encColIndex = undefined;
      } else if (colIndexType === DenseInt32Index) {
        encColIndexUType = TypedFBArray.Int32FBArray;
        encColIndex = encodeTypedArray(
          builder,
          encColIndexUType,
          df.colIndex.labels(),
        );
      } else if (colIndexType === KeyIndex) {
        encColIndexUType = TypedFBArray.JSONEncodedFBArray;
        encColIndex = encodeTypedArray(
          builder,
          encColIndexUType,
          utf8Encoder.encode(JSON.stringify(df.colIndex.labels())),
        );
      } else {
        throw new Error("Index type FBS encoding unsupported");
      }
    }
  }

  Matrix.startMatrix(builder);
  Matrix.addNRows(builder, shape[0]);
  Matrix.addNCols(builder, shape[1]);
  if (encColumns) {
    Matrix.addColumns(builder, encColumns);
  }
  if (encColIndexUType) {
    Matrix.addColIndexType(builder, encColIndexUType);
    Matrix.addColIndex(builder, encColIndex);
  }
  const root = Matrix.endMatrix(builder);
  builder.finish(root);
  return builder.asUint8Array();
}

/**
 * Decide what internal data type to use for the data returned from
 * the server.
 *
 * TODO - future optimization: not all int32/uint32 data series require
 * promotion to float64.  We COULD simply look at the data to decide.
 */
function promoteTypedArray(o: TypedArray) {
  if (isFloatTypedArray(o) || Array.isArray(o)) return o;

  let TypedArrayCtor;
  switch (o.constructor) {
    case Int8Array:
    case Uint8Array:
    case Uint8ClampedArray:
    case Int16Array:
    case Uint16Array:
      TypedArrayCtor = Float32Array;
      break;

    case Int32Array:
    case Uint32Array:
      TypedArrayCtor = Float64Array;
      break;

    default:
      throw new Error("Unexpected data type returned from server.");
  }
  if (o.constructor === TypedArrayCtor) return o;
  return new TypedArrayCtor(o);
}

/**
 * Convert array of Matrix FBS to a Dataframe.
 *
 * The application has strong assumptions that all scalar data will be
 * stored as a float32 or float64 (regardless of underlying data types).
 * For example, clipping of value ranges (eg, user-selected percentiles)
 * depends on the ability to use NaN in any numeric type.
 *
 * All float data from the server is left as is.  All non-float is promoted
 * to an appropriate float.
 */
export function matrixFBSToDataframe(
  arrayBuffers: ArrayBuffer | ArrayBuffer[],
): Dataframe {
  if (!Array.isArray(arrayBuffers)) {
    arrayBuffers = [arrayBuffers];
  }
  if (arrayBuffers.length === 0) {
    return Dataframe.empty();
  }
  const fbs = arrayBuffers.map((ab) => decodeMatrixFBS(ab, true)); // leave in place

  /* check that all FBS have same row dimensionality */
  const { nRows } = fbs[0];
  fbs.forEach((b) => {
    if (b.nRows !== nRows)
      throw new Error("FBS with inconsistent dimensionality");
  });
  const columns = fbs
    .map((fb) =>
      fb.columns.map((c) => {
        if (
          isFloatTypedArray(c) ||
          isDictEncodedTypedArray(c) ||
          Array.isArray(c)
        )
          return c;
        return promoteTypedArray(c);
      }),
    )
    .flat();

  // colIdx may be TypedArray or Array
  const colIdx = fbs
    .map((b) => {
      if (!b.colIdx) return Array.from(Array(b.nCols).keys());
      return Array.isArray(b.colIdx) ? b.colIdx : Array.from(b.colIdx);
    })
    .flat() as (string | number)[];
  const nCols = columns.length;
  const df = new Dataframe([nRows, nCols], columns, null, new KeyIndex(colIdx));
  return df;
}
