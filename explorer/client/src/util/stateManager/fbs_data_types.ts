import { DictEncoded8FBArray } from "./net-encoding/dict-encoded8-f-b-array";
import { DictEncoded16FBArray } from "./net-encoding/dict-encoded16-f-b-array";
import { DictEncoded32FBArray } from "./net-encoding/dict-encoded32-f-b-array";
import { Float32FBArray } from "./net-encoding/float32-f-b-array";
import { Float64FBArray } from "./net-encoding/float64-f-b-array";
import { Int32FBArray } from "./net-encoding/int32-f-b-array";
import { JSONEncodedFBArray } from "./net-encoding/j-s-o-n-encoded-f-b-array";
import { Uint32FBArray } from "./net-encoding/uint32-f-b-array";
import { Int16EncodedXFBArray } from "./net-encoding/int16-encoded-x-f-b-array";

export const NetEncodingDict = {
  DictEncoded8FBArray,
  DictEncoded16FBArray,
  DictEncoded32FBArray,
};

export const NetEncoding = {
  Int16EncodedXFBArray,
  Float32FBArray,
  Float64FBArray,
  Int32FBArray,
  JSONEncodedFBArray,
  Uint32FBArray,
};
