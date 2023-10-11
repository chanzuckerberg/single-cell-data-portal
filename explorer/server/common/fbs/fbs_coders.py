import json

import numpy as np
import pandas as pd
from scipy import sparse as sp

import server.common.fbs.NetEncoding.DictEncoded8FBArray as DictEncoded8FBArray
import server.common.fbs.NetEncoding.DictEncoded16FBArray as DictEncoded16FBArray
import server.common.fbs.NetEncoding.DictEncoded32FBArray as DictEncoded32FBArray
import server.common.fbs.NetEncoding.Float32FBArray as Float32FBArray
import server.common.fbs.NetEncoding.Float64FBArray as Float64FBArray
import server.common.fbs.NetEncoding.Int16EncodedXFBArray as Int16EncodedXFBArray
import server.common.fbs.NetEncoding.Int32FBArray as Int32FBArray
import server.common.fbs.NetEncoding.JSONEncodedFBArray as JSONEncodedFBArray
import server.common.fbs.NetEncoding.TypedFBArray as TypedFBArray
import server.common.fbs.NetEncoding.Uint32FBArray as Uint32FBArray
from server.common.utils.type_conversion_utils import get_encoding_dtype_of_array


class DenseNumericIntCoder:
    n_slots = 4

    def encode_array(self, array, builder, _dtype, num_bins=None):
        if num_bins is None:
            raise ValueError("num_bins must be specified for DenseNumericIntCoder")

        # convert pandas series to numpy array
        if isinstance(array, pd.Series):
            array = array.to_numpy()
        elif sp.issparse(array):
            array = array.A.flatten()

        max_val = array.max().astype("float32")
        min_val = array.min().astype("float32")
        int_coded = np.int16((array - min_val) / (max_val - min_val) * num_bins)
        vec = builder.CreateNumpyVector(int_coded)

        builder.StartObject(self.n_slots)
        builder.PrependUOffsetTRelativeSlot(0, vec, 0)
        builder.PrependFloat32Slot(1, max_val, 0)
        builder.PrependFloat32Slot(2, min_val, 0)
        builder.PrependInt32Slot(3, num_bins, 0)
        return builder.EndObject()

    def decode_array(self, u, TarType):
        arr = TarType()
        arr.Init(u.Bytes, u.Pos)
        codes = arr.CodesAsNumpy()
        max_val = arr.Max()
        min_val = arr.Min()
        num_bins = arr.Nbins()
        return (codes / num_bins * (max_val - min_val) + min_val).astype("float32")


class DenseNumericCoder:
    n_slots = 1

    def encode_array(self, array, builder, dtype, **kwargs):
        # convert pandas series to numpy array
        if isinstance(array, pd.Series):
            array = array.to_numpy()
        elif sp.issparse(array):
            array = array.A.flatten()
        # convert to the specified dtype
        if np.dtype(array.dtype).str != dtype:
            array = array.astype(dtype)

        vec = builder.CreateNumpyVector(array)
        builder.StartObject(self.n_slots)
        builder.PrependUOffsetTRelativeSlot(0, vec, 0)
        return builder.EndObject()

    def decode_array(self, u, TarType):
        arr = TarType()
        arr.Init(u.Bytes, u.Pos)
        return arr.DataAsNumpy()


class CategoricalCoder:
    n_slots = 2

    def encode_array(self, array, builder, dtype, **kwargs):
        if isinstance(array, pd.Series) and array.dtype.name == "category":
            # create the code-to-value dictionary and encode in utf-8 as a byte array
            dictionary = np.array(bytearray(json.dumps(dict(enumerate(array.cat.categories))), "utf-8"))
            codes = array.cat.codes.values

            # ensure that the dtype is able to afford the number of categories
            assert len(array.cat.categories) <= np.iinfo(dtype).max + 1

            vec_codes = builder.CreateNumpyVector(codes)
            vec_dictionary = builder.CreateNumpyVector(dictionary)
            builder.StartObject(self.n_slots)
            builder.PrependUOffsetTRelativeSlot(0, vec_codes, 0)
            builder.PrependUOffsetTRelativeSlot(1, vec_dictionary, 0)
            return builder.EndObject()
        else:
            raise ValueError("Input array must be pandas Categorical.")

    def decode_array(self, u, TarType):
        # returns a pandas Categorical
        arr = TarType()
        arr.Init(u.Bytes, u.Pos)
        codes = arr.CodesAsNumpy()
        dictionary = arr.DictAsNumpy()
        dictionary = json.loads(dictionary.tobytes().decode("utf-8"))
        return pd.Categorical.from_codes(codes, categories=list(dictionary.values()))


class PolymorphicCoder:
    n_slots = 1

    def encode_array(self, array, builder, dtype=None, **kwargs):
        # dtype is unused here as array is just getting slammed into a JSON
        if sp.issparse(array):
            array = array.A.flatten()
        array = pd.Series(array)
        as_json = array.to_json(orient="records")
        json_array = np.array(bytearray(as_json, "utf-8"))

        vec = builder.CreateNumpyVector(json_array)
        builder.StartObject(self.n_slots)
        builder.PrependUOffsetTRelativeSlot(0, vec, 0)
        return builder.EndObject()

    def decode_array(self, u, TarType):
        arr = TarType()
        arr.Init(u.Bytes, u.Pos)
        narr = arr.DataAsNumpy()
        return np.array(json.loads(narr.tobytes().decode("utf-8")))


# two-layer mapper (1) array_class, (2) encoding_dtype
# this is necessary as an int32 codes array will require a different coder
# than an int32 numeric array.
ARRAY_ENCODER = {
    "dense": {
        np.dtype(np.float64).str: (DenseNumericCoder, TypedFBArray.TypedFBArray.Float32FBArray),
        np.dtype(np.float32).str: (DenseNumericCoder, TypedFBArray.TypedFBArray.Float32FBArray),
        np.dtype(np.float16).str: (DenseNumericCoder, TypedFBArray.TypedFBArray.Float32FBArray),
        np.dtype(np.int8).str: (DenseNumericCoder, TypedFBArray.TypedFBArray.Int32FBArray),
        np.dtype(np.int16).str: (DenseNumericCoder, TypedFBArray.TypedFBArray.Int32FBArray),
        np.dtype(np.int32).str: (DenseNumericCoder, TypedFBArray.TypedFBArray.Int32FBArray),
        np.dtype(np.int64).str: (DenseNumericCoder, TypedFBArray.TypedFBArray.Int32FBArray),
        np.dtype(np.uint8).str: (DenseNumericCoder, TypedFBArray.TypedFBArray.Uint32FBArray),
        np.dtype(np.uint16).str: (DenseNumericCoder, TypedFBArray.TypedFBArray.Uint32FBArray),
        np.dtype(np.uint32).str: (DenseNumericCoder, TypedFBArray.TypedFBArray.Uint32FBArray),
        np.dtype(np.uint64).str: (DenseNumericCoder, TypedFBArray.TypedFBArray.Uint32FBArray),
    },
    "lossy": {
        np.dtype(np.float64).str: (DenseNumericIntCoder, TypedFBArray.TypedFBArray.Int16EncodedXFBArray),
        np.dtype(np.float32).str: (DenseNumericIntCoder, TypedFBArray.TypedFBArray.Int16EncodedXFBArray),
        np.dtype(np.float16).str: (DenseNumericIntCoder, TypedFBArray.TypedFBArray.Int16EncodedXFBArray),
    },
    "category": {
        np.dtype(np.int8).str: (CategoricalCoder, TypedFBArray.TypedFBArray.DictEncoded8FBArray),
        np.dtype(np.int16).str: (CategoricalCoder, TypedFBArray.TypedFBArray.DictEncoded16FBArray),
        np.dtype(np.int32).str: (CategoricalCoder, TypedFBArray.TypedFBArray.DictEncoded32FBArray),
    },
}

TYPE_MAP = {
    TypedFBArray.TypedFBArray.NONE: (None, None),
    TypedFBArray.TypedFBArray.Uint32FBArray: (DenseNumericCoder, Uint32FBArray.Uint32FBArray),
    TypedFBArray.TypedFBArray.Int32FBArray: (DenseNumericCoder, Int32FBArray.Int32FBArray),
    TypedFBArray.TypedFBArray.Float32FBArray: (DenseNumericCoder, Float32FBArray.Float32FBArray),
    TypedFBArray.TypedFBArray.Float64FBArray: (DenseNumericCoder, Float64FBArray.Float64FBArray),
    TypedFBArray.TypedFBArray.Int16EncodedXFBArray: (DenseNumericIntCoder, Int16EncodedXFBArray.Int16EncodedXFBArray),
    TypedFBArray.TypedFBArray.JSONEncodedFBArray: (PolymorphicCoder, JSONEncodedFBArray.JSONEncodedFBArray),
    TypedFBArray.TypedFBArray.DictEncoded8FBArray: (CategoricalCoder, DictEncoded8FBArray.DictEncoded8FBArray),
    TypedFBArray.TypedFBArray.DictEncoded16FBArray: (CategoricalCoder, DictEncoded16FBArray.DictEncoded16FBArray),
    TypedFBArray.TypedFBArray.DictEncoded32FBArray: (CategoricalCoder, DictEncoded32FBArray.DictEncoded32FBArray),
}


def _get_array_class(array, num_bins=None):
    # returns
    # (1) the generic type of array - category, dense, or lossy
    # this determines which coder should be used.
    # (2) the encoding data type from the corresponding attribute of the
    # source array. e.g. for pandas Categorical, the encoding data type
    # is computed from the codes array.
    if pd.api.types.is_categorical_dtype(array):
        return "category", np.dtype(array.cat.codes.values.dtype).str
    else:
        dtype = np.dtype(get_encoding_dtype_of_array(array)).str
        array_class = "lossy" if dtype.startswith("<f") and num_bins is not None else "dense"
        return array_class, dtype


def serialize_typed_array(builder, source_array, num_bins=None):
    if isinstance(source_array, pd.Index):
        source_array = source_array.to_series()

    array_class, encoding_dtype = _get_array_class(source_array, num_bins=num_bins)
    # the default coder will assume the data is polymorphic and yield a JSON encoded array
    defaultCoder = (PolymorphicCoder, TypedFBArray.TypedFBArray.JSONEncodedFBArray)
    Coder, array_type = ARRAY_ENCODER[array_class].get(encoding_dtype, defaultCoder)

    coder_obj = Coder()
    # for encoding, we require the source array, flatbuffer builder, encoding data type, and
    # number of bins for lossy integer compression
    array_value = coder_obj.encode_array(source_array, builder, encoding_dtype, num_bins=num_bins)
    return (array_type, array_value)


def deserialize_typed_array(tarr):
    (union_type, u) = tarr
    if union_type is TypedFBArray.TypedFBArray.NONE:
        return None

    Coder, TarType = TYPE_MAP.get(union_type, None)
    if TarType is None:
        raise TypeError(f"FBS contains unknown data type: {union_type}")
    # for decoding, we require the encoded column and the type of typed array
    return Coder().decode_array(u, TarType)
