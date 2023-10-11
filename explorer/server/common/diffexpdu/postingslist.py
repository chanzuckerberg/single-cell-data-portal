import itertools
import struct
from dataclasses import dataclass
from enum import IntEnum
from typing import ClassVar, List, Tuple, Type, TypeVar, Union

import bitarray  # https://github.com/ilanschnell/bitarray
import bitarray.util
import numba as nb
import numpy as np

from .utils import _blockCompress, _blockDecompress, _nonzero_bits, pairwise

"""
The format is documented in `dev_docs/diffexpdu.md`. It is implemented in both Python
and Typescript.  Implementations must maintain compatibility.
"""

PostingsList = Union[np.ndarray, List[int]]
DisjointPostingsLists = Tuple[PostingsList]
ListId = int

T = TypeVar("T")
MAGIC_NUMBER = 0xCE

headerPacker = struct.Struct("<BBH")
blockDescriptionPacker = struct.Struct("<BBHHH")
uint16 = struct.Struct("<H")


@dataclass
class _BlockDescription:
    class BlockType(IntEnum):
        BitArray = 0
        Uint16Array = 1
        Uint16ArrayInverted = 2

    block_type: BlockType
    list_id_mask: int  # 8 bit bitmask
    n_elem: int  # n elements in block
    key: int  # top two bytes (16 MSB)
    n_bytes: int  # encoded length in bytes, excluding header

    packed_length: ClassVar[int] = blockDescriptionPacker.size

    @classmethod
    def from_bytes(cls: Type[T], buf: Union[bytes, bytearray]) -> T:
        assert len(buf) >= blockDescriptionPacker.size
        block_type, list_id_mask, n_elem, key, n_bytes = blockDescriptionPacker.unpack_from(buf)
        return _BlockDescription(
            block_type=block_type, list_id_mask=list_id_mask, n_elem=n_elem + 1, key=key, n_bytes=n_bytes + 1
        )

    def to_bytes(self) -> bytes:
        return blockDescriptionPacker.pack(
            self.block_type, self.list_id_mask, self.n_elem - 1, self.key, self.n_bytes - 1
        )

    def list_ids(self) -> List[int]:
        """
        Return list IDs as list of ints
        """
        return _nonzero_bits(self.list_id_mask)


@dataclass
class _ListPartition:
    list_id: int  # 0-7, list identifier
    key: int
    start_idx: int
    end_idx: int  # exclusive range, eg, [start_idx, end_idx)
    arr: np.ndarray  # the full postings list


def packed_length(buf: Union[bytes, bytearray, memoryview], offset=0) -> int:
    """
    Return the length of the deflated posting list encoded in buf.

    Returns
    -------
    int
        Length of deflated posting list in bytes.
    """
    buf = memoryview(buf)
    if offset > 0:
        buf = buf[offset:]

    _, block_descriptions = _decode_header(buf)
    n_blocks = len(block_descriptions)
    return (
        headerPacker.size
        + n_blocks * _BlockDescription.packed_length
        + sum([block.n_bytes for block in block_descriptions])
    )


def deflate_postings_lists(disjoint_postings_lists: DisjointPostingsLists, sorted: bool = False) -> bytearray:
    """
    Efficiently pack one or more (max 8) monotonically increasing lists of unsigned integers into a byte stream.
    Lists must be disjoint (the current implementation may not fail if lists are not unique, but future
    changes may enforce/require this)
    """
    assert 0 < len(disjoint_postings_lists) <= 8

    def normalize_plist(pl):
        if isinstance(pl, list):
            pl = np.array(pl, dtype=np.uint32)
        if pl.dtype != np.uint32:
            pl = pl.astype(np.uint32)
        if not sorted:
            pl.sort()
        return pl

    disjoint_lists = tuple((normalize_plist(plist) for plist in disjoint_postings_lists))
    all_partitions = partition_lists(disjoint_lists)
    all_partitions.sort(key=lambda partition: partition.key)

    encoded_blocks = []
    for _, parts in itertools.groupby(all_partitions, key=lambda partition: partition.key):
        # encode each list as individual blocks.
        parts_tuple = tuple(parts)
        assert len(parts_tuple) <= 8
        encoded_blocks.extend(_encode_partition(part) for part in parts_tuple)

    block_count = len(encoded_blocks)
    encoded_postings = bytearray(4)
    headerPacker.pack_into(encoded_postings, 0, MAGIC_NUMBER, 0, block_count - 1)

    # encode block descriptions
    for hdr, values in encoded_blocks:
        assert hdr.n_bytes == len(values)
        encoded_postings.extend(hdr.to_bytes())

    assert len(encoded_postings) == headerPacker.size + blockDescriptionPacker.size * block_count

    # encode blocks
    for _, values in encoded_blocks:
        encoded_postings.extend(values)

    return encoded_postings


def _choose_block_type(part: _ListPartition) -> _BlockDescription.BlockType:
    """
    Given a partition, choose a block type which is likely near-optimal for
    encoding.  Heuristics live here.

    A description of the heuristics and their rationale is documented
    in the format description (`dev_docs/diffexpdu.md`)
    """
    n_elem = part.end_idx - part.start_idx
    assert 0 < n_elem <= 2**16
    assert (part.arr[part.start_idx] >> 16) == (part.arr[part.end_idx - 1] >> 16)

    interval = (part.arr[part.end_idx - 1] & 0xFFFF) + 1 - (part.arr[part.start_idx] & 0xFFFF)
    density = n_elem / interval

    minBitArrayThreshold = 2**11
    maxBitArrayThreshold = 2**16 - minBitArrayThreshold
    if (minBitArrayThreshold < n_elem < maxBitArrayThreshold) and (0.125 < density < 0.875):
        return _BlockDescription.BlockType.BitArray

    return _BlockDescription.BlockType.Uint16Array


def _encode_partition(part: _ListPartition) -> Tuple[_BlockDescription, bytearray]:
    """for a slice of the posting list, create a serialized block"""
    assert part.arr.dtype.kind == "u"
    arr = part.arr[part.start_idx : part.end_idx].astype(np.uint16)
    block_type = _choose_block_type(part)
    if block_type == _BlockDescription.BlockType.Uint16Array:
        return _encode_uint16_block(part, arr)
    if block_type == _BlockDescription.BlockType.BitArray:
        return _encode_bitarray_block(part, arr)

    raise Exception("Unknown block type prediction")


def _encode_bitarray_block(part: _ListPartition, arr: np.ndarray) -> Tuple[_BlockDescription, bytearray]:
    """
    Encode as bitarray, compressed with zlib.
    """
    bool_arr = np.zeros((2**16,), dtype=np.bool_)
    bool_arr[arr] = 1
    ba = bitarray.bitarray(endian="little")
    ba.pack(bool_arr.tobytes())

    buf = ba.tobytes()
    assert len(buf) == 2**13
    buf = _blockCompress(buf)
    key = part.arr[part.start_idx] >> 16
    header = _BlockDescription(
        block_type=_BlockDescription.BlockType.BitArray,
        list_id_mask=(1 << part.list_id),
        n_elem=len(arr),
        key=key,
        n_bytes=len(buf),
    )
    return (header, buf)


def _encode_uint16_block(part: _ListPartition, arr: np.ndarray) -> Tuple[_BlockDescription, bytearray]:
    """
    Encode lower two bytes only. Applies steps:
        * optional inversion
        * delta codes
        * byte shuffle
        * zlib.compress
    """
    key = part.arr[part.start_idx] >> 16
    n_elem = len(arr)

    if _interval_inverted_length(arr) + 2 >= n_elem:
        block_type = _BlockDescription.BlockType.Uint16Array
        arr = _byteshuffle(_delta(arr))
        buf = arr.tobytes()
        assert len(buf) == 2 * n_elem
    else:
        # invert the array if it would be shorter to encode the
        # range over which it is inverted and the elements not
        # in the original
        block_type = _BlockDescription.BlockType.Uint16ArrayInverted
        buf = bytearray(4)
        uint16.pack_into(buf, 0, int(arr[0]))  # range start
        uint16.pack_into(buf, 2, int(arr[-1]))  # range end INCLUSIVE
        arr = _byteshuffle(_delta(_interval_invert(arr)))
        buf.extend(arr.tobytes())
        assert len(buf) == 4 + len(arr)

    buf = _blockCompress(buf)
    header = _BlockDescription(
        block_type=block_type, list_id_mask=(1 << part.list_id), n_elem=n_elem, key=key, n_bytes=len(buf)
    )
    return (header, buf)


def inflate_postings_lists(buf: Union[bytes, bytearray, memoryview], offset=0) -> List[np.ndarray]:
    """
    Decode a serialized postings list, returning the original list as a ``numpy.ndarray``.
    """
    buf = memoryview(buf)
    n_lists, block_descriptions = _decode_header(buf[offset:])
    n_blocks = len(block_descriptions)
    offset += headerPacker.size + n_blocks * _BlockDescription.packed_length

    # decode each block
    sub_lists = []
    for block_description in block_descriptions:
        sub_lists.extend(_decode_block(block_description, buf[offset : offset + block_description.n_bytes]))
        offset += block_description.n_bytes

    sub_lists.sort(key=lambda sl: sl[0])  # sort by list id
    concat_lists = [
        np.concatenate([g[1] for g in grp]) for _, grp in itertools.groupby(sub_lists, key=lambda sl: sl[0])
    ]
    return concat_lists


def _decode_header(buf: memoryview) -> Tuple[int, List[_BlockDescription]]:
    magic, n_lists, n_blocks = headerPacker.unpack_from(buf)
    offset = headerPacker.size
    n_blocks += 1
    n_lists += 1
    assert magic == MAGIC_NUMBER
    assert n_blocks > 0
    assert n_lists > 0

    # decode block descriptions
    block_descriptions = [
        _BlockDescription.from_bytes(buf[b:])
        for b in range(
            offset,
            offset + n_blocks * _BlockDescription.packed_length,
            _BlockDescription.packed_length,
        )
    ]
    return n_lists, block_descriptions


def _decode_block(desc: _BlockDescription, buf: memoryview) -> List[Tuple[ListId, np.ndarray]]:
    block_type = desc.block_type
    n_elem = desc.n_elem
    buf = _blockDecompress(buf)
    list_ids = desc.list_ids()

    decoded = []  # list of (ListId, ndarray)
    if block_type == _BlockDescription.BlockType.BitArray:
        decoded.append(_decode_bitarray_block(buf, list_ids))

    elif block_type == _BlockDescription.BlockType.Uint16Array:
        decoded.append(_decode_uint16_block(buf, list_ids))

    elif block_type == _BlockDescription.BlockType.Uint16ArrayInverted:
        decoded.append(_decode_uint16_inverted_block(buf, list_ids))

    else:
        raise AssertionError("Unknown block type")

    assert len(decoded) > 0
    assert sum(len(d[1]) for d in decoded) == n_elem
    decoded = [(list_id, arr.astype(np.uint32) + (desc.key << 16)) for list_id, arr in decoded]
    return decoded


def _decode_bitarray_block(buf, list_ids):
    assert len(buf) == 8192, f"Unexpected bitarray length {len(buf)}"
    assert len(list_ids) == 1
    ba = bitarray.bitarray(endian="little")
    ba.frombytes(buf)
    arr = np.frombuffer(ba.unpack(), dtype=np.bool_).nonzero()[0]
    return (list_ids[0], arr)


def _decode_uint16_block(buf, list_ids):
    assert len(list_ids) == 1
    arr = _un_delta(_un_byteshuffle(buf, dtype=np.uint16))
    return (list_ids[0], arr)


def _decode_uint16_inverted_block(buf, list_ids):
    assert len(list_ids) == 1
    interval = (
        uint16.unpack_from(buf, 0)[0],
        uint16.unpack_from(buf, 2)[0] + 1,
    )
    arr = _interval_invert(_un_delta(_un_byteshuffle(buf[4:], dtype=np.uint16)), interval=interval)
    return (list_ids[0], arr)


def _byteshuffle(arr: np.ndarray) -> np.ndarray:
    """classic byteshuffle to improve compressibility.  Always in-place."""
    return arr.view(np.uint8).reshape((arr.itemsize, len(arr)), order="F").flatten()


def _un_byteshuffle(buf: Union[bytes, bytearray, memoryview], dtype) -> np.ndarray:
    """invert byteshuffle, in-place."""
    dtype = np.dtype(dtype)
    arr = np.frombuffer(buf, dtype=np.uint8)
    return arr.reshape(len(buf) // dtype.itemsize, dtype.itemsize, order="F").flatten().view(dtype)


def _delta(arr: np.ndarray) -> np.ndarray:
    return np.diff(arr, prepend=arr.dtype.type(0))


def _un_delta(arr: np.ndarray) -> np.ndarray:
    return np.cumsum(arr)


@nb.jit
def _interval_invert_inner(src: np.ndarray, dst: np.ndarray, interval: Tuple[int, int]) -> np.ndarray:
    src_idx = 0
    dst_idx = 0
    for i in range(interval[0], interval[1]):
        if src_idx < len(src) and src[src_idx] == i:
            src_idx += 1
            continue
        dst[dst_idx] = i
        dst_idx += 1
    return dst


def _interval_invert(arr: np.ndarray, interval=None) -> np.ndarray:
    """Invert over an interval."""
    assert arr.dtype.kind == "u"
    assert interval is not None or len(arr) > 0
    if interval is None:
        interval = (arr[0], arr[-1] + 1)
    assert len(arr) == 0 or (arr[0] >= interval[0] and arr[0] < interval[1])
    assert len(arr) == 0 or (arr[-1] >= interval[0] and arr[-1] < interval[1])
    return _interval_invert_inner(
        arr, np.empty_like(arr, shape=(_interval_inverted_length(arr, interval=interval),)), interval
    )


def _interval_inverted_length(arr: np.ndarray, interval=None) -> int:
    """
    For the given array, return the length if it was inverted across
    the interval of actual values.

    Assert: arr is an increasing monotone.
    """
    if interval is None:
        interval = (arr[0], arr[-1] + 1)
    return interval[1] - interval[0] - len(arr)


@nb.jit
def _find_partition_boundaries(arr: np.ndarray) -> List[int]:
    """
    Given an array of monotonically increasing unsigned ints, return the
    indices where the most significant 16 bits changes.

    Currently a simple scan.
    """
    partitions = []
    if len(arr) == 0:
        return partitions

    current_partition = -1
    for idx in range(len(arr)):
        val = arr[idx]
        partition = val >> 16
        if partition == current_partition:
            continue

        partitions.append(idx)
        current_partition = val >> 16

    partitions.append(len(arr))  # terminal
    return partitions


def partition_lists(plists: DisjointPostingsLists) -> List[_ListPartition]:
    assert 0 < len(plists) <= 8
    partitions = []
    for list_id in range(len(plists)):
        plist = plists[list_id]
        if isinstance(plist, list):
            plist = np.array(plist, dtype=np.uint32)
        boundaries = _find_partition_boundaries(plist)
        partitions += [
            _ListPartition(
                list_id=list_id, key=(plist[start_idx] >> 16), start_idx=start_idx, end_idx=end_idx, arr=plist
            )
            for start_idx, end_idx in pairwise(boundaries)
        ]
    return partitions
