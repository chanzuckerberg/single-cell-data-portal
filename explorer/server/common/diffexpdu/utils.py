"""
Utility/helper code for diffexpdu package
"""

import itertools
import zlib
from typing import List, Union


def _nonzero_bits(mask: int) -> List[int]:
    ids = []
    for b in range(8):
        if mask & (1 << b):
            ids.append(b)
    return ids


def pairwise(iterable):
    # Unfortunately, only available in Python >= 3.10
    # pairwise('ABCDEFG') --> AB BC CD DE EF FG
    a, b = itertools.tee(iterable)
    next(b, None)
    return zip(a, b)


"""
Compression helpers

The default zlib.compress/decompress functions include a header and checksum, for a
total of (at least) 6 bytes per compressed object. For the bulk compression methods
(JSON, uint32), this is negligible. For the incrementally compressed blocks it can
add quite a bit of overhead if there are a large number of blocks.

These two routines are optimized to compress/decompress the block contents:
* no header or checksum (assumed to be handled at a higher level of the stack)
* configuration changes to improve compression of the small integer values
"""


def _blockCompress(src: Union[bytes, memoryview], strategy=zlib.Z_DEFAULT_STRATEGY) -> bytes:
    compressor = zlib.compressobj(level=3, wbits=-15, memLevel=9, strategy=strategy)
    dst = compressor.compress(src)
    dst += compressor.flush(zlib.Z_FINISH)
    return dst


def _blockDecompress(src: memoryview) -> bytes:
    decompressor = zlib.decompressobj(wbits=-15)
    dst = decompressor.decompress(src)
    dst += decompressor.flush()
    return dst
