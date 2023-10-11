import struct
from dataclasses import dataclass
from enum import IntEnum
from typing import ClassVar, Type, TypeVar, Union

import numpy as np

from .postingslist import deflate_postings_lists, inflate_postings_lists

__all__ = ["DiffExArguments"]

# Pack/unpack support for various types
# The format is documented in `dev_docs/diffexpdu.md`
headerPacker = struct.Struct("<BB")
topNParamsPacker = struct.Struct("<H")

T = TypeVar("T")
MAGIC_NUMBER = 0xDE


@dataclass
class DiffExArguments:
    """
    Encode arguments for the differential expression (/diffex) REST API in a binary format. PDU format is documented
    in `dev_docs/diffexpdu.md`, and is shared by the Typescript implementation.

    IMPORTANT: the two cell sets *MUST* be disjoint and sorted (increasing).

    Examples
    --------
    ```
        # specify arguments
        de_args = DiffExArguments(
            mode=DiffExArguments.DiffExMode.TopN,
            params=DiffExArguments.TopNParams(N=50),
            set1=[0, 1, 2],
            set2=[3, 88, 1010],
        )
        # encode into a buffer
        encoded_postings = de_args.pack()

        # decode buffer to arguments
        de = DiffExArguments.unpack_from(buf)
        assert de.mode == DiffExArguments.DiffExMode.TopN
    ```
    """

    class DiffExMode(IntEnum):
        TopN = 0
        VarFilter = 1  # unsupported

    @dataclass
    class TopNParams:
        N: int = 50  # parameter for TopN mode
        packed_size: ClassVar[int] = topNParamsPacker.size

        @classmethod
        def unpack_from(cls: Type[T], buf: bytes, offset=0) -> T:
            (N,) = topNParamsPacker.unpack_from(buf, offset)
            return DiffExArguments.TopNParams(N=N)

        def pack(self) -> bytes:
            return topNParamsPacker.pack(self.N)

    mode: DiffExMode
    params: TopNParams
    set1: np.ndarray
    set2: np.ndarray

    @classmethod
    def unpack_from(cls: Type[T], buf: Union[bytes, bytearray, memoryview], offset=0) -> T:
        """
        Given a buffer containing encoded parameters, unpack and return an instance of DiffExArguments.
        """
        (magic, mode) = headerPacker.unpack_from(buf, offset)
        assert magic == MAGIC_NUMBER
        assert mode == DiffExArguments.DiffExMode.TopN
        offset += headerPacker.size

        params = DiffExArguments.TopNParams.unpack_from(buf, offset)
        offset += DiffExArguments.TopNParams.packed_size
        set1, set2 = inflate_postings_lists(buf, offset=offset)

        return DiffExArguments(mode=mode, params=params, set1=set1, set2=set2)

    def pack(self):
        """Pack the instance of DiffExArguments into a buffer."""
        assert self.mode == DiffExArguments.DiffExMode.TopN
        return (
            headerPacker.pack(MAGIC_NUMBER, self.mode)
            + self.params.pack()
            + deflate_postings_lists((self.set1, self.set2))
        )

    def __eq__(self, other) -> bool:
        return (
            self.mode == other.mode
            and self.params == other.params
            and np.array_equal(self.set1, other.set1)
            and np.array_equal(self.set2, other.set2)
        )
