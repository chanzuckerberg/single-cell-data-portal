import tiledb

filters_categorical = [tiledb.DictionaryFilter(), tiledb.ZstdFilter(level=+19)]
filters_numeric = [tiledb.ByteShuffleFilter(), tiledb.ZstdFilter(level=+5)]
