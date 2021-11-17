import { Collection } from "../entities";
import { TombstonedCollection } from "../queries/collections";

export function isTombstonedCollection(
  collection: TombstonedCollection | Collection | undefined
): collection is TombstonedCollection {
  return (
    !!collection && "tombstone" in collection && collection.tombstone === true
  );
}
