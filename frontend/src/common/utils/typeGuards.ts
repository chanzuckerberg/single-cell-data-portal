import { TombstonedCollection } from "../queries/collections";

export function isTombstonedCollection(
  collection: any
): collection is TombstonedCollection {
  return (
    collection && "tombstone" in collection && collection.tombstone === true
  );
}
