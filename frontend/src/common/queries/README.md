# React Queries

## Introduction

[`react-query`](https://github.com/tannerlinsley/react-query) is used for
Data Portal server state management, instead of using Redux Toolkit to do both
client state and server state management.

The benefits of this arrangement include:

1. Instead of replicating server schemas and entity tables, and trying to
   normalize entities and manage entity updates in the FE, we rely on straight
   forward query cache invalidations to signal `react-query` to atomically
   refetch server data when needed, since BE is and **should** be the single
   source of truth for server state, instead of having FE maintaining and
   managing a snapshot of server state, which is troublesome and inefficient

1. React Query comes with really awesome cache management, garbage collection,
   and convenient auto refetch [default configuration](https://react-query.tanstack.com/docs/guides/important-defaults).
   This means unused server state will be cleared from the browser memory after
   a certain time and refetch will happen in the background for us, such as
   when the window refocuses, internet reconnects, etc..

1. Most importantly, React Query provides a really nice developer experience
   (also has [React Query Devtools](https://github.com/tannerlinsley/react-query-devtools))
   and reduces even more code than what Redux Toolkit has already achieved
   compared to vanilla Redux

1. Note that React Query is **NOT** a replacement for client state management,
   which can be still done via Redux Toolkit and/or plain React Context/Hooks.

## File Structure

1. This folder `/queries` contains everything related to React Query and its
   helper functions, constants, etc.. Each entity "slice" has its own file,
   such as `collections.ts`

1. `entities.ts` is an enum of all entities provided from the server. Such as
   `collection`, `dataset`, etc.. And we will associate entities returned in a
   specific query with its `useHook`, in order to correctly invalidate stale
   queries programmatically.

   For example, `useCollection` in `/collections.ts` uses `USE_COLLECTION` object as its unique query key, which includes properties: `entities: string[]` and `id: string`.

   ```ts
    export const USE_COLLECTION = {
      entities: [ENTITIES.COLLECTION, ENTITIES.DATASET],
      id: "collection",
    };

    export function useCollection(id: string) {
      return useQuery<Collection>(
        [USE_COLLECTION, id],
        fetchCollection
      );
   }
   ```

   With `entities` key in every query key, we will be able to invalidate related queries' caches when we believe their data are now stale.

   One example is when we create a new collection, we want all queries that
   return collections to be invalidated. So we can check use the following
   code to achieve that:

   ```ts
      import { useMutation, useQueryCache } from 'react-query'

      const queryCache = useQueryCache();

      const USE_CREATE_COLLECTION = {
        entities: [ENTITIES.COLLECTION]
      };

      function createCollection() {
        // POST http call
      }

      // When this mutation succeeds, invalidate any queries with
      // `ENTITIES.COLLECTION` in its query key's `entities` array
      const [mutate] = useMutation(createCollection, {
        onSuccess: () => {
          queryCache.invalidateQueries((query) => {
            const key = query.queryKey[0];

            return hasIntersection(key.entities, USE_CREATE_COLLECTION.entities)
          })
        },
      })

      function hasIntersection(array1, array2) {
        const set2 = new Set(array2);

        const result = array1.filter(entity => set2.has(entity))

        return Boolean(result.length)
      }
   ```

## Resources

1. [React Query release tech talk](https://www.youtube.com/watch?v=seU46c6Jz7E)
1. [React Query docs](https://react-query.tanstack.com/docs/overview)
1. [My notes on React Query](https://docs.google.com/document/d/1G-pxb2NoUImo-aa9fWmz91K8jaeQTHNT_0CDXyqSeeE/edit)
1. [Invalidation from Mutations](https://react-query.tanstack.com/docs/guides/invalidations-from-mutations)
