import { useQuery } from "react-query";
import { Collection } from "src/common/entities";
import { apiTemplateToUrl } from "src/common/utils/apiTemplateToUrl";
import { API_URL } from "src/configs/configs";
import { API } from "../API";
import { ENTITIES } from "./entities";

export const USE_COLLECTIONS = {
  entities: [ENTITIES.COLLECTION],
  id: "collections",
};

export interface CollectionResponse {
  id: string;
  created_at: string;
}

async function fetchCollections(): Promise<CollectionResponse[]> {
  const json = await (await fetch(API_URL + API.COLLECTIONS)).json();

  return json.collections;
}

export function useCollections() {
  return useQuery([USE_COLLECTIONS], fetchCollections);
}

export const USE_COLLECTION = {
  entities: [ENTITIES.COLLECTION, ENTITIES.DATASET],
  id: "collection",
};

async function fetchCollection(_: unknown, id: string): Promise<Collection> {
  return (
    await fetch(apiTemplateToUrl(API_URL + API.COLLECTION, { id }))
  ).json();
}

export function useCollection(id: string) {
  return useQuery<Collection>([USE_COLLECTION, id], fetchCollection);
}
