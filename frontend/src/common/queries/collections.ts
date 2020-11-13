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

export async function createCollection(): Promise<string> {
  // DEBUG
  // DEBUG
  // DEBUG
  await new Promise((resolve) => {
    setTimeout(() => {
      resolve();
    }, 2 * 1000);
  });

  return "673637cf-dcb7-45e1-bb88-72a27c50c8ca";
}

export const formDataToObject = function (formData: FormData) {
  const payload: { [key: string]: unknown } = {};

  formData.forEach((value, key: string) => {
    const translatedKey = key.replace("-", "_");

    payload[translatedKey] = value;
  });

  return payload;
};
