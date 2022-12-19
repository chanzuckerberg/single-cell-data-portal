import { useQuery, UseQueryResult } from "react-query";
import { API_URL } from "src/configs/configs";
import { API } from "../API";
import { DEFAULT_FETCH_OPTIONS } from "./common";
import { ENTITIES } from "./entities";

export interface APIKeyResponse {
  id: string;
  key?: string;
}

export const USE_CURATOR_AUTH_KEY = {
  entities: [ENTITIES.API_KEY],
  id: "apiKey",
};

export async function generateCuratorAuthKey(): Promise<APIKeyResponse | null> {
  const response = await fetch(API_URL + API.CURATOR_AUTH_KEY, {
    ...DEFAULT_FETCH_OPTIONS,
    method: "POST",
  });

  const result = await response.json();

  if (!response.ok) {
    throw result;
  }

  return result;
}

async function fetchAuthKeyID(): Promise<APIKeyResponse | null> {
  const response = await fetch(
    API_URL + API.CURATOR_AUTH_KEY,
    DEFAULT_FETCH_OPTIONS
  );

  const result = await response.json();

  if (!response.ok) {
    throw result;
  }

  return result;
}

export function useCuratorAuthKeyID(): UseQueryResult<APIKeyResponse | null> {
  return useQuery([USE_CURATOR_AUTH_KEY], fetchAuthKeyID);
}
