import { useQuery, UseQueryResult } from "react-query";
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
  // const response = await fetch(API_URL + API.CURATOR_AUTH_KEY, {
  //   ...DEFAULT_FETCH_OPTIONS,
  //   method: "PUT",
  // });

  // DEBUG
  // DEBUG
  // DEBUG
  // DEBUG
  const response = {
    json: () => {
      return {
        id: "1234567890123456789012345678901234567890",
        key: "1234567890123456789012345678901234567890",
      };
    },
    ok: true,
  };

  const result = await response.json();

  if (!response.ok) {
    throw result;
  }

  return result;
}

async function fetchAuthKeyID(): Promise<APIKeyResponse | null> {
  // const response = await fetch(
  //   API_URL + API.CURATOR_AUTH_KEY,
  //   DEFAULT_FETCH_OPTIONS
  // );

  // DEBUG
  // DEBUG
  // DEBUG
  // DEBUG
  const response = {
    json: () => {
      return {
        id: "1234567890123456789012345678901234567890",
      };
    },
    ok: true,
  };

  const result = await response.json();

  if (!response.ok) {
    throw result;
  }

  return result;
}

export function useCuratorAuthKeyID(): UseQueryResult<APIKeyResponse | null> {
  return useQuery([USE_CURATOR_AUTH_KEY], fetchAuthKeyID);
}
