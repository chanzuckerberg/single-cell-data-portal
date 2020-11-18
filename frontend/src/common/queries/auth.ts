import { useQuery } from "react-query";
import { API_URL } from "src/configs/configs";
import { API } from "../API";
import { DEFAULT_FETCH_OPTIONS } from "./common";
import { ENTITIES } from "./entities";

export const USE_USER_INFO = {
  entities: [ENTITIES.USER_INFO],
  id: "userInfo",
};

export interface UserInfoResponse {
  email?: string;
  email_verified?: boolean;
  id?: string;
  is_authenticated?: boolean;
  name?: string;
}

async function fetchUserInfo(
  _: unknown,
  { hasAuth }: { hasAuth: boolean }
): Promise<UserInfoResponse | null> {
  if (!hasAuth) return Promise.resolve(null);

  const response = await fetch(API_URL + API.USER_INFO, DEFAULT_FETCH_OPTIONS);

  const result = await response.json();

  if (!response.ok) {
    throw result;
  }

  return result;
}

export function useUserInfo(hasAuth: boolean) {
  return useQuery([USE_USER_INFO, { hasAuth }], fetchUserInfo);
}
