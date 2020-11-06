import { useQuery } from "react-query";
import { API_URL } from "src/configs/configs";
import { API } from "../API";
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

  return (await fetch(API_URL + API.USER_INFO)).json();
}

export function useUserInfo(hasAuth: boolean) {
  return useQuery([USE_USER_INFO, { hasAuth }], fetchUserInfo);
}
