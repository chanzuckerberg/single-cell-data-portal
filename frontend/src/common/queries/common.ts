import { useAuth0 } from "@auth0/auth0-react";

export const DEFAULT_FETCH_OPTIONS: RequestInit = {
  credentials: "include",
};

export const DELETE_FETCH_OPTIONS: RequestInit = {
  credentials: "include",
  method: "DELETE",
};

export const JSON_BODY_FETCH_OPTIONS: RequestInit = {
  headers: {
    "Content-Type": "application/json",
  },
};

export function withAuthorizationHeader(
  options: RequestInit, token: string
): RequestInit {
  const headers: any = options["headers"] || {};
  headers["Authorization"] = `Bearer ${token}`;
  return options;
}

export function useAccessToken(apiCaller: CallableFunction) {
  const { getAccessTokenSilently } = useAuth0();
  return async (...args: any) => {
    const token: string = await getAccessTokenSilently();
    console.log(`Retrieved access token: ${token}`);
    return apiCaller(...args, token);
  }
}

