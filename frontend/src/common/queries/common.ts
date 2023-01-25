import { useAuth0, LocalStorageCache } from "@auth0/auth0-react";

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
  options.headers = options.headers || {} as Record<string, any>;
  options.headers["Authorization"] = `Bearer ${token}`;
  console.log("headers:");
  console.log(options.headers);
  return options;
}

export function useAccessToken(apiCaller: CallableFunction) {
  const { getAccessTokenSilently } = useAuth0();
  return async (...args: any) => {
    const token: string = await getAccessTokenSilently();
    console.log(`Retrieved access token: ${token}`);
    // Store refresh token
    const refresh_token = new LocalStorageCache;
    const key = refresh_token.allKeys().find((key: string | string[]) => key.includes('auth0spa')) || "";
    const refresh_token_value = refresh_token.get(key);
    // @ts-ignore
    const finalRefreshToken = refresh_token_value?.body?.refresh_token
    localStorage.setItem('refresh_token', finalRefreshToken);
    return apiCaller(...args, token);
  }
}

