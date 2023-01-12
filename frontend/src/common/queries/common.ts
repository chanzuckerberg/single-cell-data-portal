import { Auth0ContextInterface } from "@auth0/auth0-react";

export const DEFAULT_FETCH_OPTIONS: RequestInit = {
  credentials: "include",
};

export const DELETE_FETCH_OPTIONS: RequestInit = {
  credentials: "include",
  method: "DELETE",
};

export const CONTENT_TYPE_APPLICATION_JSON: object = {"Content-Type": "application/json"};

export function withAccessToken(
  getAccessTokenSilently: Auth0ContextInterface["getAccessTokenSilently"],
  apiCaller: CallableFunction
) {
  return async (...args: any) => {
    const token: string = await getAccessTokenSilently();
    console.log(`Retrieved access token: ${token}`);
    return apiCaller(...args, token);
  }
}

