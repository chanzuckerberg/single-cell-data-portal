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
