import { createContext, Dispatch } from "react";
import { INITIAL_STATE, PayloadAction, State } from "./reducer";

export { reducer } from "./reducer";
export type { State };
export { INITIAL_STATE };
export const DispatchContext = createContext<Dispatch<
  PayloadAction<unknown>
> | null>(null);
export const StateContext = createContext<State>(INITIAL_STATE);
