import { useReducer } from "react";
import {
  DispatchContext,
  INITIAL_STATE,
  reducer,
  StateContext,
} from "./common/store";
import Main from "./components/Main";

export default function CelCards(): JSX.Element {
  const [state, dispatch] = useReducer(reducer, INITIAL_STATE);

  return (
    <DispatchContext.Provider value={dispatch}>
      <StateContext.Provider value={state}>
        <Main />
      </StateContext.Provider>
    </DispatchContext.Provider>
  );
}
