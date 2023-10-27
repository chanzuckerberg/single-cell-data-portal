import { useReducer } from "react";

import Main from "./components/Main";
import {
  DispatchContext,
  INITIAL_STATE,
  StateContext,
  reducer,
} from "./common/store";

export default function WheresMyGene(): JSX.Element {
  const [state, dispatch] = useReducer(reducer, INITIAL_STATE);

  return (
    <DispatchContext.Provider value={dispatch}>
      <StateContext.Provider value={state}>
        <Main />
      </StateContext.Provider>
    </DispatchContext.Provider>
  );
}
