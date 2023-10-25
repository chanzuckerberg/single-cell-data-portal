import { useReducer } from "react";
import {
  DispatchContext,
  INITIAL_STATE,
  reducer,
  StateContext,
} from "../WheresMyGeneV2/common/store";
import Main from "./components/Main";

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
