import { useReducer } from "react";
import {
  DispatchContext,
  INITIAL_STATE,
  reducer,
  StateContext,
} from "./common/store";
import LandingPage from "./components/LandingPage";

export default function CellCards(): JSX.Element {
  const [state, dispatch] = useReducer(reducer, INITIAL_STATE);

  return (
    <DispatchContext.Provider value={dispatch}>
      <StateContext.Provider value={state}>
        <LandingPage />
      </StateContext.Provider>
    </DispatchContext.Provider>
  );
}
