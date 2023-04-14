// import React, { useContext } from "react";
// import { DispatchContext, StateContext } from "../../common/store";

interface Props {
  step: number;
}

export default function DeWizard({ step }: Props): JSX.Element {
  // const state = useContext(StateContext);
  // const dispatch = useContext(DispatchContext);
  return (
    <div style={{ margin: "auto", fontSize: "50px" }}>
      We are currently on step {step}!
    </div>
  );
}
