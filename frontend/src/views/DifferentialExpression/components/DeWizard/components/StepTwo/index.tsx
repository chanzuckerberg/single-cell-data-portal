// import React, { useContext } from "react";
// import { DispatchContext, StateContext } from "../../common/store";
interface Props {
  setStep: (step: number) => void;
}
export default function StepTwo({ setStep }: Props): JSX.Element {
  // const state = useContext(StateContext);
  // const dispatch = useContext(DispatchContext);
  const handleGoNext = () => {
    setStep(3);
  };
  const handleGoBack = () => {
    setStep(1);
  };
  return <div>Hello world 2!</div>;
}
