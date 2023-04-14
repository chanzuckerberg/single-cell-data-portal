// import React, { useContext } from "react";
// import { DispatchContext, StateContext } from "../../common/store";
interface Props {
  setStep: (step: number) => void;
}
export default function StepThree({ setStep }: Props): JSX.Element {
  // const state = useContext(StateContext);
  // const dispatch = useContext(DispatchContext);

  const handleGoBack = () => {
    setStep(2);
  };
  const handleStartOver = () => {
    setStep(1);
  };

  return <div>Hello world 3!</div>;
}
