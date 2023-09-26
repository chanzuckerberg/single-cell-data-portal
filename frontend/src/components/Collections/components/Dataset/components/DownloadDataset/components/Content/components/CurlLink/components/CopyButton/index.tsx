import copy from "clipboard-copy";
import { useEffect, useRef, useState } from "react";
import { CopyButton as Button } from "./style";
import { Fade } from "@mui/material";
import { ANIMATION, ANIMATION_STEP } from "./constants";

interface Props {
  curl: string;
  handleAnalytics: () => void;
}

export default function CopyButton({
  curl,
  handleAnalytics,
}: Props): JSX.Element {
  const [animationStep, setAnimationStep] = useState<ANIMATION_STEP>(
    ANIMATION_STEP.IDLE
  );
  const timeoutRef = useRef<NodeJS.Timer>();
  const animation = ANIMATION[animationStep];

  // Copy to clipboard, handle analytics, and initiate the copy animation.
  const handleCopyClick = () => {
    if (animation.isCopying) {
      return; // Copying is in progress.
    }
    copy(curl);
    handleAnalytics();
    setAnimationStep(incrementAnimationState);
  };

  // Callback fired after the "entered" status is applied.
  // Executes twice while copy animation is in progress: "COPIED_ENTER" and "COPY_ENTER".
  // Increments the animation step.
  const onEntered = () => {
    timeoutRef.current = setTimeout(() => {
      // Executes the next animation progression, after duration of the current animation is complete.
      setAnimationStep(incrementAnimationState);
    }, animation.duration);
  };

  // Callback fired after the "exited" status is applied.
  // Executes twice while copy animation is in progress: "COPY_EXIT" and "COPIED_EXIT".
  // Increments the animation step.
  const onExited = () => {
    setAnimationStep(incrementAnimationState);
  };

  // Clears timeout when unmounting.
  useEffect(() => {
    return () => {
      clearTimeout(timeoutRef.current);
    };
  }, []);

  return (
    <Button
      isAllCaps={false}
      onClick={handleCopyClick}
      sdsStyle="minimal"
      sdsType="primary"
    >
      <Fade
        appear={false}
        easing={animation.easing}
        in={animation.in}
        onEntered={onEntered}
        onExited={onExited}
        timeout={animation.timeout}
      >
        <span>{getButtonText(animationStep)}</span>
      </Fade>
    </Button>
  );
}

/**
 * Returns the button text to display.
 * @param step - Current animation step.
 * @returns button text.
 */
function getButtonText(step: number): string {
  if (
    step === ANIMATION_STEP.COPIED_ENTER ||
    step === ANIMATION_STEP.COPIED_EXIT
  ) {
    return "Copied";
  }
  return "Copy";
}

/**
 * Increments the animation step.
 * @param step - Current animation step.
 * @returns next animation step.
 */
function incrementAnimationState(step: number): number {
  if (step === ANIMATION_STEP.COPY_ENTER) {
    return ANIMATION_STEP.IDLE;
  }
  return step + 1;
}
