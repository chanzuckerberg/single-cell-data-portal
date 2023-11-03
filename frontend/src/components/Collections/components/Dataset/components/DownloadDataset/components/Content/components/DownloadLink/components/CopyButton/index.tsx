import copy from "clipboard-copy";
import { useEffect, useRef, useState } from "react";
import { CopyButton as Button } from "./style";
import { Fade } from "@mui/material";
import { ANIMATION, ANIMATION_STEP } from "./constants";

interface Props {
  downloadLink: string;
  handleAnalytics: () => void;
}

export default function CopyButton({
  downloadLink,
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
    copy(downloadLink);
    handleAnalytics();
    setAnimationStep(incrementAnimationState);
  };

  // Callback fired after the "exited" or "entered" status is applied.
  // Executes while copy animation is in progress for "COPY_EXIT", "COPIED_ENTER" and "COPIED_EXIT", "COPY_ENTER".
  // Increments the animation step.
  const onUpdateAnimationStep = () => {
    timeoutRef.current = setTimeout(() => {
      // Executes the next animation progression, after duration of the current animation is complete.
      setAnimationStep(incrementAnimationState);
    }, animation.duration);
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
        onEntered={onUpdateAnimationStep}
        onExited={onUpdateAnimationStep}
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
