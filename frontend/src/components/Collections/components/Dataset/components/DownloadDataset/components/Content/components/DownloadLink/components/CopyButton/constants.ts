import { FadeProps } from "@mui/material";

export interface Animation
  extends Pick<FadeProps, "easing" | "in" | "timeout"> {
  duration: number;
  isCopying: boolean;
}

const ANIMATION_300_EASE_OUT = {
  easing: "ease-out",
  timeout: 250,
};

const ANIMATION_200_EASE_IN = {
  easing: "ease-in",
  timeout: 150,
};

// Possible animation step values.
export enum ANIMATION_STEP {
  COPIED_ENTER = 2,
  COPIED_EXIT = 3,
  COPY_ENTER = 4,
  COPY_EXIT = 1,
  IDLE = 0,
}

// Animation properties for each step.
export const ANIMATION: Record<ANIMATION_STEP, Animation> = {
  [ANIMATION_STEP.IDLE]: {
    ...ANIMATION_300_EASE_OUT,
    in: true,
    isCopying: false,
    duration: 0,
  },
  [ANIMATION_STEP.COPY_EXIT]: {
    ...ANIMATION_300_EASE_OUT,
    in: false,
    isCopying: true,
    duration: 0,
  },
  [ANIMATION_STEP.COPIED_ENTER]: {
    ...ANIMATION_300_EASE_OUT,
    in: true,
    isCopying: true,
    duration: 1000, // Show "Copied" for 1 second.
  },
  [ANIMATION_STEP.COPIED_EXIT]: {
    ...ANIMATION_200_EASE_IN,
    in: false,
    isCopying: true,
    duration: 0,
  },
  [ANIMATION_STEP.COPY_ENTER]: {
    ...ANIMATION_200_EASE_IN,
    in: true,
    isCopying: true,
    duration: 0,
  },
};
