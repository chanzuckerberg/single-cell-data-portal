import { Intent } from "@blueprintjs/core";
import React, { useCallback, useEffect, useRef, useState } from "react";
import {
  CollectionDescription as Description,
  DescriptionText,
  DESCRIPTION_LINE_HEIGHT_PX,
  MAX_LINE_COUNT,
} from "src/components/Collection/components/CollectionDescription/style";
import { StyledPrimaryAnchorButton } from "src/components/common/Button/common/style";

enum EllipsisMode {
  "NONE" = "NONE",
  "OFF" = "OFF",
  "ON" = "ON",
}

interface Props {
  description: string;
}

export default function CollectionDescription({
  description,
}: Props): JSX.Element {
  const descriptionRef = useRef<HTMLParagraphElement>(null);
  const [ellipsisMode, setEllipsisMode] = useState<EllipsisMode>(
    EllipsisMode.NONE
  );
  const isEllipsis = ellipsisMode === EllipsisMode.ON;

  /**
   * Updates the state of ellipsis mode by toggling between "ON" and "OFF" values.
   */
  const onToggleMode = (): void => {
    setEllipsisMode((currentMode) => {
      if (currentMode === EllipsisMode.ON) {
        return EllipsisMode.OFF;
      }
      return EllipsisMode.ON;
    });
  };

  /**
   * Updates ellipsis mode state.
   */
  const onUpdateMode = useCallback((): void => {
    if (descriptionRef.current) {
      const el = descriptionRef.current;
      setEllipsisMode((currentMode) => getEllipsisMode(el, currentMode));
    }
  }, []);

  useEffect(() => {
    onUpdateMode(); /* Initialize ellipsis mode. */
    window.addEventListener("resize", onUpdateMode);
    return () => {
      document.removeEventListener("resize", onUpdateMode);
    };
  }, [onUpdateMode]);

  return (
    <Description data-test-id="collection-description">
      <DescriptionText isEllipsis={isEllipsis} ref={descriptionRef}>
        {description}
      </DescriptionText>
      {isModeActivated(ellipsisMode) && (
        <StyledPrimaryAnchorButton
          intent={Intent.PRIMARY}
          minimal
          onClick={onToggleMode}
          text={getModeText(ellipsisMode)}
        />
      )}
    </Description>
  );
}

/**
 * Returns one of the following ellipsis mode for the description paragraph:
 * - "ON" - paragraph is ellipsis view with the option to toggle view, or
 * - "OFF" - paragraph is in full view with the option to toggle view, or
 * - "NONE" - paragraph is in full view with no option to toggle view as it does not require ellipsis.
 * @param el - description paragraph element.
 * @param currentMode - current ellipsis mode.
 * @returns EllipsisMode
 */
function getEllipsisMode(
  el: HTMLParagraphElement,
  currentMode: EllipsisMode
): EllipsisMode {
  /* Grab the measurement of element height, including hidden content due to overflow. */
  const elScrollHeight = el.scrollHeight;
  /* Calculate the number of lines for the content fully rendered. */
  const elLineCount = elScrollHeight / DESCRIPTION_LINE_HEIGHT_PX;

  if (elLineCount <= MAX_LINE_COUNT) {
    /* Element does not have hidden content and line count is within allowable limit. */
    /* Mode is "NONE" - ellipsis mode not required. */
    return EllipsisMode.NONE;
  }

  if (currentMode === EllipsisMode.NONE) {
    /* Line count exceeds allowable limit, and ellipsis mode is "NONE". */
    /* Change mode to "ON" - ellipsis mode required. */
    return EllipsisMode.ON;
  }

  return currentMode;
}

/**
 * Returns applicable button text for display corresponding with current mode.
 * When the current mode is "ON" the return value will be "Show More", otherwise the return value is "Show Less".
 * @param currentMode - current ellipsis mode.
 * @returns string for display as button text.
 */
function getModeText(currentMode: EllipsisMode): string {
  if (currentMode === EllipsisMode.ON) {
    return "Show More";
  }

  return "Show Less";
}

/**
 * Ellipsis mode is "activated" when description text is long enough to be truncated.
 * i.e. the current mode will not equal "NONE".
 * @param currentMode - current ellipsis mode.
 * @returns true when mode is not "NONE".
 */
function isModeActivated(currentMode: EllipsisMode): boolean {
  return currentMode !== EllipsisMode.NONE;
}
