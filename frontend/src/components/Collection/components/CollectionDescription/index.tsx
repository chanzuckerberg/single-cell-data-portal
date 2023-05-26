import React, { useEffect, useRef, useState } from "react";
import { useResizeObserver } from "src/common/hooks/useResizeObserver";
import {
  CollectionDescription as Description,
  DESCRIPTION_LINE_HEIGHT_PX,
  DescriptionText,
  MAX_LINE_COUNT,
} from "src/components/Collection/components/CollectionDescription/style";
import { Button } from "@czi-sds/components";

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
  const [ellipsisMode, setEllipsisMode] = useState<EllipsisMode | undefined>();
  const isEllipsis = ellipsisMode === EllipsisMode.ON;
  const descriptionRect = useResizeObserver(descriptionRef);
  const { scrollHeight: descriptionScrollHeight } = descriptionRect || {};

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
   * Ellipsis mode state is initialized to "NONE" only when the component mounts,
   * so we guarantee that `descriptionRef` gets the element handle for resize observer hook.
   */
  useEffect(() => {
    setEllipsisMode(EllipsisMode.NONE);
  }, []);

  /**
   * Ellipsis mode state updates with changes to the description paragraph scroll height.
   */
  useEffect(() => {
    setEllipsisMode((currentMode) =>
      getEllipsisMode(descriptionScrollHeight, currentMode)
    );
  }, [descriptionScrollHeight]);

  return (
    <Description data-testid="collection-description">
      <DescriptionText isEllipsis={isEllipsis} ref={descriptionRef}>
        {description}
      </DescriptionText>
      {isModeActivated(ellipsisMode) && (
        <Button
          isAllCaps={false}
          onClick={onToggleMode}
          sdsStyle="minimal"
          sdsType="primary"
        >
          {getModeText(ellipsisMode)}
        </Button>
      )}
    </Description>
  );
}

/**
 * Returns one of the following ellipsis mode for the description paragraph:
 * - "ON" - paragraph is ellipsis view with the option to toggle view, or
 * - "OFF" - paragraph is in full view with the option to toggle view, or
 * - "NONE" - paragraph is in full view with no option to toggle view as it does not require ellipsis.
 * @param elScrollHeight - description paragraph scroll height.
 * @param currentMode - current ellipsis mode for the description paragraph.
 * @returns Ellipsis mode for the description paragraph.
 */
function getEllipsisMode(
  elScrollHeight: number | undefined,
  currentMode: EllipsisMode = EllipsisMode.NONE
): EllipsisMode {
  if (!elScrollHeight) {
    return currentMode;
  }

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
function getModeText(currentMode?: EllipsisMode): string {
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
function isModeActivated(currentMode?: EllipsisMode): boolean {
  return !!currentMode && currentMode !== EllipsisMode.NONE;
}
