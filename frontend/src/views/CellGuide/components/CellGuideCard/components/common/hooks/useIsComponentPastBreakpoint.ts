import { useState, useLayoutEffect } from "react";

type Condition = (containerRef: HTMLDivElement) => boolean;

function useIsComponentPastBreakpoint(condition: Condition) {
  const [containerRef, setContainerRef] = useState<HTMLDivElement | null>(null);
  const [isPastBreakpoint, setIsPastBreakpoint] = useState(false);

  useLayoutEffect(() => {
    const handleResize = () => {
      if (containerRef) {
        setIsPastBreakpoint(condition(containerRef));
      }
    };

    const resizeObserver = new ResizeObserver(handleResize);
    if (containerRef) {
      resizeObserver.observe(containerRef);
    }

    return () => {
      resizeObserver.disconnect();
    };
  }, [containerRef, condition]);

  return { isPastBreakpoint, containerRef: setContainerRef };
}

export function useIsComponentPastBreakpointWidth(breakpoint: number) {
  const condition = (containerRef: HTMLDivElement) =>
    containerRef.offsetWidth < breakpoint;
  return useIsComponentPastBreakpoint(condition);
}

export function useIsComponentPastBreakpointHeight(breakpoint: number) {
  const condition = (containerRef: HTMLDivElement) =>
    containerRef.offsetHeight >= breakpoint;
  return useIsComponentPastBreakpoint(condition);
}

export function useComponentWidth() {
  const [containerRef, setContainerRef] = useState<HTMLDivElement | null>(null);
  const [width, setWidth] = useState(0);

  useLayoutEffect(() => {
    const handleResize = () => {
      if (containerRef) {
        setWidth(containerRef.offsetWidth);
      }
    };

    const resizeObserver = new ResizeObserver(handleResize);
    if (containerRef) {
      resizeObserver.observe(containerRef);
    }

    return () => {
      resizeObserver.disconnect();
    };
  }, [containerRef]);

  return { width, containerRef: setContainerRef };
}
