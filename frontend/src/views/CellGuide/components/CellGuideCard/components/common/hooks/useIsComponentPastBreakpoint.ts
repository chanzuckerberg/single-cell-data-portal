import { useState, useLayoutEffect } from "react";

export function useIsComponentPastBreakpointWidth(breakpoint: number) {
  const [containerRef, setContainerRef] = useState<HTMLDivElement | null>(null);
  const [isPastBreakpoint, setIsPastBreakpoint] = useState(false);

  useLayoutEffect(() => {
    const handleResize = () => {
      if (containerRef) {
        setIsPastBreakpoint(containerRef.offsetWidth < breakpoint);
      }
    };

    const resizeObserver = new ResizeObserver(handleResize);
    if (containerRef) {
      resizeObserver.observe(containerRef);
    }

    return () => {
      resizeObserver.disconnect();
    };
  }, [breakpoint, containerRef]);

  return { isPastBreakpoint, containerRef: setContainerRef };
}

export function useIsComponentPastBreakpointHeight(breakpoint: number) {
  const [containerRef, setContainerRef] = useState<HTMLDivElement | null>(null);
  const [isPastBreakpoint, setIsPastBreakpoint] = useState(false);

  useLayoutEffect(() => {
    const handleResize = () => {
      if (containerRef) {
        setIsPastBreakpoint(containerRef.offsetHeight >= breakpoint);
      }
    };

    const resizeObserver = new ResizeObserver(handleResize);
    if (containerRef) {
      resizeObserver.observe(containerRef);
    }

    return () => {
      resizeObserver.disconnect();
    };
  }, [breakpoint, containerRef]);

  return { isPastBreakpoint, containerRef: setContainerRef };
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
